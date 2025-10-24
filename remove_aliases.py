#!/usr/bin/env python3
"""
Script to clean and transform SQL queries in JOB_COMPLETE.sql
Performs the following transformations:
1. Remove table aliases: "SELECT alias.column FROM table AS alias" -> "SELECT table.column FROM table"
2. Remove column aliases: "SELECT column AS alias" -> "SELECT column"
3. Convert IN clauses: "column IN (a, b, c)" -> "column = a OR column = b OR column = c"
4. Convert BETWEEN clauses: "column BETWEEN x AND y" -> "column >= x AND column <= y"
5. Remove ORDER BY clauses and ASC/DESC operators
6. Optionally remove GROUP BY clauses and aggregation functions together (use --remove-group-by flag)

Usage:
    python remove_aliases.py                    # Standard processing (keeps aggregation and GROUP BY)
    python remove_aliases.py --remove-group-by  # Remove both aggregation functions and GROUP BY clauses
"""

import re
import sys
import os

def extract_table_aliases(query):
    """Extract table aliases from a SQL query's FROM clause."""
    aliases = {}
    
    alias_pattern = r'(\w+)\s+AS\s+(\w+)'
    matches = re.findall(alias_pattern, query, re.IGNORECASE)
    
    for table_name, alias_name in matches:
        aliases[alias_name] = table_name
    
    return aliases

def remove_aliases_from_query(query, remove_group_by=False):
    """Remove aliases from a single SQL query."""
    aliases = extract_table_aliases(query)

    result = query
    
    for alias, table in aliases.items():
        pattern = rf'\b{re.escape(alias)}\.'
        replacement = f'{table}.'
        result = re.sub(pattern, replacement, result)
    
    for alias, table in aliases.items():
        pattern = rf'\b{re.escape(table)}\s+AS\s+{re.escape(alias)}\b'
        replacement = table
        result = re.sub(pattern, replacement, result, flags=re.IGNORECASE)
    
    result = re.sub(r'\s+AS\s+\w+', '', result, flags=re.IGNORECASE)
    
    result = convert_in_clauses(result)
    
    result = convert_between_clauses(result)
    
    if remove_group_by:
        result = remove_aggregation_functions(result)
        result = remove_group_by_clauses(result)
    
    result = remove_order_by_clauses(result)
    
    result = remove_null_and_not_conditions(result)
    
    return result

def convert_in_clauses(query):
    """Convert IN clauses to OR conditions."""
    pattern = r'(\w+(?:\.\w+)?)\s+IN\s+\('
    
    def replace_in_clause(match):
        column = match.group(1)
        start_pos = match.end()
        
        paren_count = 1
        in_quotes = False
        quote_char = None
        i = start_pos
        
        while i < len(query) and paren_count > 0:
            char = query[i]
            
            if char in ('"', "'") and not in_quotes:
                in_quotes = True
                quote_char = char
            elif char == quote_char and in_quotes:
                if i + 1 < len(query) and query[i + 1] == quote_char:
                    i += 1 
                else:
                    in_quotes = False
                    quote_char = None
            elif not in_quotes:
                if char == '(':
                    paren_count += 1
                elif char == ')':
                    paren_count -= 1
            
            i += 1
        
        if paren_count > 0:
            return match.group(0)
        
        values_str = query[start_pos:i-1].strip()
        
        values = []
        current_value = ""
        in_quotes = False
        quote_char = None
        paren_depth = 0
        
        j = 0
        while j < len(values_str):
            char = values_str[j]
            
            if char in ('"', "'") and not in_quotes:
                in_quotes = True
                quote_char = char
                current_value += char
            elif char == quote_char and in_quotes:
                if j + 1 < len(values_str) and values_str[j + 1] == quote_char:
                    current_value += char + char
                    j += 1
                else:
                    in_quotes = False
                    quote_char = None
                    current_value += char
            elif char == '(' and not in_quotes:
                paren_depth += 1
                current_value += char
            elif char == ')' and not in_quotes:
                paren_depth -= 1
                current_value += char
            elif char == ',' and not in_quotes and paren_depth == 0:
                if current_value.strip():
                    values.append(current_value.strip())
                current_value = ""
            else:
                current_value += char
            
            j += 1
        
        if current_value.strip():
            values.append(current_value.strip())
        
        or_conditions = []
        for value in values:
            value = value.strip()
            if value: 
                or_conditions.append(f'{column} = {value}')
        
        return '(' + ' OR '.join(or_conditions) + ')'
    
    while True:
        match = re.search(pattern, query, flags=re.IGNORECASE)
        if not match:
            break
        
        replacement = replace_in_clause(match)
        start = match.start()
        
        paren_count = 1
        in_quotes = False
        quote_char = None
        i = match.end()
        
        while i < len(query) and paren_count > 0:
            char = query[i]
            
            if char in ('"', "'") and not in_quotes:
                in_quotes = True
                quote_char = char
            elif char == quote_char and in_quotes:
                if i + 1 < len(query) and query[i + 1] == quote_char:
                    i += 1 
                else:
                    in_quotes = False
                    quote_char = None
            elif not in_quotes:
                if char == '(':
                    paren_count += 1
                elif char == ')':
                    paren_count -= 1
            
            i += 1
        
        query = query[:start] + replacement + query[i:]
    
    return query

def convert_between_clauses(query):
    """Convert BETWEEN clauses to >= AND <= conditions."""
    pattern = r'(\w+(?:\.\w+)?)\s+BETWEEN\s+([^\s]+)\s+AND\s+([^\s]+)'
    
    def replace_between_clause(match):
        column = match.group(1)
        value1 = match.group(2)
        value2 = match.group(3)
        
        return f'{column} >= {value1} AND {column} <= {value2}'
    
    result = re.sub(pattern, replace_between_clause, query, flags=re.IGNORECASE)
    return result

def remove_aggregation_functions(query):
    """Remove aggregation functions (MIN, MAX, SUM) and keep only the content within parentheses."""
    pattern = r'\b(?:MIN|MAX|SUM)\s*\(([^)]+)\)'
    
    def replace_aggregation(match):
        return match.group(1)
    
    result = re.sub(pattern, replace_aggregation, query, flags=re.IGNORECASE)
    return result

def remove_order_by_clauses(query):
    """Remove ORDER BY clauses and ASC/DESC operators from the query."""
    # Remove entire ORDER BY clauses (including everything until semicolon or end of query)
    # Pattern matches ORDER BY followed by any content until semicolon or end of string
    pattern = r'\s+ORDER\s+BY\s+[^;]*(?=;|$)'
    result = re.sub(pattern, '', query, flags=re.IGNORECASE)
    
    # Also remove standalone ASC/DESC keywords that might remain
    result = re.sub(r'\s+(?:ASC|DESC)\b', '', result, flags=re.IGNORECASE)
    
    return result

def remove_group_by_clauses(query):
    """Remove GROUP BY clauses from the query."""
    # Remove GROUP BY clauses (everything from GROUP BY until ORDER BY, semicolon, or end of query)
    pattern = r'\s+GROUP\s+BY\s+[^;]*?(?=\s+ORDER\s+BY|;|$)'
    result = re.sub(pattern, '', query, flags=re.IGNORECASE)
    
    return result

def remove_null_and_not_conditions(query):
    """Remove IS NULL, IS NOT NULL, and other NOT conditions from the query."""
    lines = query.split('\n')
    processed_lines = []
    
    for line in lines:
        original_line = line
        
        line = re.sub(r'\s+AND\s+\w+(?:\.\w+)?\s+IS\s+NULL', '', line, flags=re.IGNORECASE)
        line = re.sub(r'\w+(?:\.\w+)?\s+IS\s+NULL\s+AND\s+', '', line, flags=re.IGNORECASE)
        line = re.sub(r'\w+(?:\.\w+)?\s+IS\s+NULL', '', line, flags=re.IGNORECASE)
        
        line = re.sub(r'\s+AND\s+\w+(?:\.\w+)?\s+IS\s+NOT\s+NULL', '', line, flags=re.IGNORECASE)
        line = re.sub(r'\w+(?:\.\w+)?\s+IS\s+NOT\s+NULL\s+AND\s+', '', line, flags=re.IGNORECASE)
        line = re.sub(r'\w+(?:\.\w+)?\s+IS\s+NOT\s+NULL', '', line, flags=re.IGNORECASE)
        
        line = re.sub(r'\s+AND\s+\w+(?:\.\w+)?\s+NOT\s+LIKE\s+["\'][^"\']*["\']', '', line, flags=re.IGNORECASE)
        line = re.sub(r'\w+(?:\.\w+)?\s+NOT\s+LIKE\s+["\'][^"\']*["\']\s+AND\s+', '', line, flags=re.IGNORECASE)
        line = re.sub(r'\w+(?:\.\w+)?\s+NOT\s+LIKE\s+["\'][^"\']*["\']', '', line, flags=re.IGNORECASE)
        
        line = re.sub(r'\s+AND\s+\w+(?:\.\w+)?\s+!=\s*["\'][^"\']*["\']', '', line, flags=re.IGNORECASE)
        line = re.sub(r'\w+(?:\.\w+)?\s+!=\s*["\'][^"\']*["\']\s+AND\s+', '', line, flags=re.IGNORECASE)
        line = re.sub(r'\w+(?:\.\w+)?\s+!=\s*["\'][^"\']*["\']', '', line, flags=re.IGNORECASE)
        
        line = re.sub(r'\s+AND\s+AND\s+', ' AND ', line, flags=re.IGNORECASE)
        line = re.sub(r'\s+OR\s+OR\s+', ' OR ', line, flags=re.IGNORECASE)
        
        line = re.sub(r'\bWHERE\s+AND\s+', 'WHERE ', line, flags=re.IGNORECASE)
        
        line = re.sub(r'\s+(?:AND|OR)\s*;?\s*$', ';' if original_line.strip().endswith(';') else '', line, flags=re.IGNORECASE)
        
        line = re.sub(r'\bWHERE\s*;', ';', line, flags=re.IGNORECASE)
        
        processed_lines.append(line)
    
    return '\n'.join(processed_lines)

def process_sql_file(input_file, output_file, remove_group_by=False):
    """Process the entire SQL file to remove aliases."""
    with open(input_file, 'r') as f:
        content = f.read()
    
    lines = content.split('\n')
    processed_lines = []
    
    current_query = ""
    in_query = False
    
    for line in lines:
        stripped = line.strip()
        
        if stripped.upper().startswith('SELECT'):
            if current_query:
                processed_query = remove_aliases_from_query(current_query, remove_group_by)
                processed_lines.append(processed_query)
                current_query = ""
            current_query = line
            in_query = True
        elif in_query and stripped.endswith(';'):
            current_query += "\n" + line
            processed_query = remove_aliases_from_query(current_query, remove_group_by)
            processed_lines.append(processed_query)
            current_query = ""
            in_query = False
        elif in_query:
            current_query += "\n" + line
        else:
            processed_lines.append(line)
    
    if current_query:
        processed_query = remove_aliases_from_query(current_query, remove_group_by)
        processed_lines.append(processed_query)
    
    with open(output_file, 'w') as f:
        f.write('\n'.join(processed_lines))

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Clean and transform SQL queries')
    parser.add_argument('--remove-group-by', action='store_true', 
                       help='Remove GROUP BY clauses from queries')
    args = parser.parse_args()
    
    input_files = ['/Users/oliver/TU_BERLIN/MASTER/mutable/benchmark/job/JOB_COMPLETE.sql', '/Users/oliver/TU_BERLIN/MASTER/mutable/benchmark/ssb/ssb_queries.sql']

    for input_file in input_files:
        try:
            if not input_file:
                continue
            base, ext = os.path.splitext(input_file)
            suffix = "_filtered_no_group" if args.remove_group_by else "_filtered"
            output_file = f"{base}{suffix}{ext}"
            process_sql_file(input_file, output_file, args.remove_group_by)
            print(f"Successfully processed {input_file}")
            print(f"Output written to {output_file}")
            
            print("\nSample transformations:")
            with open(input_file, 'r') as f:
                original_lines = f.readlines()
            with open(output_file, 'r') as f:
                processed_lines = f.readlines()
                
            for i, line in enumerate(original_lines):
                if line.strip().upper().startswith('SELECT'):
                    print(f"\nOriginal query {i+1}:")
                    query_lines = []
                    for j in range(i, min(i+5, len(original_lines))):
                        query_lines.append(original_lines[j].rstrip())
                        if original_lines[j].strip().endswith(';'):
                            break
                    print(''.join(line + '\n' for line in query_lines))
                    
                    print(f"Transformed query {i+1}:")
                    query_lines = []
                    for j in range(i, min(i+5, len(processed_lines))):
                        query_lines.append(processed_lines[j].rstrip())
                        if processed_lines[j].strip().endswith(';'):
                            break
                    print(''.join(line + '\n' for line in query_lines))
                    break
                    
        except Exception as e:
            print(f"Error processing file: {e}")
            sys.exit(1)

if __name__ == "__main__":
    main()