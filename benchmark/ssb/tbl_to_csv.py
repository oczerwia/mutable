#!/usr/bin/env python3
"""
Convert .tbl files to .csv files in a directory.
Usage: python tbl_to_csv.py <directory>
"""

import csv
import sys
from pathlib import Path

def main():
    if len(sys.argv) != 2:
        print("Usage: python tbl_to_csv.py <directory>")
        sys.exit(1)
    
    directory = Path(sys.argv[1])
    if not directory.exists() or not directory.is_dir():
        print(f"Error: {directory} does not exist or is not a directory")
        sys.exit(1)
    
    for tbl_file in directory.glob("*.tbl"):
        csv_file = tbl_file.with_suffix(".csv")
        
        with tbl_file.open("r") as inf, csv_file.open("w", newline="") as outf:
            writer = csv.writer(outf)
            for line in inf:
                line = line.rstrip("\n\r")
                if line.endswith("|"):
                    line = line[:-1]
                row = line.split("|")
                writer.writerow(row)
        
        print(f"Converted: {tbl_file} -> {csv_file}")

if __name__ == "__main__":
    main()