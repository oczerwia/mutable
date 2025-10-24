#!/usr/bin/env python3
# filepath: /Users/oliver/TU_BERLIN/MASTER/mutable/run_parameter_sweep_fixed.py

import sys
import os
sys.path.append('benchmark')

from database_connectors.mutable import MutableConnector
from benchmark_utils import parse_schema_from_yaml, load_queries_from_yaml
import csv
import time
from pathlib import Path
import yaml

# Configuration
BASE_OUTPUT_DIR = "z_output_data_"
DATA_DIR = f"{BASE_OUTPUT_DIR}/data"
LOG_DIR = f"{BASE_OUTPUT_DIR}/logs"

# Create directories
Path(DATA_DIR).mkdir(parents=True, exist_ok=True)
Path(LOG_DIR).mkdir(parents=True, exist_ok=True)

# Parameters (same as your shell script)
QUERY_FILES = [
    "benchmark/job-light/ALL_QUERIES.sql",
    "benchmark/job/JOB_COMPLETE_filtered.sql"
]

PLAN_ENUMERATORS = ["RangeGOO", "DPsizeOptRange"]
COLLAPSE_FUNCTIONS = ["Mean", "Upper", "Uncertainty"]
TOPK_SIZES = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768]
LEARN_MODES = ["no-learn", "learn"]
UNCERTAINTY_IMPACTS = [0.1, 0.2, 0.3, 0.4]
REPETITIONS = 5

def parse_queries_from_sql_file(filepath):
    """Parse individual queries from SQL file (like your shell script does)"""
    queries = []
    current_query_lines = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Skip comments and empty lines
            if not line or line.startswith('--'):
                if current_query_lines and line.startswith('-- Query'):
                    # Save previous query
                    if current_query_lines:
                        query_text = ' '.join(current_query_lines).strip()
                        if query_text and not query_text.upper().startswith(('CREATE', 'USE', 'IMPORT')):
                            queries.append(query_text)
                    current_query_lines = []
                continue
            
            # Skip database setup commands
            if line.upper().startswith(('CREATE DATABASE', 'USE ', 'IMPORT INTO')):
                continue
                
            current_query_lines.append(line)
    
    # Add last query
    if current_query_lines:
        query_text = ' '.join(current_query_lines).strip()
        if query_text and not query_text.upper().startswith(('CREATE', 'USE', 'IMPORT')):
            queries.append(query_text)
    
    return queries

def run_single_experiment(config):
    """Run a single experiment configuration"""
    query_file, plan_enumerator, collapse_func, topk_size, learn_mode, uncertainty_impact, rep = config
    
    # Determine dataset name
    if "job-light" in query_file:
        dataset = "job-light"
    elif "job" in query_file:
        dataset = "job-complete"
    else:
        dataset = "unknown"
    
    # Create experiment name
    if uncertainty_impact is not None:
        exp_name = f"{dataset}_{plan_enumerator}_{collapse_func}_topk{topk_size}_{learn_mode}_impact{uncertainty_impact}_run{rep}"
    else:
        exp_name = f"{dataset}_{plan_enumerator}_{collapse_func}_topk{topk_size}_{learn_mode}_run{rep}"
    
    print(f"Running experiment: {exp_name}")
    
    # Build command arguments (same as your shell script)
    mutable_args = [
        "build/debug_shared/bin/shell",
        "--data-layout", "Row",
        "--plan-enumerator", plan_enumerator,
        "--cardinality-estimator", "ExperimentalRange",
        "--backend", "Interpreter",
        "--cardinality-csv", f"{DATA_DIR}/{exp_name}.csv",
        "--collapse-function", collapse_func,
        "--topk-size", str(topk_size)
    ]
    
    if uncertainty_impact is not None:
        mutable_args.extend(["--uncertainty-impact", str(uncertainty_impact)])
    
    if learn_mode == "no-learn":
        mutable_args.append("--no-learn-cardinalities")
    
    mutable_args.append(query_file)
    
    try:
        # Create connector with specific parameters for this experiment
        connector = MutableConnector(
            database_name="job",
            schema_name="job", 
            timeout=1800,  # 30 minute timeout
            mutable_args=mutable_args  # Pass the full argument list
        )
        
        start_time = time.time()
        
        # The connector will handle:
        # 1. Starting mutable with the right parameters
        # 2. Loading the schema and data
        # 3. Executing all queries in the file
        result = connector.run_benchmark()
        
        end_time = time.time()
        
        # Log results
        with open(f"{LOG_DIR}/{exp_name}.log", 'w') as f:
            f.write(f"Experiment: {exp_name}\n")
            f.write(f"Total time: {end_time - start_time:.2f}s\n")
            f.write(f"Result: {result}\n")
        
        print(f"SUCCESS: {exp_name} ({end_time - start_time:.1f}s)")
        
    except Exception as e:
        print(f"ERROR: {exp_name} - {str(e)}")
        with open(f"{LOG_DIR}/{exp_name}_error.log", 'w') as f:
            f.write(f"Error in experiment: {exp_name}\n")
            f.write(f"Error: {str(e)}\n")
    
    finally:
        # Ensure cleanup
        if 'connector' in locals():
            try:
                connector.cleanup()
            except:
                pass

def main():
    # Generate all experiment configurations (same as your shell script)
    experiments = []
    
    for query_file in QUERY_FILES:
        for plan_enumerator in PLAN_ENUMERATORS:
            for collapse_func in COLLAPSE_FUNCTIONS:
                for topk_size in TOPK_SIZES:
                    for learn_mode in LEARN_MODES:
                        for rep in range(1, REPETITIONS + 1):
                            if collapse_func == "Uncertainty":
                                for uncertainty_impact in UNCERTAINTY_IMPACTS:
                                    experiments.append((
                                        query_file, plan_enumerator, collapse_func,
                                        topk_size, learn_mode, uncertainty_impact, rep
                                    ))
                            else:
                                experiments.append((
                                    query_file, plan_enumerator, collapse_func,
                                    topk_size, learn_mode, None, rep
                                ))
    
    print(f"Total experiments to run: {len(experiments)}")
    
    # Calculate expected counts (same as your shell script)
    uncertainty_count = len(QUERY_FILES) * len(PLAN_ENUMERATORS) * 1 * len(TOPK_SIZES) * len(LEARN_MODES) * len(UNCERTAINTY_IMPACTS) * REPETITIONS
    other_count = len(QUERY_FILES) * len(PLAN_ENUMERATORS) * 2 * len(TOPK_SIZES) * len(LEARN_MODES) * REPETITIONS
    print(f"Expected: {uncertainty_count} Uncertainty + {other_count} Mean/Upper = {uncertainty_count + other_count} total")
    
    # Run experiments sequentially (each needs its own mutable process)
    for i, config in enumerate(experiments):
        print(f"Progress: {i+1}/{len(experiments)}")
        run_single_experiment(config)
        
        # Short delay between experiments to ensure clean process shutdown
        time.sleep(0.5)
    
    print("All experiments completed!")
    print(f"Results stored in: {BASE_OUTPUT_DIR}")

if __name__ == "__main__":
    main()