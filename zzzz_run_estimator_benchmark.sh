#!/bin/bash
# filepath: /Users/oliver/TU_BERLIN/MASTER/mutable/zzzz_run_estimator_benchmark.sh

# Create output directory if it doesn't exist
mkdir -p z_output_data_

# Set script variables
SHELL_PATH="./build/debug_shared/bin/shell"
SQL_FILE="benchmark/job-light/job-light_1.sql"
LOG_FILE="z_output_data_/benchmark_run.log"

# Remove existing log file
rm -f "$LOG_FILE"

echo "Starting benchmark run at $(date)" | tee -a "$LOG_FILE"
echo "Results will be saved in z_output_data_/" | tee -a "$LOG_FILE"

# Function to run a benchmark with specific estimator
run_benchmark() {
    local estimator=$1
    local csv_file=$2
    local extra_args=$3
    
    # Remove existing CSV file
    rm -f "$csv_file"
    
    echo "----------------------------------------" | tee -a "$LOG_FILE"
    echo "Running benchmark with $estimator estimator" | tee -a "$LOG_FILE"
    echo "Output CSV: $csv_file" | tee -a "$LOG_FILE"
    echo "Started at: $(date)" | tee -a "$LOG_FILE"
    
    # Run the benchmark
    $SHELL_PATH $extra_args \
      --cardinality-csv "$csv_file" \
      --plan-enumerator GOO \
      --cardinality-estimator "$estimator" \
      --backend Interpreter \
      "$SQL_FILE" 2>&1 | tee -a "$LOG_FILE"
    
    echo "Finished at: $(date)" | tee -a "$LOG_FILE"
    echo "----------------------------------------" | tee -a "$LOG_FILE"
    echo "" | tee -a "$LOG_FILE"
}

# Run CartesianProduct estimator
run_benchmark "CartesianProduct" "z_output_data_/CartesianProduct.csv" "--no-statistics"

# Run SelectivityBased estimator
run_benchmark "Selectivitybased" "z_output_data_/Selectivity.csv" ""

# Run Histogram estimator  
run_benchmark "Histogram" "z_output_data_/Histogram.csv" ""

echo "All benchmarks completed at $(date)" | tee -a "$LOG_FILE"
echo "Log saved to $LOG_FILE"