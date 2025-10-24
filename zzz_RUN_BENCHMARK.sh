#!/bin/bash

BASE_OUTPUT_DIR="z_output_data_"
DATA_DIR="${BASE_OUTPUT_DIR}/data"
LOG_DIR="${BASE_OUTPUT_DIR}/logs"

if [ -d "$BASE_OUTPUT_DIR" ]; then
  rm -rf "$BASE_OUTPUT_DIR"
fi
mkdir -p "$DATA_DIR"
mkdir -p "$LOG_DIR"

PROGRAM="${PWD}/build/debug_shared/bin/shell"

# Query files
QUERY_FILES=(
    "benchmark/job-light/ALL_QUERIES.sql"
    "benchmark/job/JOB_COMPLETE_filtered.sql"
)

# Plan enumerators
PLAN_ENUMERATORS=("RangeGOO" "DPsizeOptRange")

COLLAPSE_FUNCTIONS=("Mean" "Upper" "Uncertainty")

TOPK_SIZES=(2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384 32768)

LEARN_MODES=("no-learn" "learn")

# Uncertainty impact values (only used with Uncertainty collapse function)
UNCERTAINTY_IMPACTS=(0.1 0.2 0.3 0.4)

REPETITIONS=5

MAX_JOBS=6

echo "Starting comprehensive ExperimentalRange benchmark..."
echo "Data will be saved to: $DATA_DIR"
echo "Logs will be saved to: $LOG_DIR"
echo "Running with up to $MAX_JOBS parallel processes"
echo "Each experiment will be repeated $REPETITIONS times"

declare -a experiments=()
for query_file in "${QUERY_FILES[@]}"; do
    if [[ "$query_file" == *"job-light"* ]]; then
        dataset="job-light"
    elif [[ "$query_file" == *"job"* ]]; then
        dataset="job-complete"
    else
        dataset="unknown"
    fi
    
    for plan_enumerator in "${PLAN_ENUMERATORS[@]}"; do
        for collapse_func in "${COLLAPSE_FUNCTIONS[@]}"; do
            for topk_size in "${TOPK_SIZES[@]}"; do
                for learn_mode in "${LEARN_MODES[@]}"; do
                    if [ "$collapse_func" == "Uncertainty" ]; then
                        for uncertainty_impact in "${UNCERTAINTY_IMPACTS[@]}"; do
                            for rep in $(seq 1 $REPETITIONS); do
                                exp_name="${dataset}_${plan_enumerator}_${collapse_func}_topk${topk_size}_${learn_mode}_impact${uncertainty_impact}_run${rep}"
                                experiments+=("$query_file|$plan_enumerator|$collapse_func|$topk_size|$learn_mode|$uncertainty_impact|$exp_name")
                            done
                        done
                    else
                        for rep in $(seq 1 $REPETITIONS); do
                            exp_name="${dataset}_${plan_enumerator}_${collapse_func}_topk${topk_size}_${learn_mode}_run${rep}"
                            experiments+=("$query_file|$plan_enumerator|$collapse_func|$topk_size|$learn_mode||$exp_name")
                        done
                    fi
                done
            done
        done
    done
done

UNCERTAINTY_BASE=$((${#QUERY_FILES[@]} * ${#PLAN_ENUMERATORS[@]} * 1 * ${#TOPK_SIZES[@]} * ${#LEARN_MODES[@]} * ${#UNCERTAINTY_IMPACTS[@]} * $REPETITIONS))
OTHER_BASE=$((${#QUERY_FILES[@]} * ${#PLAN_ENUMERATORS[@]} * 2 * ${#TOPK_SIZES[@]} * ${#LEARN_MODES[@]} * $REPETITIONS))

echo "Total experiments: ${#experiments[@]} ($UNCERTAINTY_BASE Uncertainty experiments + $OTHER_BASE Mean/Upper experiments)"

run_experiment() {
    local config="$1"
    IFS='|' read -r query_file plan_enumerator collapse_func topk_size learn_mode uncertainty_impact exp_name <<< "$config"
    
    args=(
        "--data-layout" "Row"
        "--plan-enumerator" "$plan_enumerator"
        "--cardinality-estimator" "ExperimentalRange"
        "--backend" "Interpreter"
        "--cardinality-csv" "${DATA_DIR}/${exp_name}.csv"
        "--collapse-function" "$collapse_func"
        "--topk-size" "$topk_size"
    )
    
    if [ -n "$uncertainty_impact" ]; then
        args+=("--uncertainty-impact" "$uncertainty_impact")
    fi
    
    if [ "$learn_mode" == "no-learn" ]; then
        args+=("--no-learn-cardinalities")
    fi
    
    args+=("$query_file")
    
    echo "[$$] Running: $exp_name"
    
    if $PROGRAM "${args[@]}" > "${LOG_DIR}/${exp_name}.log" 2>&1; then
        echo "[$$] SUCCESS: $exp_name"
    else
        echo "[$$] ERROR: $exp_name (check ${exp_name}.log)"
    fi
}

export -f run_experiment
export PROGRAM DATA_DIR LOG_DIR

if command -v parallel >/dev/null 2>&1; then
    echo "Using GNU parallel for execution..."
    printf '%s\n' "${experiments[@]}" | parallel -j "$MAX_JOBS" run_experiment
elif command -v xargs >/dev/null 2>&1; then
    echo "Using xargs for parallel execution..."
    printf '%s\n' "${experiments[@]}" | xargs -n 1 -P "$MAX_JOBS" -I {} bash -c 'run_experiment "$@"' _ {}
else
    echo "Neither GNU parallel nor xargs found. Running sequentially..."
    for experiment in "${experiments[@]}"; do
        run_experiment "$experiment"
    done
fi

echo ""
echo "All experiments completed!"
echo "Results stored in: $BASE_OUTPUT_DIR"
echo ""
echo "Summary:"
echo "  CSV files: $(ls -1 ${DATA_DIR}/*.csv 2>/dev/null | wc -l)"
echo "  Log files: $(ls -1 ${LOG_DIR}/*.log 2>/dev/null | wc -l)"
echo ""
echo "Expected files per configuration: $REPETITIONS runs"
UNCERTAINTY_CONFIGS=$((${#QUERY_FILES[@]} * ${#PLAN_ENUMERATORS[@]} * 1 * ${#TOPK_SIZES[@]} * ${#LEARN_MODES[@]} * ${#UNCERTAINTY_IMPACTS[@]}))
OTHER_CONFIGS=$((${#QUERY_FILES[@]} * ${#PLAN_ENUMERATORS[@]} * 2 * ${#TOPK_SIZES[@]} * ${#LEARN_MODES[@]}))
echo "Total unique configurations: $((UNCERTAINTY_CONFIGS + OTHER_CONFIGS)) ($UNCERTAINTY_CONFIGS Uncertainty + $OTHER_CONFIGS Mean/Upper)"
echo ""
echo "Example files generated:"
echo "Data files:"
ls -1 ${DATA_DIR}/*.csv 2>/dev/null | head -5
echo "Log files:"
ls -1 ${LOG_DIR}/*.log 2>/dev/null | head -5
echo "..."