#!/bin/bash

OUTPUT_DIR="Z_MICRO_BENCHMARK"
if [ -d "$OUTPUT_DIR" ]; then
  rm -rf "$OUTPUT_DIR"
fi
mkdir -p "$OUTPUT_DIR"

PROGRAM="${PWD}/build/debug_shared/bin/shell"

QUERY_FILE="${PWD}/benchmark/job-light/ALL_QUERIES.sql"

echo "Running all configurations from launch.json..."

# echo "Running CartesianProduct..."
# $PROGRAM \
#     --no-statistics \
#     --no-learn-cardinalities  \
#     --cardinality-csv ${OUTPUT_DIR}/CartesianProduct.csv \
#     --plan-enumerator GOO \
#     --cardinality-estimator CartesianProduct \
#     --backend Interpreter \
#     "$QUERY_FILE" > "${OUTPUT_DIR}/CartesianProduct_raw.log" 2>&1
# echo "Output saved to ${OUTPUT_DIR}/CartesianProduct_raw.log"

# echo "Running CartesianProduct..."
# $PROGRAM \
#     --no-statistics \
#     --cardinality-csv ${OUTPUT_DIR}/CartesianProduct.csv \
#     --plan-enumerator GOO \
#     --cardinality-estimator CartesianProduct \
#     --backend Interpreter \
#     "$QUERY_FILE" > "${OUTPUT_DIR}/CartesianProduct_learned.log" 2>&1
# echo "Output saved to ${OUTPUT_DIR}/CartesianProduct_learned.log"


echo "Running SelectivityBased..."
$PROGRAM \
    --no-learn-cardinalities \
    --cardinality-csv ${OUTPUT_DIR}/Selecvitity.csv \
    --plan-enumerator GOO \
    --cardinality-estimator Selectivitybased \
    --backend Interpreter \
    --sample-size 100000 \
    "$QUERY_FILE" > "${OUTPUT_DIR}/SelectivityBased_raw.log" 2>&1
echo "Output saved to ${OUTPUT_DIR}/SelectivityBased_raw.log"

echo "Running SelectivityBased..."
$PROGRAM \
    --cardinality-csv ${OUTPUT_DIR}/Selecvitity.csv \
    --plan-enumerator GOO \
    --cardinality-estimator Selectivitybased \
    --backend Interpreter \
    --sample-size 100000 \
    "$QUERY_FILE" > "${OUTPUT_DIR}/SelectivityBased_learned.log" 2>&1
echo "Output saved to ${OUTPUT_DIR}/SelectivityBased_learned.log"

# Run a single dry run
echo "Running ExperimentalRangeEstimator..."
$PROGRAM \
    --plan-enumerator RangeGOO \
    --cardinality-estimator ExperimentalRange \
    --backend Interpreter \
    --sample-size 100000 \
    --cardinality-csv ${OUTPUT_DIR}/ExperimentalRangeEstimator_dry.csv \
    --collapse-function Mean \
    --no-learn-cardinalities \
    "$QUERY_FILE" > "${OUTPUT_DIR}/ExperimentalRangeEstimator_dry.log" 2>&1
echo "Output saved to ${OUTPUT_DIR}/ExperimentalRangeEstimator_dry.log"

echo "Running ExperimentalRangeEstimator..."
$PROGRAM \
    --plan-enumerator RangeGOO \
    --cardinality-estimator ExperimentalRange \
    --backend Interpreter \
    --sample-size 100000 \
    --cardinality-csv ${OUTPUT_DIR}/ExperimentalRangeEstimator_MEAN_0.csv \
    --collapse-function Mean \
    "$QUERY_FILE" > "${OUTPUT_DIR}/ExperimentalRangeEstimator_COLLAPSE_MEAN_0.log" 2>&1
echo "Output saved to ${OUTPUT_DIR}/ExperimentalRangeEstimator_COLLAPSE_MEAN_0.log"

echo "Running ExperimentalRangeEstimator..."
$PROGRAM \
    --plan-enumerator RangeGOO \
    --cardinality-estimator ExperimentalRange \
    --backend Interpreter \
    --sample-size 100000 \
    --cardinality-csv ${OUTPUT_DIR}/ExperimentalRangeEstimator_MEAN_1.csv \
    --collapse-function Mean \
    "$QUERY_FILE" > "${OUTPUT_DIR}/ExperimentalRangeEstimator_COLLAPSE_MEAN_1.log" 2>&1
echo "Output saved to ${OUTPUT_DIR}/ExperimentalRangeEstimator_COLLAPSE_MEAN_1.log"

echo "Running ExperimentalRangeEstimator..."
$PROGRAM \
    --plan-enumerator RangeGOO \
    --cardinality-estimator ExperimentalRange \
    --backend Interpreter \
    --sample-size 100000 \
    --cardinality-csv ${OUTPUT_DIR}/ExperimentalRangeEstimator_MEAN_2.csv \
    --collapse-function Mean \
    "$QUERY_FILE" > "${OUTPUT_DIR}/ExperimentalRangeEstimator_COLLAPSE_MEAN_2.log" 2>&1
echo "Output saved to ${OUTPUT_DIR}/ExperimentalRangeEstimator_COLLAPSE_MEAN_2.log"


echo "Running ExperimentalRangeEstimator..."
$PROGRAM \
    --plan-enumerator RangeGOO \
    --cardinality-estimator ExperimentalRange \
    --backend Interpreter \
    --sample-size 100000 \
    --cardinality-csv ${OUTPUT_DIR}/ExperimentalRangeEstimator_UPPER_0.csv \
    --collapse-function upperBound \
    "$QUERY_FILE" > "${OUTPUT_DIR}/ExperimentalRangeEstimator_COLLAPSE_UPPER_0.log" 2>&1
echo "Output saved to ${OUTPUT_DIR}/ExperimentalRangeEstimator_COLLAPSE_UPPER_0.log"

echo "Running ExperimentalRangeEstimator..."
$PROGRAM \
    --plan-enumerator RangeGOO \
    --cardinality-estimator ExperimentalRange \
    --backend Interpreter \
    --sample-size 100000 \
    --cardinality-csv ${OUTPUT_DIR}/ExperimentalRangeEstimator_UPPER_0.csv \
    --collapse-function upperBound \
    "$QUERY_FILE" > "${OUTPUT_DIR}/ExperimentalRangeEstimator_COLLAPSE_UPPER_1.log" 2>&1
echo "Output saved to ${OUTPUT_DIR}/ExperimentalRangeEstimator_COLLAPSE_UPPER_1.log"

echo "Running ExperimentalRangeEstimator..."
$PROGRAM \
    --plan-enumerator RangeGOO \
    --cardinality-estimator ExperimentalRange \
    --backend Interpreter \
    --sample-size 100000 \
    --cardinality-csv ${OUTPUT_DIR}/ExperimentalRangeEstimator_UPPER_0.csv \
    --collapse-function upperBound \
    "$QUERY_FILE" > "${OUTPUT_DIR}/ExperimentalRangeEstimator_COLLAPSE_UPPER_2.log" 2>&1
echo "Output saved to ${OUTPUT_DIR}/ExperimentalRangeEstimator_COLLAPSE_UPPER_2.log"


# echo "Running RangeGOO..."
# $PROGRAM \
#     --no-learn-cardinalities  \
#     --plan-enumerator RangeGOO \
#     --cardinality-estimator RangeCartesianProduct \
#     --backend Interpreter \
#     --cardinality-csv ${OUTPUT_DIR}/range_first_try.csv \
#     --collapse-function UpperBound \
#     "$QUERY_FILE" > "${OUTPUT_DIR}/RangeGOO_output.log" 2>&1
# echo "Output saved to ${OUTPUT_DIR}/RangeGOO_output.log"

# echo "Running Histogram..."
# $PROGRAM \
#     --no-learn-cardinalities  \
#     --cardinality-csv ${OUTPUT_DIR}/histogram.csv \
#     --plan-enumerator GOO \
#     --cardinality-estimator Histogram \
#     --backend Interpreter \
#     "$QUERY_FILE" > "${OUTPUT_DIR}/Histogram_output.log" 2>&1
# echo "Output saved to ${OUTPUT_DIR}/Histogram_output.log"

# echo "Running DPsizeOpt..."
# $PROGRAM \
#     --no-learn-cardinalities  \
#     --plan-enumerator DPsizeOpt \
#     --backend Interpreter \
#     --cardinality-csv ${OUTPUT_DIR}/dp_size_opt.csv \
#     --cardinality-estimator CartesianProduct \
#     "$QUERY_FILE" > "${OUTPUT_DIR}/DPsizeOpt_output.log" 2>&1
# echo "Output saved to ${OUTPUT_DIR}/DPsizeOpt_output.log"

# echo "Running DPsizeOpt..."
# $PROGRAM \
#     --no-learn-cardinalities  \
#     --plan-enumerator DPsizeOpt \
#     --backend Interpreter \
#     --cardinality-csv ${OUTPUT_DIR}/dp_size_opt.csv \
#     --cardinality-estimator Selectivitybased \
#     "$QUERY_FILE" > "${OUTPUT_DIR}/DPsizeOpt_selectivity.log" 2>&1
# echo "Output saved to ${OUTPUT_DIR}/DPsizeOpt_selectivity.log"

echo "All configurations completed. Check the ${OUTPUT_DIR} directory for output logs."