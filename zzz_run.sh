#!/bin/bash
# filepath: /Users/oliver/TU_BERLIN/MASTER/mutable/extract_errors.sh

# Remove existing log files
rm -f zzz.log zzz_error.log

# Run the test script and redirect output to zzz.log
./test_run_full.sh >> zzz.log

# Extract error lines with 2 lines of context before each into zzz_error.log
grep -B2 "error:" zzz.log > zzz_error.log

# Also include other common error indicators
grep -B2 "undefined" zzz.log >> zzz_error.log
grep -B2 "FAILED" zzz.log >> zzz_error.log

# Remove duplicate separator lines to make output more readable
sed -i '' '/^--$/d' zzz_error.log

echo "Errors have been extracted to zzz_error.log"