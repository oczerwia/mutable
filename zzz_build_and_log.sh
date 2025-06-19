#!/bin/bash

rm zzz.log && ./test_run_full.sh >> zzz.log

grep -i -A 2 -B 1 "error" /Users/oliver/TU_BERLIN/MASTER/mutable/zzz.log >> zzz_error.log