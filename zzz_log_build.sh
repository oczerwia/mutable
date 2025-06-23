#!/bin/bash

chmod +x ./make_and_build.sh

rm -f zzz.log && ./make_and_build.sh >> zzz.log

rm -f zzz_error.log && grep -B 2 -i error zzz.log > zzz_error.log
