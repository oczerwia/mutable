#!/bin/bash

pipenv sync # Install python


# fake the versioning
git tag v0.0.0
git push origin v0.0.0

# Fake the tbl file
echo "X(GIT_REV, a828cd0a6c5a92966417d4eddf2ca52409dab2e7)" > /Users/oliver/TU_BERLIN/MASTER/mutable/include/mutable/gitversion.tbl
echo "X(GIT_BRANCH, main)" >> /Users/oliver/TU_BERLIN/MASTER/mutable/include/mutable/gitversion.tbl
echo "X(SEM_VERSION, v0.0.0)" >> /Users/oliver/TU_BERLIN/MASTER/mutable/include/mutable/gitversion.tbl

# Use ccache for debugging
export CC="ccache $(brew --prefix llvm@17)/bin/clang"
export CXX="ccache $(brew --prefix llvm@17)/bin/clang++"


# Use csvclean to clean out faulty lines in job-light benchmark

touch benchmark/job-ligh/data/movie_companies_cleaned.csv
touch benchmark/job-ligh/data/movie_info_idx_cleaned.csv

csvclean -a benchmark/job-light/data/movie_companies.csv >> benchmark/job-ligh/data/movie_companies_cleaned.csv
csvclean -a benchmark/job-light/data/movie_info_idx.csv >> benchmark/job-ligh/data/movie_info_idx_cleaned.csv