#!/bin/bash
# filepath: /Users/oliver/TU_BERLIN/MASTER/mutable/zzz_build_server.sh
set -e

export CCACHE_DIR="$HOME/.ccache"
mkdir -p "$CCACHE_DIR"

ccache --set-config sloppiness=file_macro,locale,time_macros,include_file_ctime,include_file_mtime
ccache --set-config compression=false
ccache --set-config compression_level=6
ccache --set-config max_size=5G
echo "Using ccache config:"
ccache -p | grep -E 'sloppiness|compression'

BOOST_PATH="$HOME/boost_1_89_0"
BUILD_DIR="build/debug_shared"
BUILD_LOG="build_full.log"
ERROR_LOG="build_errors.log"

EXTRA_CMAKE_ARGS="$@"

export CC="clang"
export CXX="clang++"
echo "Using system compilers: $(which clang) and $(which clang++)"
echo "Clang version: $(clang --version | head -1)"

if [ -f "$BUILD_DIR/CMakeCache.txt" ]; then
    rm -f "$BUILD_DIR/CMakeCache.txt"
fi

echo "Starting build... Full log: $BUILD_LOG, Errors: $ERROR_LOG"

# Run cmake and capture output
cmake -S . -B "$BUILD_DIR" \
    -G Ninja \
    -DCMAKE_C_COMPILER="$CC" \
    -DCMAKE_CXX_COMPILER="$CXX" \
    -DCMAKE_C_COMPILER_LAUNCHER=ccache \
    -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
    -DBUILD_SHARED_LIBS=ON \
    -DWITH_V8=OFF \
    -DENABLE_SANITIZERS=OFF \
    -DENABLE_SANITY_FIELDS=OFF \
    -DBOOST_ROOT="$BOOST_PATH" \
    -DCMAKE_CXX_FLAGS="-stdlib=libstdc++" \
    -DCMAKE_EXE_LINKER_FLAGS="-stdlib=libstdc++" \
    $EXTRA_CMAKE_ARGS 2>&1 | tee "$BUILD_LOG"

ninja -C "$BUILD_DIR" -v 2>&1 | tee -a "$BUILD_LOG"
echo "Extracting errors with context to $ERROR_LOG..."
grep -B 5 -i "error\|failed\|fatal" "$BUILD_LOG" > "$ERROR_LOG" 2>/dev/null || echo "No errors found in build log" > "$ERROR_LOG"

echo ""
echo "Build complete!"
echo "Full build log: $BUILD_LOG"
echo "Error summary: $ERROR_LOG"

if [ -s "$ERROR_LOG" ]; then
    echo ""
    echo "=== BUILD ERRORS FOUND ==="
    cat "$ERROR_LOG"
    echo "=========================="
fi

ccache --show-stats --verbose