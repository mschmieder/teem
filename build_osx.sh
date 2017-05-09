#!/usr/bin/env bash
SCRIPT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
mkdir -p build
cd build

OPTS=()
OPTS+=("--source-directory=${SCRIPT_DIR}")
OPTS+=("--build-directory=${SCRIPT_DIR}/build/")
#OPTS+=("--clean")
OPTS+=("--target-system=macos")
OPTS+=("--compiler=clang-3.8")
OPTS+=("--target-architecture=x86_64")
OPTS+=("--build-type=release")
OPTS+=("--link-type=shared")
#OPTS+=("--analyze")
#OPTS+=("--add-analyzer cppcheck")
OPTS+=("--coverage")
OPTS+=("--gcovr-output-formatter=html")
OPTS+=("--build")
OPTS+=("--test")
#OPTS+=("--install")
#OPTS+=("--install-directory=${SCRIPT_DIR}/build/install")
OPTS+=("--num-threads=4")
OPTS+=("--verbose")

../travis/build.sh "teem" "native" ${OPTS[@]}
