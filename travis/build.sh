#!/usr/bin/env bash

SCRIPT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
PROJECT_NAME=$1
OS=$2
shift
shift
ARGS=$@

set -e

case $OS in
docker-*)
    CPP_BASE_BUILDSYSTEM_TAG=$(echo $OS | sed 's/docker-//g') 
    python3 ${SCRIPT_DIR}/../build-tools/cmakew_docker --container-name ${PROJECT_NAME} --docker-image-tag ${CPP_BASE_BUILDSYSTEM_TAG} ${ARGS}
    exit $?
    ;;
macos | linux | windows | native )
    bash ${SCRIPT_DIR}/../build-tools/cmakew ${ARGS}
    exit $?
    ;;
*)
    echo "Unsupported OS: $OS"
    exit 1
esac
