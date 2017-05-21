#!/usr/bin/env bash

# delete travis_env_paths file if it already exists
if [ -f "travis_env_paths" ]; then
    rm travis_env_paths
fi

set -e

install_brew_package() {
  if brew list -1 | grep -q "^$1\$"; then
    # Package is installed, upgrade if needed
    brew outdated "$1" || brew upgrade "$@"
  else
    # Package not installed yet, install.
    # If there are conflicts, try overwriting the files (these are in /usr/local anyway so it should be ok).
    brew install "$@" || brew link --overwrite gcc49    
  fi
}

brew update

# For md5sum
install_brew_package md5sha1sum
# For `timeout'
install_brew_package coreutils
# For analysis
install_brew_package cppcheck
# build script needs bash version >= 4.0
install_brew_package bash

sudo chsh -s /usr/local/bin/bash

# build script needs gnu getopt
install_brew_package gnu-getopt
brew link --force gnu-getopt

if [[ "${INSTALL_VALGRIND}" == "1" ]]
then
    install_brew_package valgrind
fi

which cmake &>/dev/null || install_brew_package cmake

case "${COMPILER}" in
gcc-4.8)
    install_brew_package gcc@4.8 
    echo "$(brew --prefix llvm@3.9)/bin" >> travis_env_paths
    ;;
gcc-4.9)
    install_brew_package gcc@4.9
    echo "$(brew --prefix gcc@4.9)/bin" >> travis_env_paths
    ;;
gcc-5)
    install_brew_package gcc@5 
    echo "$(brew --prefix gcc@5)/bin" >> travis_env_paths
    ;;
gcc-6)         
    install_brew_package gcc@6 
    echo "$(brew --prefix gcc@6)/bin" >> travis_env_paths
    ;;
clang-default) ;;
clang-3.7)     
    install_brew_package llvm@3.7 --with-clang --with-libcxx
    echo "$(brew --prefix llvm@3.7)/bin" >> travis_env_paths
    ;;
clang-3.8)     
    install_brew_package llvm@3.8 --with-clang --with-libcxx
    echo "$(brew --prefix llvm@3.8)/bin" >> travis_env_paths
    ;;
clang-3.9)     
    install_brew_package llvm@3.9 --with-clang --with-libcxx
    echo "$(brew --prefix llvm@3.9)/bin" >> travis_env_paths
    ;;
clang-4.0)     
    install_brew_package llvm     --with-clang --with-libcxx
    echo "$(brew --prefix llvm)/bin" >> travis_env_paths
    ;;
*) echo "Compiler not supported: ${COMPILER}. See travis_ci_install_osx.sh"; exit 1 ;;
esac

brew install python3
pip3 install pytest
pip3 install pytest-xdist
pip3 install sh
