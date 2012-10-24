#!/bin/bash
# rm all the temporary files created by cmake
# usage : ./cmake-clean.sh build
REP=$1
cd $REP
make clean
rm -rf CMakeFiles/ CMakeCache.txt Makefile DartConfiguration.tcl CTestTestfile.cmake cmake_install.cmake Testing
cd ..