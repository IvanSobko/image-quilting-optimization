#!/bin/bash

# Note: this script will not run on Windows.
# It should run on MacOS and Linux if bash is installed.

GAL_DIR="../../../gallery"
RES_FILE="../../../results/timing_results_-other_implementation_O3-ffast-math-march=native.txt"

cmake .
cmake --build .
./ImageQuilting --input=$GAL_DIR/input0_24x24.png --output=$GAL_DIR/output.png --blockW=12 --blockH=12 --width=48 --height=48 >> $RES_FILE
./ImageQuilting --input=$GAL_DIR/input0_48x48.png --output=$GAL_DIR/output.png --blockW=24 --blockH=24 --width=96 --height=96 >> $RES_FILE
./ImageQuilting --input=$GAL_DIR/input0_96x96.png --output=$GAL_DIR/output.png --blockW=48 --blockH=48 --width=192 --height=192 >> $RES_FILE
./ImageQuilting --input=$GAL_DIR/input0_192x192.png --output=$GAL_DIR/output.png --blockW=96 --blockH=96 --width=384 --height=384 >> $RES_FILE
./ImageQuilting --input=$GAL_DIR/input0_384x384.png --output=$GAL_DIR/output.png --blockW=192 --blockH=192 --width=768 --height=768 >> $RES_FILE
./ImageQuilting --input=$GAL_DIR/input0_768x768.png --output=$GAL_DIR/output.png --blockW=384 --blockH=384 --width=1536 --height=1536 >> $RES_FILE

rm -rf CMakeCache.txt CMakeFiles/ cmake_install.cmake ImageQuilting Makefile 