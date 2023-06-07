#!/bin/bash

# Note: this script will not run on Windows.
# It should run on MacOS and Linux if bash is installed.

GAL_DIR="../../../timing/input"
RES_FILE_HIGH="../../../timing/results/baseline_-O3-ffast-math-march=native.txt"
RES_FILE_MID="../../../timing/results/baseline_-O3-fno-tree-vectorize.txt"
RES_FILE_LOW="../../../timing/results/baseline_-O1.txt"

<<"COMMENT" # Those values raise errors (too low ?)
./build_low/ImageQuilting --input=$GAL_DIR/input0_32.png --output=output.png --blockW=8 --blockH=8 --width=64 --height=64 > $RES_FILE_LOW
./build_mid/ImageQuilting --input=$GAL_DIR/input0_32.png --output=output.png --blockW=8 --blockH=8 --width=64 --height=64 > $RES_FILE_MID
./build_high/ImageQuilting --input=$GAL_DIR/input0_32.png --output=output.png --blockW=8 --blockH=8 --width=64 --height=64 > $RES_FILE_HIGH

./build_low/ImageQuilting --input=$GAL_DIR/input0_42.png --output=output.png --blockW=10 --blockH=10 --width=84 --height=84 >> $RES_FILE_LOW
./build_mid/ImageQuilting --input=$GAL_DIR/input0_42.png --output=output.png --blockW=10 --blockH=10 --width=84 --height=84 >> $RES_FILE_MID
./build_high/ImageQuilting --input=$GAL_DIR/input0_42.png --output=output.png --blockW=10 --blockH=10 --width=84 --height=84 >> $RES_FILE_HIGH
COMMENT

./build_low/ImageQuilting --input=$GAL_DIR/input0_55.png --output=output.png --blockW=13 --blockH=13 --width=110 --height=110 > $RES_FILE_LOW
./build_mid/ImageQuilting --input=$GAL_DIR/input0_55.png --output=output.png --blockW=13 --blockH=13 --width=110 --height=110 > $RES_FILE_MID
./build_high/ImageQuilting --input=$GAL_DIR/input0_55.png --output=output.png --blockW=13 --blockH=13 --width=110 --height=110 > $RES_FILE_HIGH

./build_low/ImageQuilting --input=$GAL_DIR/input0_72.png --output=output.png --blockW=18 --blockH=18 --width=144 --height=144 >> $RES_FILE_LOW
./build_mid/ImageQuilting --input=$GAL_DIR/input0_72.png --output=output.png --blockW=18 --blockH=18 --width=144 --height=144 >> $RES_FILE_MID
./build_high/ImageQuilting --input=$GAL_DIR/input0_72.png --output=output.png --blockW=18 --blockH=18 --width=144 --height=144 >> $RES_FILE_HIGH

./build_low/ImageQuilting --input=$GAL_DIR/input0_95.png --output=output.png --blockW=23 --blockH=23 --width=190 --height=190 >> $RES_FILE_LOW
./build_mid/ImageQuilting --input=$GAL_DIR/input0_95.png --output=output.png --blockW=23 --blockH=23 --width=190 --height=190 >> $RES_FILE_MID
./build_high/ImageQuilting --input=$GAL_DIR/input0_95.png --output=output.png --blockW=23 --blockH=23 --width=190 --height=190 >> $RES_FILE_HIGH

./build_low/ImageQuilting --input=$GAL_DIR/input0_124.png --output=output.png --blockW=31 --blockH=31 --width=248 --height=248 >> $RES_FILE_LOW
./build_mid/ImageQuilting --input=$GAL_DIR/input0_124.png --output=output.png --blockW=31 --blockH=31 --width=248 --height=248 >> $RES_FILE_MID
./build_high/ImageQuilting --input=$GAL_DIR/input0_124.png --output=output.png --blockW=31 --blockH=31 --width=248 --height=248 >> $RES_FILE_HIGH

./build_low/ImageQuilting --input=$GAL_DIR/input0_162.png --output=output.png --blockW=40 --blockH=40 --width=324 --height=324 >> $RES_FILE_LOW
./build_mid/ImageQuilting --input=$GAL_DIR/input0_162.png --output=output.png --blockW=40 --blockH=40 --width=324 --height=324 >> $RES_FILE_MID
./build_high/ImageQuilting --input=$GAL_DIR/input0_162.png --output=output.png --blockW=40 --blockH=40 --width=324 --height=324 >> $RES_FILE_HIGH

./build_low/ImageQuilting --input=$GAL_DIR/input0_212.png --output=output.png --blockW=53 --blockH=53 --width=424 --height=424 >> $RES_FILE_LOW
./build_mid/ImageQuilting --input=$GAL_DIR/input0_212.png --output=output.png --blockW=53 --blockH=53 --width=424 --height=424 >> $RES_FILE_MID
./build_high/ImageQuilting --input=$GAL_DIR/input0_212.png --output=output.png --blockW=53 --blockH=53 --width=424 --height=424 >> $RES_FILE_HIGH

./build_low/ImageQuilting --input=$GAL_DIR/input0_277.png --output=output.png --blockW=69 --blockH=69 --width=554 --height=554 >> $RES_FILE_LOW
./build_mid/ImageQuilting --input=$GAL_DIR/input0_277.png --output=output.png --blockW=69 --blockH=69 --width=554 --height=554 >> $RES_FILE_MID
./build_high/ImageQuilting --input=$GAL_DIR/input0_277.png --output=output.png --blockW=69 --blockH=69 --width=554 --height=554 >> $RES_FILE_HIGH

./build_low/ImageQuilting --input=$GAL_DIR/input0_363.png --output=output.png --blockW=90 --blockH=90 --width=726 --height=726 >> $RES_FILE_LOW
./build_mid/ImageQuilting --input=$GAL_DIR/input0_363.png --output=output.png --blockW=90 --blockH=90 --width=726 --height=726 >> $RES_FILE_MID
./build_high/ImageQuilting --input=$GAL_DIR/input0_363.png --output=output.png --blockW=90 --blockH=90 --width=726 --height=726 >> $RES_FILE_HIGH

rm -rf CMakeCache.txt CMakeFiles/ cmake_install.cmake ImageQuilting Makefile output.png