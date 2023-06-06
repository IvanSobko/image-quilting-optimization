#!/bin/bash

# Set the different compile flag sets
flag_sets=(
  "-O1"
  "-O3 -fno-tree-vectorize"
  "-O3 -ffast-math -march=native"
)

# Set the different compile flag set names
flag_set_names=(
  "-O1"
  "-O3-fno-tree-vectorize"
  "-O3-ffast-math-march=native"
)

# Loop through the index of flag sets
for index in "${!flag_sets[@]}"
do
  # cd into the build directory
  cd build

  # Get the current compile flag set and name
  flags="${flag_sets[index]}"
  name="${flag_set_names[index]}"

  # Echo state
  echo "Building project with \"$flags\":"

  # Configure CMake with the specific compile flags and define _CompileFlags accordingly
  if [ "$index" != 2 ]; then
    echo "cmake" "$project_dir" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="\"$flags\" -D_CompileFlags=\"\\\""$name"\\\"\"" ".."
    cmake"$project_dir" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="\"$flags\" -D_CompileFlags=\"\\\""$name"\\\"\"" ".."
  else
    echo "cmake" "$project_dir" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="\"$flags\" -D_CompileFlags=\"\\\""$name"\\\"\" -DVECTORIZED" ".."
    cmake"$project_dir" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="\"$flags\" -D_CompileFlags=\"\\\""$name"\\\"\" -DVECTORIZED" ".."
  fi


  # Build the project
  cmake --build . -j

  # Run the project with the functional timing argument
  echo
  cd ..
  ./build/team19 --timingFunctional
done
