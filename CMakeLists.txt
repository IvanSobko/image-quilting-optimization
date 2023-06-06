cmake_minimum_required(VERSION 3.18)
project(team19)

set(CMAKE_FIND_FRAMEWORK LAST)
set(CMAKE_CXX_STANDARD 17)

include_directories(src)
add_executable(team19
        main.cpp
        src/PngReader.cpp
        src/PngReader.h
        src/ImageQuilting.cpp
        src/ImageQuilting.h
        src/AdvanceAlgOptimiz.cpp
        src/AdvanceAlgOptimiz.h
        src/CompOverlapOptimiz.cpp
        src/CompOverlapOptimiz.h
        src/ImgData.h
        src/Testing.cpp
        src/Testing.h
        src/benchmarking/tsc_x86.h
        src/benchmarking/timing.h
        src/benchmarking/timing.cpp
        src/Blocking.h
        src/Blocking.cpp)

# We can use different compiler flags
#set(CompileFlags "-g" "-O3" "-ffast-math" "-march=native")
#add_definitions(-D_CompileFlags="-O3-ffast-math-march=native") # pass variable _CompileFlags to c++ code

#set(CompileFlags "-g" "-O3" "-fno-tree-vectorize")
#add_definitions(-D_CompileFlags="-O3-fno-tree-vectorize")

#set(CompileFlags "-g" "-O1")
#add_definitions(-D_CompileFlags="-O1")

#target_compile_options(team19 PRIVATE ${CompileFlags})

find_package(PNG 1.6 REQUIRED)
include_directories(${PNG_INCLUDE_DIR})
target_link_libraries(team19 ${PNG_LIBRARY})

