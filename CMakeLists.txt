cmake_minimum_required(VERSION 3.18)
project(team19)

set(CMAKE_FIND_FRAMEWORK LAST)
set(CMAKE_CXX_STANDARD 14)

include_directories(src)
add_executable(team19
        main.cpp
        src/PngReader.cpp
        src/PngReader.h
        src/ImageQuilting.cpp
        src/ImageQuilting.h
        src/ImgData.h

        src/benchmarking/tsc_x86.h
        src/benchmarking/timing.h
        src/benchmarking/timing.cpp)

find_package(PNG 1.6 REQUIRED)
include_directories(${PNG_INCLUDE_DIR})
target_link_libraries(team19 ${PNG_LIBRARY})

