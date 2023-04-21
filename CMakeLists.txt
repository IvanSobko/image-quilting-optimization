cmake_minimum_required(VERSION 3.23)
project(team19)

set(CMAKE_CXX_STANDARD 14)

add_executable(team19
        main.cpp
        src/PngReader.cpp
        src/PngReader.h
        src/ImageQuilting.cpp
        src/ImageQuilting.h
        src/ImgData.h)

find_package(PNG 1.6.39 EXACT REQUIRED)
include_directories(${PNG_INCLUDE_DIR})
target_link_libraries(team19 ${PNG_LIBRARY})

