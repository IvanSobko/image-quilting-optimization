cmake_minimum_required(VERSION 3.23)
project(team19)

set(CMAKE_CXX_STANDARD 14)

add_executable(team19
        main.cpp
        src/pngReader.cpp
        src/pngReader.h)

find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})
target_link_libraries(team19 ${PNG_LIBRARY})