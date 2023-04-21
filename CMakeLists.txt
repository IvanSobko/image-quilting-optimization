cmake_minimum_required(VERSION 3.23)
project(team19)

set(CMAKE_CXX_STANDARD 14)

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/libs/libpng/include/libpng16/)
#
#link_directories(${CMAKE_CURRENT_SOURCE_DIR}/libs/libpng/lib)

add_executable(team19
        main.cpp
        src/PngReader.cpp
        src/PngReader.h
        src/ImageQuilting.cpp
        src/ImageQuilting.h
        src/ImgData.h)

list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/libs/libpng")
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/libs/zlib")

find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})
target_link_libraries(team19 ${PNG_LIBRARY})



#target_link_libraries (team19 libpng16)
#
#target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libs/libpng/include/libpng16/)
#target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libs/libpng/lib)
#target_link_libraries(${PROJECT_NAME} PUBLIC
#        ${CMAKE_CURRENT_SOURCE_DIR}/libs/libpng/lib/libpng16.lib)
