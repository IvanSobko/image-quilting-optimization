cmake_minimum_required(VERSION 3.18)
project(ImageQuilting)

add_subdirectory(lib/agz-utils)

set(CMAKE_FIND_FRAMEWORK LAST)
set(CMAKE_CXX_STANDARD 17)

include_directories(src)
add_executable(ImageQuilting
        src/main.cpp
        src/imageQuilting.cpp
        src/imageQuilting.h
        src/tsc_x86.h)

# We can use different compiler flags
set(CompileFlags "-g" "-O3" "-ffast-math" "-march=native")
add_definitions(-D_CompileFlags="-O3-ffast-math-march=native") # pass variable _CompileFlags to c++ code

#set(CompileFlags "-g" "-O3" "-fno-tree-vectorize")
#add_definitions(-D_CompileFlags="-O3-fno-tree-vectorize")

#set(CompileFlags "-g" "-O1")
#add_definitions(-D_CompileFlags="-O1")

target_compile_definitions(AGZUtils PUBLIC AGZ_UTILS_SSE)
target_compile_options(ImageQuilting PRIVATE ${CompileFlags})

include_directories(ImageQuilting PUBLIC "${PROJECT_SOURCE_DIR}/lib/cxxopts")
target_link_libraries(ImageQuilting PUBLIC AGZUtils)
