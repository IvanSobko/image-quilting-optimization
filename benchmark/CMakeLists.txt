CMAKE_MINIMUM_REQUIRED(VERSION 3.10)

PROJECT(IMAGE-QUILTING)

ADD_SUBDIRECTORY(lib/agz-utils)
TARGET_COMPILE_DEFINITIONS(AGZUtils PUBLIC AGZ_UTILS_SSE)

FILE(GLOB_RECURSE SRC
		"${PROJECT_SOURCE_DIR}/src/*.cpp"
		"${PROJECT_SOURCE_DIR}/src/*.h")

ADD_EXECUTABLE(ImageQuilting ${SRC})

FOREACH(_SRC IN ITEMS ${SRC})
    GET_FILENAME_COMPONENT(SRC "${_SRC}" PATH)
    STRING(REPLACE "${PROJECT_SOURCE_DIR}/src" "src" _GRP_PATH "${SRC}")
    STRING(REPLACE "/" "\\" _GRP_PATH "${_GRP_PATH}")
    SOURCE_GROUP("${_GRP_PATH}" FILES "${_SRC}")
ENDFOREACH()

# We can use different compiler flags
set(CompileFlags "-O3" "-ffast-math" "-march=native")
add_definitions(-D_CompileFlags="-O3-ffast-math-march=native") # pass variable _CompileFlags to c++ code

#set(CompileFlags "-O3" "-fno-tree-vectorize")
#add_definitions(-D_CompileFlags="-O3-fno-tree-vectorize")

#set(CompileFlags "-O1")
#add_definitions(-D_CompileFlags="-O1")

SET_PROPERTY(TARGET ImageQuilting PROPERTY CXX_STANDARD 17)
SET_PROPERTY(TARGET ImageQuilting PROPERTY CXX_STANDARD_REQUIRED ON)

TARGET_INCLUDE_DIRECTORIES(ImageQuilting PUBLIC "${PROJECT_SOURCE_DIR}/lib/cxxopts")
TARGET_LINK_LIBRARIES(ImageQuilting PUBLIC AGZUtils)
