# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.25

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/baptiste/team19/src/benchmarking/other_implementation

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/baptiste/team19/src/benchmarking/other_implementation/build_high

# Include any dependencies generated for this target.
include CMakeFiles/ImageQuilting.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/ImageQuilting.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/ImageQuilting.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/ImageQuilting.dir/flags.make

CMakeFiles/ImageQuilting.dir/src/main.cpp.o: CMakeFiles/ImageQuilting.dir/flags.make
CMakeFiles/ImageQuilting.dir/src/main.cpp.o: /home/baptiste/team19/src/benchmarking/other_implementation/src/main.cpp
CMakeFiles/ImageQuilting.dir/src/main.cpp.o: CMakeFiles/ImageQuilting.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/baptiste/team19/src/benchmarking/other_implementation/build_high/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/ImageQuilting.dir/src/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/ImageQuilting.dir/src/main.cpp.o -MF CMakeFiles/ImageQuilting.dir/src/main.cpp.o.d -o CMakeFiles/ImageQuilting.dir/src/main.cpp.o -c /home/baptiste/team19/src/benchmarking/other_implementation/src/main.cpp

CMakeFiles/ImageQuilting.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ImageQuilting.dir/src/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/baptiste/team19/src/benchmarking/other_implementation/src/main.cpp > CMakeFiles/ImageQuilting.dir/src/main.cpp.i

CMakeFiles/ImageQuilting.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ImageQuilting.dir/src/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/baptiste/team19/src/benchmarking/other_implementation/src/main.cpp -o CMakeFiles/ImageQuilting.dir/src/main.cpp.s

CMakeFiles/ImageQuilting.dir/src/imageQuilting.cpp.o: CMakeFiles/ImageQuilting.dir/flags.make
CMakeFiles/ImageQuilting.dir/src/imageQuilting.cpp.o: /home/baptiste/team19/src/benchmarking/other_implementation/src/imageQuilting.cpp
CMakeFiles/ImageQuilting.dir/src/imageQuilting.cpp.o: CMakeFiles/ImageQuilting.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/baptiste/team19/src/benchmarking/other_implementation/build_high/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/ImageQuilting.dir/src/imageQuilting.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/ImageQuilting.dir/src/imageQuilting.cpp.o -MF CMakeFiles/ImageQuilting.dir/src/imageQuilting.cpp.o.d -o CMakeFiles/ImageQuilting.dir/src/imageQuilting.cpp.o -c /home/baptiste/team19/src/benchmarking/other_implementation/src/imageQuilting.cpp

CMakeFiles/ImageQuilting.dir/src/imageQuilting.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ImageQuilting.dir/src/imageQuilting.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/baptiste/team19/src/benchmarking/other_implementation/src/imageQuilting.cpp > CMakeFiles/ImageQuilting.dir/src/imageQuilting.cpp.i

CMakeFiles/ImageQuilting.dir/src/imageQuilting.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ImageQuilting.dir/src/imageQuilting.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/baptiste/team19/src/benchmarking/other_implementation/src/imageQuilting.cpp -o CMakeFiles/ImageQuilting.dir/src/imageQuilting.cpp.s

# Object files for target ImageQuilting
ImageQuilting_OBJECTS = \
"CMakeFiles/ImageQuilting.dir/src/main.cpp.o" \
"CMakeFiles/ImageQuilting.dir/src/imageQuilting.cpp.o"

# External object files for target ImageQuilting
ImageQuilting_EXTERNAL_OBJECTS =

ImageQuilting: CMakeFiles/ImageQuilting.dir/src/main.cpp.o
ImageQuilting: CMakeFiles/ImageQuilting.dir/src/imageQuilting.cpp.o
ImageQuilting: CMakeFiles/ImageQuilting.dir/build.make
ImageQuilting: lib/agz-utils/libAGZUtils.a
ImageQuilting: CMakeFiles/ImageQuilting.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/baptiste/team19/src/benchmarking/other_implementation/build_high/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable ImageQuilting"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ImageQuilting.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ImageQuilting.dir/build: ImageQuilting
.PHONY : CMakeFiles/ImageQuilting.dir/build

CMakeFiles/ImageQuilting.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ImageQuilting.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ImageQuilting.dir/clean

CMakeFiles/ImageQuilting.dir/depend:
	cd /home/baptiste/team19/src/benchmarking/other_implementation/build_high && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/baptiste/team19/src/benchmarking/other_implementation /home/baptiste/team19/src/benchmarking/other_implementation /home/baptiste/team19/src/benchmarking/other_implementation/build_high /home/baptiste/team19/src/benchmarking/other_implementation/build_high /home/baptiste/team19/src/benchmarking/other_implementation/build_high/CMakeFiles/ImageQuilting.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/ImageQuilting.dir/depend
