# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.22.0/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.22.0/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/scarrot/Dropbox/Developes/FPFHE_V1/FPFHE_mac

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/scarrot/Dropbox/Developes/FPFHE_V1/FPFHE_mac/build

# Include any dependencies generated for this target.
include CMakeFiles/Correct.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/Correct.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/Correct.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Correct.dir/flags.make

CMakeFiles/Correct.dir/testCorrect.cpp.o: CMakeFiles/Correct.dir/flags.make
CMakeFiles/Correct.dir/testCorrect.cpp.o: ../testCorrect.cpp
CMakeFiles/Correct.dir/testCorrect.cpp.o: CMakeFiles/Correct.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/scarrot/Dropbox/Developes/FPFHE_V1/FPFHE_mac/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Correct.dir/testCorrect.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/Correct.dir/testCorrect.cpp.o -MF CMakeFiles/Correct.dir/testCorrect.cpp.o.d -o CMakeFiles/Correct.dir/testCorrect.cpp.o -c /Users/scarrot/Dropbox/Developes/FPFHE_V1/FPFHE_mac/testCorrect.cpp

CMakeFiles/Correct.dir/testCorrect.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Correct.dir/testCorrect.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/scarrot/Dropbox/Developes/FPFHE_V1/FPFHE_mac/testCorrect.cpp > CMakeFiles/Correct.dir/testCorrect.cpp.i

CMakeFiles/Correct.dir/testCorrect.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Correct.dir/testCorrect.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/scarrot/Dropbox/Developes/FPFHE_V1/FPFHE_mac/testCorrect.cpp -o CMakeFiles/Correct.dir/testCorrect.cpp.s

# Object files for target Correct
Correct_OBJECTS = \
"CMakeFiles/Correct.dir/testCorrect.cpp.o"

# External object files for target Correct
Correct_EXTERNAL_OBJECTS =

Correct: CMakeFiles/Correct.dir/testCorrect.cpp.o
Correct: CMakeFiles/Correct.dir/build.make
Correct: /usr/local/lib/libPALISADEpke.1.11.5.dylib
Correct: /usr/local/lib/libPALISADEbinfhe.1.11.5.dylib
Correct: /usr/local/lib/libPALISADEfpfhe.1.11.5.dylib
Correct: /usr/local/lib/libPALISADEcore.1.11.5.dylib
Correct: /usr/local/lib/libhexl.1.2.3.dylib
Correct: CMakeFiles/Correct.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/scarrot/Dropbox/Developes/FPFHE_V1/FPFHE_mac/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Correct"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Correct.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Correct.dir/build: Correct
.PHONY : CMakeFiles/Correct.dir/build

CMakeFiles/Correct.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Correct.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Correct.dir/clean

CMakeFiles/Correct.dir/depend:
	cd /Users/scarrot/Dropbox/Developes/FPFHE_V1/FPFHE_mac/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/scarrot/Dropbox/Developes/FPFHE_V1/FPFHE_mac /Users/scarrot/Dropbox/Developes/FPFHE_V1/FPFHE_mac /Users/scarrot/Dropbox/Developes/FPFHE_V1/FPFHE_mac/build /Users/scarrot/Dropbox/Developes/FPFHE_V1/FPFHE_mac/build /Users/scarrot/Dropbox/Developes/FPFHE_V1/FPFHE_mac/build/CMakeFiles/Correct.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Correct.dir/depend

