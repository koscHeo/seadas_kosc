# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
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
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jaemoo/seadas-7.4/ocssw/build

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jaemoo/seadas-7.4/ocssw/build/cbuild

# Include any dependencies generated for this target.
include src/l3bin/CMakeFiles/l3binmerge.dir/depend.make

# Include the progress variables for this target.
include src/l3bin/CMakeFiles/l3binmerge.dir/progress.make

# Include the compile flags for this target's objects.
include src/l3bin/CMakeFiles/l3binmerge.dir/flags.make

src/l3bin/CMakeFiles/l3binmerge.dir/l3binmerge.cpp.o: src/l3bin/CMakeFiles/l3binmerge.dir/flags.make
src/l3bin/CMakeFiles/l3binmerge.dir/l3binmerge.cpp.o: ../src/l3bin/l3binmerge.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/jaemoo/seadas-7.4/ocssw/build/cbuild/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/l3bin/CMakeFiles/l3binmerge.dir/l3binmerge.cpp.o"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l3bin && /usr/bin/g++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/l3binmerge.dir/l3binmerge.cpp.o -c /home/jaemoo/seadas-7.4/ocssw/build/src/l3bin/l3binmerge.cpp

src/l3bin/CMakeFiles/l3binmerge.dir/l3binmerge.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/l3binmerge.dir/l3binmerge.cpp.i"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l3bin && /usr/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/jaemoo/seadas-7.4/ocssw/build/src/l3bin/l3binmerge.cpp > CMakeFiles/l3binmerge.dir/l3binmerge.cpp.i

src/l3bin/CMakeFiles/l3binmerge.dir/l3binmerge.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/l3binmerge.dir/l3binmerge.cpp.s"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l3bin && /usr/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/jaemoo/seadas-7.4/ocssw/build/src/l3bin/l3binmerge.cpp -o CMakeFiles/l3binmerge.dir/l3binmerge.cpp.s

src/l3bin/CMakeFiles/l3binmerge.dir/l3binmerge.cpp.o.requires:
.PHONY : src/l3bin/CMakeFiles/l3binmerge.dir/l3binmerge.cpp.o.requires

src/l3bin/CMakeFiles/l3binmerge.dir/l3binmerge.cpp.o.provides: src/l3bin/CMakeFiles/l3binmerge.dir/l3binmerge.cpp.o.requires
	$(MAKE) -f src/l3bin/CMakeFiles/l3binmerge.dir/build.make src/l3bin/CMakeFiles/l3binmerge.dir/l3binmerge.cpp.o.provides.build
.PHONY : src/l3bin/CMakeFiles/l3binmerge.dir/l3binmerge.cpp.o.provides

src/l3bin/CMakeFiles/l3binmerge.dir/l3binmerge.cpp.o.provides.build: src/l3bin/CMakeFiles/l3binmerge.dir/l3binmerge.cpp.o

src/l3bin/CMakeFiles/l3binmerge.dir/l3bin_input.c.o: src/l3bin/CMakeFiles/l3binmerge.dir/flags.make
src/l3bin/CMakeFiles/l3binmerge.dir/l3bin_input.c.o: ../src/l3bin/l3bin_input.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/jaemoo/seadas-7.4/ocssw/build/cbuild/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object src/l3bin/CMakeFiles/l3binmerge.dir/l3bin_input.c.o"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l3bin && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/l3binmerge.dir/l3bin_input.c.o   -c /home/jaemoo/seadas-7.4/ocssw/build/src/l3bin/l3bin_input.c

src/l3bin/CMakeFiles/l3binmerge.dir/l3bin_input.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/l3binmerge.dir/l3bin_input.c.i"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l3bin && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/jaemoo/seadas-7.4/ocssw/build/src/l3bin/l3bin_input.c > CMakeFiles/l3binmerge.dir/l3bin_input.c.i

src/l3bin/CMakeFiles/l3binmerge.dir/l3bin_input.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/l3binmerge.dir/l3bin_input.c.s"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l3bin && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/jaemoo/seadas-7.4/ocssw/build/src/l3bin/l3bin_input.c -o CMakeFiles/l3binmerge.dir/l3bin_input.c.s

src/l3bin/CMakeFiles/l3binmerge.dir/l3bin_input.c.o.requires:
.PHONY : src/l3bin/CMakeFiles/l3binmerge.dir/l3bin_input.c.o.requires

src/l3bin/CMakeFiles/l3binmerge.dir/l3bin_input.c.o.provides: src/l3bin/CMakeFiles/l3binmerge.dir/l3bin_input.c.o.requires
	$(MAKE) -f src/l3bin/CMakeFiles/l3binmerge.dir/build.make src/l3bin/CMakeFiles/l3binmerge.dir/l3bin_input.c.o.provides.build
.PHONY : src/l3bin/CMakeFiles/l3binmerge.dir/l3bin_input.c.o.provides

src/l3bin/CMakeFiles/l3binmerge.dir/l3bin_input.c.o.provides.build: src/l3bin/CMakeFiles/l3binmerge.dir/l3bin_input.c.o

# Object files for target l3binmerge
l3binmerge_OBJECTS = \
"CMakeFiles/l3binmerge.dir/l3binmerge.cpp.o" \
"CMakeFiles/l3binmerge.dir/l3bin_input.c.o"

# External object files for target l3binmerge
l3binmerge_EXTERNAL_OBJECTS =

src/l3bin/l3binmerge: src/l3bin/CMakeFiles/l3binmerge.dir/l3binmerge.cpp.o
src/l3bin/l3binmerge: src/l3bin/CMakeFiles/l3binmerge.dir/l3bin_input.c.o
src/l3bin/l3binmerge: src/l3bin/CMakeFiles/l3binmerge.dir/build.make
src/l3bin/l3binmerge: src/libbin++/libbin++.a
src/l3bin/l3binmerge: src/libbin/libbin.a
src/l3bin/l3binmerge: src/libgenutils/libgenutils.a
src/l3bin/l3binmerge: src/libdfutils/libdfutils.a
src/l3bin/l3binmerge: src/libhdf4utils/libhdf4utils.a
src/l3bin/l3binmerge: src/libnetcdfutils/libnetcdfutils.a
src/l3bin/l3binmerge: src/libhdf5utils/libhdf5utils.a
src/l3bin/l3binmerge: src/libpiutils/libpiutils.a
src/l3bin/l3binmerge: src/libgenutils/libgenutils.a
src/l3bin/l3binmerge: src/libtimeutils/libtimeutils.a
src/l3bin/l3binmerge: src/l3bin/CMakeFiles/l3binmerge.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable l3binmerge"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l3bin && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/l3binmerge.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/l3bin/CMakeFiles/l3binmerge.dir/build: src/l3bin/l3binmerge
.PHONY : src/l3bin/CMakeFiles/l3binmerge.dir/build

src/l3bin/CMakeFiles/l3binmerge.dir/requires: src/l3bin/CMakeFiles/l3binmerge.dir/l3binmerge.cpp.o.requires
src/l3bin/CMakeFiles/l3binmerge.dir/requires: src/l3bin/CMakeFiles/l3binmerge.dir/l3bin_input.c.o.requires
.PHONY : src/l3bin/CMakeFiles/l3binmerge.dir/requires

src/l3bin/CMakeFiles/l3binmerge.dir/clean:
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l3bin && $(CMAKE_COMMAND) -P CMakeFiles/l3binmerge.dir/cmake_clean.cmake
.PHONY : src/l3bin/CMakeFiles/l3binmerge.dir/clean

src/l3bin/CMakeFiles/l3binmerge.dir/depend:
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jaemoo/seadas-7.4/ocssw/build /home/jaemoo/seadas-7.4/ocssw/build/src/l3bin /home/jaemoo/seadas-7.4/ocssw/build/cbuild /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l3bin /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l3bin/CMakeFiles/l3binmerge.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/l3bin/CMakeFiles/l3binmerge.dir/depend
