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
include src/libhdf4utils/CMakeFiles/hdf4utils.dir/depend.make

# Include the progress variables for this target.
include src/libhdf4utils/CMakeFiles/hdf4utils.dir/progress.make

# Include the compile flags for this target's objects.
include src/libhdf4utils/CMakeFiles/hdf4utils.dir/flags.make

src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf_utils.c.o: src/libhdf4utils/CMakeFiles/hdf4utils.dir/flags.make
src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf_utils.c.o: ../src/libhdf4utils/hdf_utils.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/jaemoo/seadas-7.4/ocssw/build/cbuild/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf_utils.c.o"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libhdf4utils && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/hdf4utils.dir/hdf_utils.c.o   -c /home/jaemoo/seadas-7.4/ocssw/build/src/libhdf4utils/hdf_utils.c

src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf_utils.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/hdf4utils.dir/hdf_utils.c.i"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libhdf4utils && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/jaemoo/seadas-7.4/ocssw/build/src/libhdf4utils/hdf_utils.c > CMakeFiles/hdf4utils.dir/hdf_utils.c.i

src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf_utils.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/hdf4utils.dir/hdf_utils.c.s"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libhdf4utils && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/jaemoo/seadas-7.4/ocssw/build/src/libhdf4utils/hdf_utils.c -o CMakeFiles/hdf4utils.dir/hdf_utils.c.s

src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf_utils.c.o.requires:
.PHONY : src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf_utils.c.o.requires

src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf_utils.c.o.provides: src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf_utils.c.o.requires
	$(MAKE) -f src/libhdf4utils/CMakeFiles/hdf4utils.dir/build.make src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf_utils.c.o.provides.build
.PHONY : src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf_utils.c.o.provides

src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf_utils.c.o.provides.build: src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf_utils.c.o

src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf4_utils.c.o: src/libhdf4utils/CMakeFiles/hdf4utils.dir/flags.make
src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf4_utils.c.o: ../src/libhdf4utils/hdf4_utils.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/jaemoo/seadas-7.4/ocssw/build/cbuild/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf4_utils.c.o"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libhdf4utils && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/hdf4utils.dir/hdf4_utils.c.o   -c /home/jaemoo/seadas-7.4/ocssw/build/src/libhdf4utils/hdf4_utils.c

src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf4_utils.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/hdf4utils.dir/hdf4_utils.c.i"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libhdf4utils && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/jaemoo/seadas-7.4/ocssw/build/src/libhdf4utils/hdf4_utils.c > CMakeFiles/hdf4utils.dir/hdf4_utils.c.i

src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf4_utils.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/hdf4utils.dir/hdf4_utils.c.s"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libhdf4utils && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/jaemoo/seadas-7.4/ocssw/build/src/libhdf4utils/hdf4_utils.c -o CMakeFiles/hdf4utils.dir/hdf4_utils.c.s

src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf4_utils.c.o.requires:
.PHONY : src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf4_utils.c.o.requires

src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf4_utils.c.o.provides: src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf4_utils.c.o.requires
	$(MAKE) -f src/libhdf4utils/CMakeFiles/hdf4utils.dir/build.make src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf4_utils.c.o.provides.build
.PHONY : src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf4_utils.c.o.provides

src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf4_utils.c.o.provides.build: src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf4_utils.c.o

# Object files for target hdf4utils
hdf4utils_OBJECTS = \
"CMakeFiles/hdf4utils.dir/hdf_utils.c.o" \
"CMakeFiles/hdf4utils.dir/hdf4_utils.c.o"

# External object files for target hdf4utils
hdf4utils_EXTERNAL_OBJECTS =

src/libhdf4utils/libhdf4utils.a: src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf_utils.c.o
src/libhdf4utils/libhdf4utils.a: src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf4_utils.c.o
src/libhdf4utils/libhdf4utils.a: src/libhdf4utils/CMakeFiles/hdf4utils.dir/build.make
src/libhdf4utils/libhdf4utils.a: src/libhdf4utils/CMakeFiles/hdf4utils.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libhdf4utils.a"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libhdf4utils && $(CMAKE_COMMAND) -P CMakeFiles/hdf4utils.dir/cmake_clean_target.cmake
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libhdf4utils && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/hdf4utils.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/libhdf4utils/CMakeFiles/hdf4utils.dir/build: src/libhdf4utils/libhdf4utils.a
.PHONY : src/libhdf4utils/CMakeFiles/hdf4utils.dir/build

src/libhdf4utils/CMakeFiles/hdf4utils.dir/requires: src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf_utils.c.o.requires
src/libhdf4utils/CMakeFiles/hdf4utils.dir/requires: src/libhdf4utils/CMakeFiles/hdf4utils.dir/hdf4_utils.c.o.requires
.PHONY : src/libhdf4utils/CMakeFiles/hdf4utils.dir/requires

src/libhdf4utils/CMakeFiles/hdf4utils.dir/clean:
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libhdf4utils && $(CMAKE_COMMAND) -P CMakeFiles/hdf4utils.dir/cmake_clean.cmake
.PHONY : src/libhdf4utils/CMakeFiles/hdf4utils.dir/clean

src/libhdf4utils/CMakeFiles/hdf4utils.dir/depend:
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jaemoo/seadas-7.4/ocssw/build /home/jaemoo/seadas-7.4/ocssw/build/src/libhdf4utils /home/jaemoo/seadas-7.4/ocssw/build/cbuild /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libhdf4utils /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libhdf4utils/CMakeFiles/hdf4utils.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/libhdf4utils/CMakeFiles/hdf4utils.dir/depend
