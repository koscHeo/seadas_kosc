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
include src/smigen/CMakeFiles/smigen.dir/depend.make

# Include the progress variables for this target.
include src/smigen/CMakeFiles/smigen.dir/progress.make

# Include the compile flags for this target's objects.
include src/smigen/CMakeFiles/smigen.dir/flags.make

src/smigen/CMakeFiles/smigen.dir/smigen.cpp.o: src/smigen/CMakeFiles/smigen.dir/flags.make
src/smigen/CMakeFiles/smigen.dir/smigen.cpp.o: ../src/smigen/smigen.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/jaemoo/seadas-7.4/ocssw/build/cbuild/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/smigen/CMakeFiles/smigen.dir/smigen.cpp.o"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/smigen && /usr/bin/g++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/smigen.dir/smigen.cpp.o -c /home/jaemoo/seadas-7.4/ocssw/build/src/smigen/smigen.cpp

src/smigen/CMakeFiles/smigen.dir/smigen.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/smigen.dir/smigen.cpp.i"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/smigen && /usr/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/jaemoo/seadas-7.4/ocssw/build/src/smigen/smigen.cpp > CMakeFiles/smigen.dir/smigen.cpp.i

src/smigen/CMakeFiles/smigen.dir/smigen.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/smigen.dir/smigen.cpp.s"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/smigen && /usr/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/jaemoo/seadas-7.4/ocssw/build/src/smigen/smigen.cpp -o CMakeFiles/smigen.dir/smigen.cpp.s

src/smigen/CMakeFiles/smigen.dir/smigen.cpp.o.requires:
.PHONY : src/smigen/CMakeFiles/smigen.dir/smigen.cpp.o.requires

src/smigen/CMakeFiles/smigen.dir/smigen.cpp.o.provides: src/smigen/CMakeFiles/smigen.dir/smigen.cpp.o.requires
	$(MAKE) -f src/smigen/CMakeFiles/smigen.dir/build.make src/smigen/CMakeFiles/smigen.dir/smigen.cpp.o.provides.build
.PHONY : src/smigen/CMakeFiles/smigen.dir/smigen.cpp.o.provides

src/smigen/CMakeFiles/smigen.dir/smigen.cpp.o.provides.build: src/smigen/CMakeFiles/smigen.dir/smigen.cpp.o

src/smigen/CMakeFiles/smigen.dir/smigen_input.c.o: src/smigen/CMakeFiles/smigen.dir/flags.make
src/smigen/CMakeFiles/smigen.dir/smigen_input.c.o: ../src/smigen/smigen_input.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/jaemoo/seadas-7.4/ocssw/build/cbuild/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object src/smigen/CMakeFiles/smigen.dir/smigen_input.c.o"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/smigen && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/smigen.dir/smigen_input.c.o   -c /home/jaemoo/seadas-7.4/ocssw/build/src/smigen/smigen_input.c

src/smigen/CMakeFiles/smigen.dir/smigen_input.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/smigen.dir/smigen_input.c.i"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/smigen && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/jaemoo/seadas-7.4/ocssw/build/src/smigen/smigen_input.c > CMakeFiles/smigen.dir/smigen_input.c.i

src/smigen/CMakeFiles/smigen.dir/smigen_input.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/smigen.dir/smigen_input.c.s"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/smigen && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/jaemoo/seadas-7.4/ocssw/build/src/smigen/smigen_input.c -o CMakeFiles/smigen.dir/smigen_input.c.s

src/smigen/CMakeFiles/smigen.dir/smigen_input.c.o.requires:
.PHONY : src/smigen/CMakeFiles/smigen.dir/smigen_input.c.o.requires

src/smigen/CMakeFiles/smigen.dir/smigen_input.c.o.provides: src/smigen/CMakeFiles/smigen.dir/smigen_input.c.o.requires
	$(MAKE) -f src/smigen/CMakeFiles/smigen.dir/build.make src/smigen/CMakeFiles/smigen.dir/smigen_input.c.o.provides.build
.PHONY : src/smigen/CMakeFiles/smigen.dir/smigen_input.c.o.provides

src/smigen/CMakeFiles/smigen.dir/smigen_input.c.o.provides.build: src/smigen/CMakeFiles/smigen.dir/smigen_input.c.o

src/smigen/CMakeFiles/smigen.dir/put_smi.cpp.o: src/smigen/CMakeFiles/smigen.dir/flags.make
src/smigen/CMakeFiles/smigen.dir/put_smi.cpp.o: ../src/smigen/put_smi.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/jaemoo/seadas-7.4/ocssw/build/cbuild/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/smigen/CMakeFiles/smigen.dir/put_smi.cpp.o"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/smigen && /usr/bin/g++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/smigen.dir/put_smi.cpp.o -c /home/jaemoo/seadas-7.4/ocssw/build/src/smigen/put_smi.cpp

src/smigen/CMakeFiles/smigen.dir/put_smi.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/smigen.dir/put_smi.cpp.i"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/smigen && /usr/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/jaemoo/seadas-7.4/ocssw/build/src/smigen/put_smi.cpp > CMakeFiles/smigen.dir/put_smi.cpp.i

src/smigen/CMakeFiles/smigen.dir/put_smi.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/smigen.dir/put_smi.cpp.s"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/smigen && /usr/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/jaemoo/seadas-7.4/ocssw/build/src/smigen/put_smi.cpp -o CMakeFiles/smigen.dir/put_smi.cpp.s

src/smigen/CMakeFiles/smigen.dir/put_smi.cpp.o.requires:
.PHONY : src/smigen/CMakeFiles/smigen.dir/put_smi.cpp.o.requires

src/smigen/CMakeFiles/smigen.dir/put_smi.cpp.o.provides: src/smigen/CMakeFiles/smigen.dir/put_smi.cpp.o.requires
	$(MAKE) -f src/smigen/CMakeFiles/smigen.dir/build.make src/smigen/CMakeFiles/smigen.dir/put_smi.cpp.o.provides.build
.PHONY : src/smigen/CMakeFiles/smigen.dir/put_smi.cpp.o.provides

src/smigen/CMakeFiles/smigen.dir/put_smi.cpp.o.provides.build: src/smigen/CMakeFiles/smigen.dir/put_smi.cpp.o

# Object files for target smigen
smigen_OBJECTS = \
"CMakeFiles/smigen.dir/smigen.cpp.o" \
"CMakeFiles/smigen.dir/smigen_input.c.o" \
"CMakeFiles/smigen.dir/put_smi.cpp.o"

# External object files for target smigen
smigen_EXTERNAL_OBJECTS =

src/smigen/smigen: src/smigen/CMakeFiles/smigen.dir/smigen.cpp.o
src/smigen/smigen: src/smigen/CMakeFiles/smigen.dir/smigen_input.c.o
src/smigen/smigen: src/smigen/CMakeFiles/smigen.dir/put_smi.cpp.o
src/smigen/smigen: src/smigen/CMakeFiles/smigen.dir/build.make
src/smigen/smigen: src/libbin++/libbin++.a
src/smigen/smigen: src/libbin/libbin.a
src/smigen/smigen: src/libdfutils/libdfutils.a
src/smigen/smigen: src/libhdf4utils/libhdf4utils.a
src/smigen/smigen: src/libnetcdfutils/libnetcdfutils.a
src/smigen/smigen: src/libhdf5utils/libhdf5utils.a
src/smigen/smigen: src/libpiutils/libpiutils.a
src/smigen/smigen: src/libgenutils/libgenutils.a
src/smigen/smigen: src/libtimeutils/libtimeutils.a
src/smigen/smigen: src/smigen/CMakeFiles/smigen.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable smigen"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/smigen && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/smigen.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/smigen/CMakeFiles/smigen.dir/build: src/smigen/smigen
.PHONY : src/smigen/CMakeFiles/smigen.dir/build

src/smigen/CMakeFiles/smigen.dir/requires: src/smigen/CMakeFiles/smigen.dir/smigen.cpp.o.requires
src/smigen/CMakeFiles/smigen.dir/requires: src/smigen/CMakeFiles/smigen.dir/smigen_input.c.o.requires
src/smigen/CMakeFiles/smigen.dir/requires: src/smigen/CMakeFiles/smigen.dir/put_smi.cpp.o.requires
.PHONY : src/smigen/CMakeFiles/smigen.dir/requires

src/smigen/CMakeFiles/smigen.dir/clean:
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/smigen && $(CMAKE_COMMAND) -P CMakeFiles/smigen.dir/cmake_clean.cmake
.PHONY : src/smigen/CMakeFiles/smigen.dir/clean

src/smigen/CMakeFiles/smigen.dir/depend:
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jaemoo/seadas-7.4/ocssw/build /home/jaemoo/seadas-7.4/ocssw/build/src/smigen /home/jaemoo/seadas-7.4/ocssw/build/cbuild /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/smigen /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/smigen/CMakeFiles/smigen.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/smigen/CMakeFiles/smigen.dir/depend

