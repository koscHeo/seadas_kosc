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
include src/l2bin/CMakeFiles/l2bin.dir/depend.make

# Include the progress variables for this target.
include src/l2bin/CMakeFiles/l2bin.dir/progress.make

# Include the compile flags for this target's objects.
include src/l2bin/CMakeFiles/l2bin.dir/flags.make

src/l2bin/CMakeFiles/l2bin.dir/l2bin.c.o: src/l2bin/CMakeFiles/l2bin.dir/flags.make
src/l2bin/CMakeFiles/l2bin.dir/l2bin.c.o: ../src/l2bin/l2bin.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/jaemoo/seadas-7.4/ocssw/build/cbuild/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object src/l2bin/CMakeFiles/l2bin.dir/l2bin.c.o"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l2bin && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/l2bin.dir/l2bin.c.o   -c /home/jaemoo/seadas-7.4/ocssw/build/src/l2bin/l2bin.c

src/l2bin/CMakeFiles/l2bin.dir/l2bin.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/l2bin.dir/l2bin.c.i"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l2bin && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/jaemoo/seadas-7.4/ocssw/build/src/l2bin/l2bin.c > CMakeFiles/l2bin.dir/l2bin.c.i

src/l2bin/CMakeFiles/l2bin.dir/l2bin.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/l2bin.dir/l2bin.c.s"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l2bin && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/jaemoo/seadas-7.4/ocssw/build/src/l2bin/l2bin.c -o CMakeFiles/l2bin.dir/l2bin.c.s

src/l2bin/CMakeFiles/l2bin.dir/l2bin.c.o.requires:
.PHONY : src/l2bin/CMakeFiles/l2bin.dir/l2bin.c.o.requires

src/l2bin/CMakeFiles/l2bin.dir/l2bin.c.o.provides: src/l2bin/CMakeFiles/l2bin.dir/l2bin.c.o.requires
	$(MAKE) -f src/l2bin/CMakeFiles/l2bin.dir/build.make src/l2bin/CMakeFiles/l2bin.dir/l2bin.c.o.provides.build
.PHONY : src/l2bin/CMakeFiles/l2bin.dir/l2bin.c.o.provides

src/l2bin/CMakeFiles/l2bin.dir/l2bin.c.o.provides.build: src/l2bin/CMakeFiles/l2bin.dir/l2bin.c.o

src/l2bin/CMakeFiles/l2bin.dir/l2bin_input.c.o: src/l2bin/CMakeFiles/l2bin.dir/flags.make
src/l2bin/CMakeFiles/l2bin.dir/l2bin_input.c.o: ../src/l2bin/l2bin_input.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/jaemoo/seadas-7.4/ocssw/build/cbuild/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object src/l2bin/CMakeFiles/l2bin.dir/l2bin_input.c.o"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l2bin && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/l2bin.dir/l2bin_input.c.o   -c /home/jaemoo/seadas-7.4/ocssw/build/src/l2bin/l2bin_input.c

src/l2bin/CMakeFiles/l2bin.dir/l2bin_input.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/l2bin.dir/l2bin_input.c.i"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l2bin && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/jaemoo/seadas-7.4/ocssw/build/src/l2bin/l2bin_input.c > CMakeFiles/l2bin.dir/l2bin_input.c.i

src/l2bin/CMakeFiles/l2bin.dir/l2bin_input.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/l2bin.dir/l2bin_input.c.s"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l2bin && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/jaemoo/seadas-7.4/ocssw/build/src/l2bin/l2bin_input.c -o CMakeFiles/l2bin.dir/l2bin_input.c.s

src/l2bin/CMakeFiles/l2bin.dir/l2bin_input.c.o.requires:
.PHONY : src/l2bin/CMakeFiles/l2bin.dir/l2bin_input.c.o.requires

src/l2bin/CMakeFiles/l2bin.dir/l2bin_input.c.o.provides: src/l2bin/CMakeFiles/l2bin.dir/l2bin_input.c.o.requires
	$(MAKE) -f src/l2bin/CMakeFiles/l2bin.dir/build.make src/l2bin/CMakeFiles/l2bin.dir/l2bin_input.c.o.provides.build
.PHONY : src/l2bin/CMakeFiles/l2bin.dir/l2bin_input.c.o.provides

src/l2bin/CMakeFiles/l2bin.dir/l2bin_input.c.o.provides.build: src/l2bin/CMakeFiles/l2bin.dir/l2bin_input.c.o

src/l2bin/CMakeFiles/l2bin.dir/query_disc_wrapper.f90.o: src/l2bin/CMakeFiles/l2bin.dir/flags.make
src/l2bin/CMakeFiles/l2bin.dir/query_disc_wrapper.f90.o: ../src/l2bin/query_disc_wrapper.f90
	$(CMAKE_COMMAND) -E cmake_progress_report /home/jaemoo/seadas-7.4/ocssw/build/cbuild/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building Fortran object src/l2bin/CMakeFiles/l2bin.dir/query_disc_wrapper.f90.o"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l2bin && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_FLAGS) -c /home/jaemoo/seadas-7.4/ocssw/build/src/l2bin/query_disc_wrapper.f90 -o CMakeFiles/l2bin.dir/query_disc_wrapper.f90.o

src/l2bin/CMakeFiles/l2bin.dir/query_disc_wrapper.f90.o.requires:
.PHONY : src/l2bin/CMakeFiles/l2bin.dir/query_disc_wrapper.f90.o.requires

src/l2bin/CMakeFiles/l2bin.dir/query_disc_wrapper.f90.o.provides: src/l2bin/CMakeFiles/l2bin.dir/query_disc_wrapper.f90.o.requires
	$(MAKE) -f src/l2bin/CMakeFiles/l2bin.dir/build.make src/l2bin/CMakeFiles/l2bin.dir/query_disc_wrapper.f90.o.provides.build
.PHONY : src/l2bin/CMakeFiles/l2bin.dir/query_disc_wrapper.f90.o.provides

src/l2bin/CMakeFiles/l2bin.dir/query_disc_wrapper.f90.o.provides.build: src/l2bin/CMakeFiles/l2bin.dir/query_disc_wrapper.f90.o

src/l2bin/CMakeFiles/l2bin.dir/dataday.c.o: src/l2bin/CMakeFiles/l2bin.dir/flags.make
src/l2bin/CMakeFiles/l2bin.dir/dataday.c.o: ../src/l2bin/dataday.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/jaemoo/seadas-7.4/ocssw/build/cbuild/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object src/l2bin/CMakeFiles/l2bin.dir/dataday.c.o"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l2bin && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/l2bin.dir/dataday.c.o   -c /home/jaemoo/seadas-7.4/ocssw/build/src/l2bin/dataday.c

src/l2bin/CMakeFiles/l2bin.dir/dataday.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/l2bin.dir/dataday.c.i"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l2bin && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/jaemoo/seadas-7.4/ocssw/build/src/l2bin/dataday.c > CMakeFiles/l2bin.dir/dataday.c.i

src/l2bin/CMakeFiles/l2bin.dir/dataday.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/l2bin.dir/dataday.c.s"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l2bin && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/jaemoo/seadas-7.4/ocssw/build/src/l2bin/dataday.c -o CMakeFiles/l2bin.dir/dataday.c.s

src/l2bin/CMakeFiles/l2bin.dir/dataday.c.o.requires:
.PHONY : src/l2bin/CMakeFiles/l2bin.dir/dataday.c.o.requires

src/l2bin/CMakeFiles/l2bin.dir/dataday.c.o.provides: src/l2bin/CMakeFiles/l2bin.dir/dataday.c.o.requires
	$(MAKE) -f src/l2bin/CMakeFiles/l2bin.dir/build.make src/l2bin/CMakeFiles/l2bin.dir/dataday.c.o.provides.build
.PHONY : src/l2bin/CMakeFiles/l2bin.dir/dataday.c.o.provides

src/l2bin/CMakeFiles/l2bin.dir/dataday.c.o.provides.build: src/l2bin/CMakeFiles/l2bin.dir/dataday.c.o

# Object files for target l2bin
l2bin_OBJECTS = \
"CMakeFiles/l2bin.dir/l2bin.c.o" \
"CMakeFiles/l2bin.dir/l2bin_input.c.o" \
"CMakeFiles/l2bin.dir/query_disc_wrapper.f90.o" \
"CMakeFiles/l2bin.dir/dataday.c.o"

# External object files for target l2bin
l2bin_EXTERNAL_OBJECTS =

src/l2bin/l2bin: src/l2bin/CMakeFiles/l2bin.dir/l2bin.c.o
src/l2bin/l2bin: src/l2bin/CMakeFiles/l2bin.dir/l2bin_input.c.o
src/l2bin/l2bin: src/l2bin/CMakeFiles/l2bin.dir/query_disc_wrapper.f90.o
src/l2bin/l2bin: src/l2bin/CMakeFiles/l2bin.dir/dataday.c.o
src/l2bin/l2bin: src/l2bin/CMakeFiles/l2bin.dir/build.make
src/l2bin/l2bin: src/libl2/libl2.a
src/l2bin/l2bin: src/libbin/libbin.a
src/l2bin/l2bin: src/libnav/libnav.a
src/l2bin/l2bin: src/libseawifs/libseawifs.a
src/l2bin/l2bin: src/libnav/libnav.a
src/l2bin/l2bin: src/libdfutils/libdfutils.a
src/l2bin/l2bin: src/libhdf4utils/libhdf4utils.a
src/l2bin/l2bin: src/libnetcdfutils/libnetcdfutils.a
src/l2bin/l2bin: src/libhdf5utils/libhdf5utils.a
src/l2bin/l2bin: src/libpiutils/libpiutils.a
src/l2bin/l2bin: src/libgenutils/libgenutils.a
src/l2bin/l2bin: src/libtimeutils/libtimeutils.a
src/l2bin/l2bin: src/l2bin/CMakeFiles/l2bin.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking Fortran executable l2bin"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l2bin && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/l2bin.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/l2bin/CMakeFiles/l2bin.dir/build: src/l2bin/l2bin
.PHONY : src/l2bin/CMakeFiles/l2bin.dir/build

src/l2bin/CMakeFiles/l2bin.dir/requires: src/l2bin/CMakeFiles/l2bin.dir/l2bin.c.o.requires
src/l2bin/CMakeFiles/l2bin.dir/requires: src/l2bin/CMakeFiles/l2bin.dir/l2bin_input.c.o.requires
src/l2bin/CMakeFiles/l2bin.dir/requires: src/l2bin/CMakeFiles/l2bin.dir/query_disc_wrapper.f90.o.requires
src/l2bin/CMakeFiles/l2bin.dir/requires: src/l2bin/CMakeFiles/l2bin.dir/dataday.c.o.requires
.PHONY : src/l2bin/CMakeFiles/l2bin.dir/requires

src/l2bin/CMakeFiles/l2bin.dir/clean:
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l2bin && $(CMAKE_COMMAND) -P CMakeFiles/l2bin.dir/cmake_clean.cmake
.PHONY : src/l2bin/CMakeFiles/l2bin.dir/clean

src/l2bin/CMakeFiles/l2bin.dir/depend:
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jaemoo/seadas-7.4/ocssw/build /home/jaemoo/seadas-7.4/ocssw/build/src/l2bin /home/jaemoo/seadas-7.4/ocssw/build/cbuild /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l2bin /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l2bin/CMakeFiles/l2bin.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/l2bin/CMakeFiles/l2bin.dir/depend

