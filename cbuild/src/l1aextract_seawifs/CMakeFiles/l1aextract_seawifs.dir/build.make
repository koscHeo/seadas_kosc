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
include src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/depend.make

# Include the progress variables for this target.
include src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/progress.make

# Include the compile flags for this target's objects.
include src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/flags.make

src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/main_l1aextract.c.o: src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/flags.make
src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/main_l1aextract.c.o: ../src/l1aextract_seawifs/main_l1aextract.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/jaemoo/seadas-7.4/ocssw/build/cbuild/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/main_l1aextract.c.o"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l1aextract_seawifs && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/l1aextract_seawifs.dir/main_l1aextract.c.o   -c /home/jaemoo/seadas-7.4/ocssw/build/src/l1aextract_seawifs/main_l1aextract.c

src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/main_l1aextract.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/l1aextract_seawifs.dir/main_l1aextract.c.i"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l1aextract_seawifs && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/jaemoo/seadas-7.4/ocssw/build/src/l1aextract_seawifs/main_l1aextract.c > CMakeFiles/l1aextract_seawifs.dir/main_l1aextract.c.i

src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/main_l1aextract.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/l1aextract_seawifs.dir/main_l1aextract.c.s"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l1aextract_seawifs && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/jaemoo/seadas-7.4/ocssw/build/src/l1aextract_seawifs/main_l1aextract.c -o CMakeFiles/l1aextract_seawifs.dir/main_l1aextract.c.s

src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/main_l1aextract.c.o.requires:
.PHONY : src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/main_l1aextract.c.o.requires

src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/main_l1aextract.c.o.provides: src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/main_l1aextract.c.o.requires
	$(MAKE) -f src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/build.make src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/main_l1aextract.c.o.provides.build
.PHONY : src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/main_l1aextract.c.o.provides

src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/main_l1aextract.c.o.provides.build: src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/main_l1aextract.c.o

src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/extract_sub.c.o: src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/flags.make
src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/extract_sub.c.o: ../src/l1aextract_seawifs/extract_sub.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/jaemoo/seadas-7.4/ocssw/build/cbuild/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/extract_sub.c.o"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l1aextract_seawifs && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/l1aextract_seawifs.dir/extract_sub.c.o   -c /home/jaemoo/seadas-7.4/ocssw/build/src/l1aextract_seawifs/extract_sub.c

src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/extract_sub.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/l1aextract_seawifs.dir/extract_sub.c.i"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l1aextract_seawifs && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/jaemoo/seadas-7.4/ocssw/build/src/l1aextract_seawifs/extract_sub.c > CMakeFiles/l1aextract_seawifs.dir/extract_sub.c.i

src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/extract_sub.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/l1aextract_seawifs.dir/extract_sub.c.s"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l1aextract_seawifs && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/jaemoo/seadas-7.4/ocssw/build/src/l1aextract_seawifs/extract_sub.c -o CMakeFiles/l1aextract_seawifs.dir/extract_sub.c.s

src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/extract_sub.c.o.requires:
.PHONY : src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/extract_sub.c.o.requires

src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/extract_sub.c.o.provides: src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/extract_sub.c.o.requires
	$(MAKE) -f src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/build.make src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/extract_sub.c.o.provides.build
.PHONY : src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/extract_sub.c.o.provides

src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/extract_sub.c.o.provides.build: src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/extract_sub.c.o

# Object files for target l1aextract_seawifs
l1aextract_seawifs_OBJECTS = \
"CMakeFiles/l1aextract_seawifs.dir/main_l1aextract.c.o" \
"CMakeFiles/l1aextract_seawifs.dir/extract_sub.c.o"

# External object files for target l1aextract_seawifs
l1aextract_seawifs_EXTERNAL_OBJECTS =

src/l1aextract_seawifs/l1aextract_seawifs: src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/main_l1aextract.c.o
src/l1aextract_seawifs/l1aextract_seawifs: src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/extract_sub.c.o
src/l1aextract_seawifs/l1aextract_seawifs: src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/build.make
src/l1aextract_seawifs/l1aextract_seawifs: src/libseawifs/libseawifs.a
src/l1aextract_seawifs/l1aextract_seawifs: src/libhdf4utils/libhdf4utils.a
src/l1aextract_seawifs/l1aextract_seawifs: src/libpiutils/libpiutils.a
src/l1aextract_seawifs/l1aextract_seawifs: src/libgenutils/libgenutils.a
src/l1aextract_seawifs/l1aextract_seawifs: src/libtimeutils/libtimeutils.a
src/l1aextract_seawifs/l1aextract_seawifs: src/libnav/libnav.a
src/l1aextract_seawifs/l1aextract_seawifs: src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable l1aextract_seawifs"
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l1aextract_seawifs && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/l1aextract_seawifs.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/build: src/l1aextract_seawifs/l1aextract_seawifs
.PHONY : src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/build

src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/requires: src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/main_l1aextract.c.o.requires
src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/requires: src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/extract_sub.c.o.requires
.PHONY : src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/requires

src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/clean:
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l1aextract_seawifs && $(CMAKE_COMMAND) -P CMakeFiles/l1aextract_seawifs.dir/cmake_clean.cmake
.PHONY : src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/clean

src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/depend:
	cd /home/jaemoo/seadas-7.4/ocssw/build/cbuild && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jaemoo/seadas-7.4/ocssw/build /home/jaemoo/seadas-7.4/ocssw/build/src/l1aextract_seawifs /home/jaemoo/seadas-7.4/ocssw/build/cbuild /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l1aextract_seawifs /home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/l1aextract_seawifs/CMakeFiles/l1aextract_seawifs.dir/depend
