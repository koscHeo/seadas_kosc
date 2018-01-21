#! /bin/sh
#
# Copyright by The HDF Group.
# Copyright by the Board of Trustees of the University of Illinois.
# All rights reserved.
#
# This file is part of HDF5.  The full HDF5 copyright notice, including
# terms governing use, modification, and redistribution, is contained in
# the files COPYING and Copyright.html.  COPYING can be found at the root
# of the source code distribution tree; Copyright.html can be found at the
# root level of an installed copy of the electronic HDF5 document set and
# is linked from the top-level documents page.  It can also be found at
# http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have
# access to either file, you may request a copy from help@hdfgroup.org.
#

# Scripts for running testphdf5 program with a variety of parameters
top_srcdir=..
top_builddir=..
srcdir=.

TEST_APP=testphdf5               # The tool name
TEST_APP_BIN=`pwd`/$TEST_APP    # The path of the tool binary

nerrors=0
verbose=yes

# Print a line-line message left justified in a field of 70 characters
# beginning with the word "Testing".
#
TESTING() {
   SPACES="                                                               "
   echo "Testing $* $SPACES" | cut -c1-70 | tr -d '\012'
}

# Run a test.  If a test fails then increment the `nerrors' global variable.
#
TOOLTEST() {
    # Run test.
    echo $RUNPARALLEL $TEST_APP_BIN "$@"
    eval $RUNPARALLEL $TEST_APP_BIN "$@"

    # Check if the command failed and increment nerrors if so.
    if test $? -ne 0 ; then
        nerrors="`expr $nerrors + 1`"
    fi
}

##############################################################################
##############################################################################
###			  T H E   T E S T S                                ###
##############################################################################
##############################################################################

# testphdf5 test using the MPI-POSIX VFL driver
TOOLTEST -p

# Temporary patch:
# Run t_shapesame this way. Need more permanent solution.
echo $RUNPARALLEL ./t_shapesame -p
eval $RUNPARALLEL ./t_shapesame -p

# Check if the command failed and increment nerrors if so.
if test $? -ne 0 ; then
    nerrors="`expr $nerrors + 1`"
fi
# Temporary patch ended.

# Emit message about testing status
if test $nerrors -eq 0 ; then
   echo "All $TEST_APP tests passed."
else
   echo "ERROR! One or more $TEST_APP tests failed."
fi

# Propagate a useful exit code
exit $nerrors
