# ===========================================================================
#       Begin Program Specific Block
# ===========================================================================

# Makefile name
MAKEFILE= aq_sca_main.mk

# Progam to make
EXE     = aq_sca_main

# Object modules for EXE

OBJ     = \
	aq_sca_main.o \
	tb2vsm.o

# Include file locations
INCLUDE = -I$(LIB3_INC)

# Library locations
LIBS 	= -L$(LIB3_LIB)

# Libraries required to link
LLIBS 	= -lhdf5_fortran -lhdf5 -lz

LOCAL_FFLAGS = -O0

include $(MAKEFILE_APP_TEMPLATE)
