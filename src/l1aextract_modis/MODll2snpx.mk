# ===========================================================================
#       Begin Program Specific Block
# ===========================================================================

# Makefile name
MAKEFILE = MODll2snpx.mk

# Progam to make
EXE     = MODll2snpx

# Object modules for EXE
OBJ     = MODll2snpx.o

# Include file locations
INCLUDE	= -I$(LIB3_INC)

# LIBS needed to compile EXE
LIBS    = -L$(LIB3_LIB)

# Libraries required to link
LLIBS 	= -lmfhdf -ldf -ljpeg -lz

# Program specific additions to compiler flags
LOCAL_CFLAGS =

include $(MAKEFILE_APP_TEMPLATE)

