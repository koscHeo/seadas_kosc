# ===========================================================================
#       Begin Program Specific Block
# ===========================================================================

# Makefile name
MAKEFILE = l1aextract_modis.mk

# Progam to make
EXE     = l1aextract_modis

# Object modules for EXE
OBJ     = l1aextract_modis.o

# Include file locations
INCLUDE	= -I$(INCDIR)/utils -I$(LIB3_INC)

# LIBS needed to compile EXE
LIBS    = -L$(LIB3_LIB)

# Libraries required to link
LLIBS 	= -lmfhdf -ldf -ljpeg -lz

# Program specific additions to compiler flags
LOCAL_CFLAGS =

include $(MAKEFILE_APP_TEMPLATE)

