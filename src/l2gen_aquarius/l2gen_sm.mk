# ===========================================================================
#       Begin Program Specific Block
# ===========================================================================

# Makefile name
MAKEFILE= l2gen_sm.mk

# Progam to make
EXE     = l2gen_sm

# Object modules for EXE

OBJ     = \
	l2gen_sm.o

# Include file locations
INCLUDE = -I$(LIB3_INC) -I$(INCDIR)/swfinc

# Library locations
LIBS 	= -L$(LIBDIR) -L$(LIB3_LIB)

# Libraries required to link
LLIBS 	= -laquarius -lhdfutils -lhdf5 -lgenutils -lz

include $(MAKEFILE_APP_TEMPLATE)
