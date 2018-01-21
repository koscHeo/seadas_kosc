# ===========================================================================
#       Begin Program Specific Block
# ===========================================================================

# Makefile name
MAKEFILE= l2gen_scatter.mk

# Progam to make
EXE     = l2gen_scatter

# Object modules for EXE

OBJ     = \
	l2gen_scatter.o

# Include file locations
INCLUDE = -I$(LIB3_INC) -I$(LIB3_LIB) -I$(INCDIR)/swfinc -I$(INCDIR)/aquarius -I$(INCDIR)

# Library locations
LIBS 	= -L$(LIBDIR) -L$(LIB3_LIB)

# Libraries required to link
LLIBS 	= -laqscatter -lscatterutils -laqsss \
	-lhdfutils \
	-lhdf5_hl -lhdf5_fortran -lhdf5 -lgenutils -lmfhdf -ldf -ljpeg -lz

LOCAL_FFLAGS = -O0

include $(MAKEFILE_APP_TEMPLATE)
