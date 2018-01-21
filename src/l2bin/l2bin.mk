# ===========================================================================
#	Begin Program Specific Block
# ===========================================================================

# Makefile name
MAKEFILE= Makefile

# Progam to make
EXE	= l2bin

# Object modules for EXE

OBJ    	= l2bin.o \
	  l2bin_input.o


# Library locations
LIBS 	= -L$(LIBDIR) -L$(LIB3_LIB)

# LIBS needed to compile EXE
LLIBS 	= \
	  -ll2 \
	  -lbin \
          -lgenutils \
          -lseawifs \
          -lhdfutils \
          -lnav \
          -lmfhdf \
          -ldf \
          -ljpeg \
          -lz

# Include file locations
INCLUDE = -I$(INCDIR)/swfinc -I$(INCDIR)/utils -I$(LIB3_INC)

#LD = gfortran

include $(MAKEFILE_APP_TEMPLATE)
