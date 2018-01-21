# ===========================================================================
#	Begin Program Specific Block
# ===========================================================================

# Makefile name
MAKEFILE= Makefile

# Progam to make
EXE	= l3bin

# Object modules for EXE

OBJ    	= l3bin.o \
	  l3bin_input.o

# Library locations
LIBS 	= -L$(LIBDIR) -L$(LIB3_LIB)

# LIBS needed to compile EXE
LLIBS 	= \
	  -lbin++ \
	  -lbin \
          -lgenutils \
          -lhdfutils \
          -lmfhdf \
          -ldf \
          -lhdf5 \
          -ljpeg \
          -lz

# Include file locations
INCLUDE = -I$(INCDIR)/swfinc -I$(INCDIR)/utils -I$(LIB3_INC)

LD = $(CC)

include $(MAKEFILE_APP_TEMPLATE)

ifdef ppcDarwinArchitecture
  LOCAL_CFLAGS = -O0
  LOCAL_LDFLAGS = -Wl,-stack_size -Wl,0x40000000 -Wl,-stack_addr -Wl,0x70000000
endif
