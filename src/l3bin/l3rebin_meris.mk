# ===========================================================================
#       Begin Program Specific Block
# ===========================================================================

# Makefile name
MAKEFILE= l3rebin_meris.mk

# Progam to make
EXE     = l3rebin_meris

# Object modules for EXE

OBJ     = l3rebin_meris.o

# Include file locations
INCLUDE = -I$(INCDIR)/swfinc -I$(LIB3_INC)

# Library locations
LIBS    = -L$(LIBDIR) -L$(LIB3_LIB)


# Libraries required to link
LLIBS   = \
          -lbin++ \
          -lbin \
          -lgenutils \
          -lhdfutils \
          -lmfhdf \
          -ldf \
          -lnetcdf \
          -lhdf5_hl \
          -lhdf5 \
          -ljpeg \
          -lz


include $(MAKEFILE_APP_TEMPLATE)


