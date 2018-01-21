# ===========================================================================
#       Begin Program Specific Block
# ===========================================================================

# Makefile name
MAKEFILE= l2gen_aquarius.mk

# Progam to make
EXE     = l2gen_aquarius

# Object modules for EXE

MOD     = \
	l2gen_aquarius_module.mod

OBJ     = \
	l2gen_aquarius.o \
	get_prelaunch_calibration_coeffs.o \
	count_to_ta.o \
	l2gen_aquarius_module.o \
	get_rfi_auxiliary_data.o \
	rfi_detection.o \
	SolarFluxFlag.o

# Include file locations
INCLUDE = -I$(LIB3_INC) -I$(INCDIR)/swfinc -I$(INCDIR)

# Library locations
LIBS 	= -L$(LIBDIR) -L$(LIB3_LIB)

# Libraries required to link
LLIBS 	= -laqgeo -laqsss -laquarius -laqscatter -lscatterutils \
	-lgsl -lgslcblas -lhdfutils \
	-lhdf5_hl -lhdf5_fortran -lhdf5 -lgenutils -lmfhdf -ldf -ljpeg -lz

include $(MAKEFILE_APP_TEMPLATE)
