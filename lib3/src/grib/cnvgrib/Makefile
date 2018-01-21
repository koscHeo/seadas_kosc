# ===========================================================================
#       Begin Program Specific Block
# ===========================================================================

# Makefile name
MAKEFILE= Makefile

# Progam to make
EXE     = cnvgrib

# Object modules for EXE

OBJ     = \
	gds2gdt.o \
	pds2pdt.o \
	pds2pdtens.o \
	cnvgrib.o \
	cnv12.o \
	cnv22.o \
	cnv21.o \
	gdt2gds.o \
	makepds.o \
	putgbexn.o \
	makepdsens.o \
	setbit.o


#INCLUDE = -I$(INCDIR)/swfinc -I../libg2
INCLUDE = -I$(INCDIR)/swfinc -I../libg2

# LIBS needed to compile EXE
LIBS    = -L$(LIBDIR)

# Libraries required to link
LLIBS 	= -lg2 -lw3

#FC = gfortran
#LD = gfortran

include $(MAKEFILE_APP_TEMPLATE)

