# ===========================================================================
#       Begin Program Specific Block
# ===========================================================================

MAKEFILE = Makefile

# Progam to make
EXE     = wgrib2

# Object modules for EXE
OBJ     = Alias.o \
CodeTable.o \
Code_Values.o \
Data.o \
Dump.o \
Ensemble.o \
File.o \
FlagTable.o \
Flush.o \
Grads.o \
Grib.o \
Help.o \
Inv.o \
Latlon.o \
Level.o \
Lola.o \
Macro.o \
Match.o \
Misc.o \
Netcdf.o \
Range.o \
Scan.o \
Sec0.o \
Sec1.o \
Sec3.o \
Sec4.o \
Sec5.o \
Sec6.o \
Sec_len.o \
Tigge.o \
Undefine.o \
VerfTime.o \
addtime.o \
cname.o \
earth.o \
fatal_error.o \
flt2ieee.o \
fnlist.o \
geo.o \
gribtab.o \
ieee2flt.o \
int8.o \
intpower.o \
missing.o \
rd_grib2_msg.o \
seekgrib2.o \
wgrib2.o \
wrt_grib.o \
wrtieee.o \
g2_free.o \
g2_getfld.o \
gbits.o \
g2_unpack1.o \
g2_unpack2.o \
g2_unpack3.o \
g2_unpack4.o \
g2_unpack5.o \
g2_unpack6.o \
g2_unpack7.o \
gridtemplates.o \
pdstemplates.o \
rdieee.o \
drstemplates.o \
simunpack.o \
comunpack.o \
specunpack.o \
int_power.o \
jpcunpack.o \
dec_jpeg2000.o \
jas_stream.o \
jas_image.o \
jas_seq.o \
jpc_dec.o \
jas_malloc.o \
jas_cm.o \
jas_string.o \
jpc_t1dec.o \
jas_debug.o \
jpc_tsfb.o \
jpc_mct.o \
jpc_cs.o \
jpc_t2cod.o \
jas_tvp.o \
jpc_t1cod.o \
jpc_bs.o \
jpc_mqdec.o \
jpc_tagtree.o \
jpc_t2dec.o \
jas_icc.o \
jpc_qmfb.o \
jpc_mqcod.o \
jas_iccdata.o \
jpc_math.o

# Include file locations
INCLUDE = \
-I$(LIB3_INC) \
-I.

# Library locations
LIBS    = \
-L$(LIB3_LIB)

# Libraries required to link
LLIBS   = \
-lmfhdf \
-lz

LOCAL_CFLAGS = -DUSE_TIGGE -DUSE_JPEG2000 -DUSE_REGEX

LD = $(CC)

include $(MAKEFILE_APP_TEMPLATE)
