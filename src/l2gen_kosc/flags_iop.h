#ifndef _FLAGS_IOP_H
#define _FLAGS_IOP_H

#include "l12_parms.h"
#include "hdf.h"

/* flag bit settings */
#define IOPF_ISMASKED 0x0001   // pixel masked
#define IOPF_FAILED   0x0002   // algorithm signals failure
#define IOPF_MAXITER  0x0004   // max iterations reached
#define IOPF_BADRRS   0x0008   // insufficient valid Rrs
#define IOPF_NAN      0x0010   // inversion returned NAN
#define IOPF_RRSDIFF  0x0020   // mean Rrs diff exceeded threshold
#define IOPF_ALO      0x0040   // retrieved a below threshold
#define IOPF_AHI      0x0080   // retrieved a above threshold
#define IOPF_APHLO    0x0100   // retrieved aph below threshold
#define IOPF_APHHI    0x0200   // retrieved aph above threshold
#define IOPF_ADGLO    0x0400   // retrieved adg below threshold
#define IOPF_ADGHI    0x0800   // retrieved adg above threshold
#define IOPF_BBLO     0x1000   // retrieved bb below threshold
#define IOPF_BBHI     0x2000   // retrieved abb above threshold
#define IOPF_BBPLO    0x4000   // retrieved bbp below threshold
#define IOPF_BBPHI    0x8000   // retrieved bbp above threshold

static const char *giop_flag_lname[NGIOPFLAGS] =     {"ISMASKED",
                                                        "FAILED",
                                                        "MAXITER",
                                                        "BADRRS",
                                                        "NAN",
                                                        "RRSDIFF",
                                                        "ALO",
                                                        "AHI",
                                                        "APHLO",
                                                        "APHHI",
                                                        "ADGLO",
                                                        "ADGHI",
                                                        "BBLO",
                                                        "BBHI",
                                                        "BBPLO",
                                                        "BBPHI"
};
typedef struct iopflagctl_struc {

    int32 nwave;
    int32 maxiter;

    float32 *a_lo  ;
    float32 *a_hi  ;
    float32 *a_on  ;

    float32 *aph_lo;
    float32 *aph_hi;
    float32 *aph_on;

    float32 *adg_lo;
    float32 *adg_hi;
    float32 *adg_on;

    float32 *bb_lo ;
    float32 *bb_hi ;
    float32 *bb_on ;

    float32 *bbp_lo;
    float32 *bbp_hi;
    float32 *bbp_on;

} iopfstr;

void set_iop_flag(float32 wave[], int32 nwave, 
                  float32 a[], float32 aph[], float32 adg[], 
                  float32 bb[], float32 bbp[],int16 *flag);

#endif
