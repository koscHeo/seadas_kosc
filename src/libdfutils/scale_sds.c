#include <math.h>
#include <stdio.h>
#include <dfutils.h>
#include <genutils.h>

#define ROUND(x) ((((x) >= 0)?0.5:-0.5) + (x))

int16 *float2int16(float32 fbuf[], int32_t spix, int32_t npix, int incr, 
        float slope, float offset)
{
    static int32_t numPixels = 0;
    static  int16 *ibuf = NULL;
    float32 fval;
    int32_t    i;

    double maxval = slope*  32767  + offset;
    double minval = slope*(-32766) + offset;

    /*
     * -32766 is used for minval to allow -32767 to remain the sentinal bad value
     * Yes, this does leave one digit hanging out down low, but meh....
     */

    if(npix > numPixels) {
        numPixels = npix;
        if(ibuf)
            free(ibuf);
        if ( (ibuf = calloc(numPixels,sizeof(int16))) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: Unable to allocate buffer space.\n",
                    __FILE__,__LINE__);
            exit(1);
        }
    }

    for (i=0; i<npix; i++) {
        fval = fbuf[spix + i*incr];
        if (fval == BAD_FLT)
            ibuf[i] = BAD_INT;
        else if (fval >= maxval)
            ibuf[i] =  32767;
        else if (fval <= minval)
            ibuf[i] = -32766;
        else
            ibuf[i] = ROUND((fval - offset)/slope);
    }

    return(ibuf);
}

uint16 *float2uint16(float32 fbuf[], int32_t spix, int32_t npix, int incr,
        float slope, float offset)
{
    static int32_t numPixels = 0;
    static  uint16 *ibuf = NULL;
    float32 fval;
    int32_t    i;

    double maxval = slope*  65534  + offset;
    double minval = offset;

    if(npix > numPixels) {
        numPixels = npix;
        if(ibuf)
            free(ibuf);
        if ( (ibuf = calloc(numPixels,sizeof(uint16))) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: Unable to allocate buffer space.\n",
                    __FILE__,__LINE__);
            exit(1);
        }
    }

    for (i=0; i<npix; i++) {
        fval = fbuf[spix + i*incr];
        if (fval == BAD_FLT)
            ibuf[i] = BAD_UINT;
        else if (fval >= maxval)
            ibuf[i] =  65534;
        else if (fval <= minval)
            ibuf[i] = 0;
        else
            ibuf[i] = ROUND((fval - offset)/slope);
    }

    return(ibuf);
}

uint8_t *float2uint8(float32 fbuf[], int32_t spix, int32_t npix, int incr,
        float slope, float offset)
{
    static int32_t numPixels = 0;
    static  uint8_t *ibuf = NULL;
    float32 fval;
    int32_t    i;

    double maxval = slope*254 + offset;
    double minval = slope*  0 + offset;

    if(npix > numPixels) {
        numPixels = npix;
        if(ibuf)
            free(ibuf);
        if ( (ibuf = calloc(numPixels,sizeof(uint8_t))) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: Unable to allocate buffer space.\n",
                    __FILE__,__LINE__);
            exit(1);
        }
    }

    for (i=0; i<npix; i++) {
        fval = fbuf[spix + i*incr];
        if (fval == BAD_FLT)
            ibuf[i] = BAD_UBYTE;
        if (fval >= maxval)
            ibuf[i] = 254;
        else if (fval <= minval)
            ibuf[i] = 0;
        else
            ibuf[i] = ROUND((fval - offset)/slope);
    }

    return(ibuf);
}

int8_t *float2int8(float32 fbuf[], int32_t spix, int32_t npix, int incr,
        float slope, float offset)
{
    static int32_t numPixels = 0;
    static  int8_t *ibuf = NULL;
    float32 fval;
    int32_t    i;

    double maxval = slope*128 + offset;
    double minval = slope*-127 + offset;

    if(npix > numPixels) {
        numPixels = npix;
        if(ibuf)
            free(ibuf);
        if ( (ibuf = calloc(numPixels,sizeof(int8))) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: Unable to allocate buffer space.\n",
                    __FILE__,__LINE__);
            exit(1);
        }
    }

    for (i=0; i<npix; i++) {
        fval = fbuf[spix + i*incr];
        if (fval == BAD_FLT)
            ibuf[i] = BAD_BYTE;
        if (fval >= maxval)
            ibuf[i] = 128;
        else if (fval <= minval)
            ibuf[i] = -127;
        else
            ibuf[i] = ROUND((fval - offset)/slope);
    }

    return(ibuf);
}

VOIDP scale_sds(float *data, l2prodstr *p, int nbands)
{
    VOIDP pbuf;
    int   band = MIN(MAX(p->prod_ix,0),nbands-1);
    int   npix = p->dim[1];
    switch (p->datatype) {
    case DFNT_INT16:
        pbuf = (VOIDP) float2int16(data,band,npix,nbands,p->slope,p->offset);
        break;
    case DFNT_UINT16:
        pbuf = (VOIDP) float2uint16(data,band,npix,nbands,p->slope,p->offset);
        break;
    case DFNT_INT8:
        pbuf = (VOIDP) float2int8(data,band,npix,nbands,p->slope,p->offset);
        break;
    case DFNT_UINT8:
        pbuf = (VOIDP) float2uint8(data,band,npix,nbands,p->slope,p->offset);
        break;
    case DFNT_FLOAT32:
        fprintf(stderr,
               "-W- %s Line %d: Stubbornly refusing to scale floating point data: \n\t%s\n",__FILE__,__LINE__,p->title);
        pbuf = data;
        break;
    default:
        fprintf(stderr,
                "-W- %s Line %d: Unknown data type %d.\n",__FILE__,__LINE__,p->datatype);
        pbuf = data;
        break;
    }

    return (pbuf);
}

float *unscale_sds(VOIDP *data, l2prodstr *p, int32_t spix, int32_t npix, int incr)
{
    static int32_t numPixels = 0;
    static float *fbuf = NULL;

    float fval, *fptr;
    int16   ival, *iptr;
    uint8_t   bval, *bptr;
    int32_t    i;

    if(npix > numPixels) {
        numPixels = npix;
        if(fbuf)
            free(fbuf);
        if ( (fbuf = calloc(numPixels,sizeof(float))) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: Unable to allocate buffer space.\n",
                    __FILE__,__LINE__);
            exit(1);
        }
    }

    switch (p->datatype) {
    case DFNT_INT8:
    case DFNT_UINT8:
        bptr = (uint8_t *) data;
        for (i=0; i<npix; i++) {
            bval = bptr[spix + i*incr];
            fbuf[i] = bval*p->slope+p->offset;
        }
        break;
    case DFNT_INT16:
    case DFNT_UINT16:
        iptr = (int16 *) data;
        for (i=0; i<npix; i++) {
            ival = iptr[spix + i*incr];
            fbuf[i] = ival*p->slope+p->offset;
        }
        break;
    case DFNT_FLOAT32:
        fptr = (float32 *) data;
        for (i=0; i<npix; i++) {
            fval = fptr[spix + i*incr];
            fbuf[i] = fval*p->slope+p->offset;
        }
        break;
    default:
        fprintf(stderr,"-W- %s Line %d: Unknown data type %d product %d.\n",
                __FILE__,__LINE__,p->datatype,p->cat_ix);
        break;
    }

    return (fbuf);
}