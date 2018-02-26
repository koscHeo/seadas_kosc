#include "l12_proto.h"
/*  W. Robinson, SAIC, 12 Aug 2011 add VIIRS scn_fmt, margin_s  */

void cpl1rec(l1str *new, l1str *old)
{
    new->sensorID      = old->sensorID;
    new->npix          = old->npix;
    new->nscans        = old->nscans;
    new->spix          = old->spix;
    new->epix          = old->epix;
    new->dpix          = old->dpix;
    new->sscan         = old->sscan;
    new->escan         = old->escan;
    new->dscan         = old->dscan;
    new->nbands        = old->nbands;
    new->nbandsir      = old->nbandsir;
    new->ndets         = old->ndets;
    new->bindx         = old->bindx;
    new->length        = old->length;
    new->iscan         = old->iscan;
    new->detnum        = old->detnum;
    new->mside         = old->mside;
    new->n_inprods     = old->n_inprods;

    new->scn_fmt       = old->scn_fmt;
    new->margin_s      = old->margin_s;

    new->input         = old->input;

    new->landMaskOn    = old->landMaskOn;
    new->bathMaskOn    = old->bathMaskOn;
    new->cloudMaskOn   = old->cloudMaskOn;
    new->glintMaskOn   = old->glintMaskOn;
    new->hiltMaskOn    = old->hiltMaskOn;
    new->stlightMaskOn = old->stlightMaskOn;
    new->senzMaskOn    = old->senzMaskOn;
    new->solzMaskOn    = old->solzMaskOn;

    new->fsol          = old->fsol;
    new->tilt          = old->tilt;

    memcpy(new->iwave, old->iwave, old->nbands*sizeof(int32_t ));
    memcpy(new->fwave, old->fwave, old->nbands*sizeof(float));
    memcpy(new->Fo,    old->Fo,    old->nbands*sizeof(float));
    memcpy(new->Fobar, old->Fobar, old->nbands*sizeof(float));
    memcpy(new->Fonom, old->Fonom, old->nbands*sizeof(float));
    memcpy(new->Tau_r, old->Tau_r, old->nbands*sizeof(float));
    memcpy(new->k_oz,  old->k_oz,  old->nbands*sizeof(float));
    memcpy(new->aw,    old->aw,    old->nbands*sizeof(float));
    memcpy(new->bbw,   old->bbw,   old->nbands*sizeof(float));

    memcpy(new->data,old->data,new->length);
}    

