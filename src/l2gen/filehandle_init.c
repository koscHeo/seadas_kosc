#include <string.h>
#include "filehandle.h"

void filehandle_init(filehandle *file)
{
    int32_t i;

    file->format        = -1;
    file->sensorID      = -1;
    file->subsensorID   = -1;
    strcpy(file->spatialResolution, "");
    file->modis_subset_compat = 0;

    file->length        =  0;
    file->spix          =  0;
    file->epix          = -1;
    file->dpix          =  1;
    file->npix          =  0;
    file->ndets         =  1;
    file->ctl_pt_incr   =  1;
    file->nscan         =  0;
    file->nbands        =  0;
    file->nbandsir      =  0;
    file->sd_id         =  0;
    file->tot_prod      =  0;
    file->percent_cloud =  0;
    file->percent_land  =  0;
    file->percent_water =  0;

    file->terrain_corrected = 0;
    file->sv_with_moon = 0;

    file->ocrvc_opt = 0;

    file->mode          = READ;
    file->bindx         = NULL;
    file->pro_control   = NULL;
    file->input_parms   = NULL;
    file->input         = NULL;
    file->calfile       = NULL;
    file->geofile       = NULL;
    file->mask_names    = NULL;

    strcpy(file->name,         "");
    strcpy(file->l2prod,       "");
    strcpy(file->def_l2prod,   "");

    for (i=0; i<MAXPROD; i++)
        strcpy(file->l2_prod_names[i],"");

    file->prodptr = NULL;

    strcpy(file->node_crossing_time,   "");
    file->orbit_number       = 0;
    file->orbit_node_lon     = -999.0;
    file->n_inprods = 0;

    memset(file->flag_cnt,0,NFLAGS*sizeof(int32_t));

    file->private_data = NULL;

    return;
}
