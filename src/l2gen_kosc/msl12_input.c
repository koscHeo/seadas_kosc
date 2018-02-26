#include "l12_proto.h"
#include <string.h>
#include <clo.h>
#include <sensorInfo.h>

#include <assert.h>
#include <strings.h>
#include "version.h"

static int32_t numTauas = -1;
static int32_t numGains = -1;
static int32_t numGainUncs = -1;
static int32_t numOffsets = -1;

static int enableFileDecending = 1;

int defnaermodels = 80;
char defaermodels[][32] = {"r30f95v01", "r30f80v01", "r30f50v01", "r30f30v01", "r30f20v01", "r30f10v01", "r30f05v01", "r30f02v01", "r30f01v01", "r30f00v01",
        "r50f95v01", "r50f80v01", "r50f50v01", "r50f30v01", "r50f20v01", "r50f10v01", "r50f05v01", "r50f02v01", "r50f01v01", "r50f00v01",
        "r70f95v01", "r70f80v01", "r70f50v01", "r70f30v01", "r70f20v01", "r70f10v01", "r70f05v01", "r70f02v01", "r70f01v01", "r70f00v01",
        "r75f95v01", "r75f80v01", "r75f50v01", "r75f30v01", "r75f20v01", "r75f10v01", "r75f05v01", "r75f02v01", "r75f01v01", "r75f00v01",
        "r80f95v01", "r80f80v01", "r80f50v01", "r80f30v01", "r80f20v01", "r80f10v01", "r80f05v01", "r80f02v01", "r80f01v01", "r80f00v01",
        "r85f95v01", "r85f80v01", "r85f50v01", "r85f30v01", "r85f20v01", "r85f10v01", "r85f05v01", "r85f02v01", "r85f01v01", "r85f00v01",
        "r90f95v01", "r90f80v01", "r90f50v01", "r90f30v01", "r90f20v01", "r90f10v01", "r90f05v01", "r90f02v01", "r90f01v01", "r90f00v01",
        "r95f95v01", "r95f80v01", "r95f50v01", "r95f30v01", "r95f20v01", "r95f10v01", "r95f05v01", "r95f02v01", "r95f01v01", "r95f00v01"};

// need a place to store the default product lists
static char default_l2prod[MAX_OFILES][PRODSTRLEN];

static char default_flaguse[1024];

// store the name of program we are running.
static char mainProgramName[50];

static char *l2gen_optionKeys[] = {
    "-help",
    "-version",
    "-dump_options",
    "-dump_options_paramfile",
    "-dump_options_xmlfile",
    "par",
    "pversion",
    "suite",
    "eval",
    "ifile",
    "ilist",
    "geofile",
    "ofile",
    "oformat",
    "il2file",
    "tgtfile",
    "aerfile",
    "metafile",
    "l2prod",
    "spixl",
    "epixl",
    "dpixl",
    "sline",
    "eline",
    "dline",
    "ctl_pt_incr",
    "proc_ocean",
    "proc_land",
    "proc_sst",
    "atmocor",
    "mode",
    "aer_opt",
    "aer_wave_short",
    "aer_wave_long",
    "aer_swir_short",
    "aer_swir_long",
    "aer_rrs_short",
    "aer_rrs_long",
    "aermodmin",
    "aermodmax",
    "aermodrat",
    "aer_angstrom",
    "aer_iter_max",
    "mumm_alpha",
    "mumm_gamma",
    "mumm_epsilon",
    "absaer_opt",
    "glint_opt",
    "outband_opt",
    "oxaband_opt",
    "cirrus_opt",
    "filter_opt",
    "filter_file",
    "brdf_opt",
    "fqfile",
    "parfile",
    "gas_opt",
    "atrem_opt",
    "atrem_full",
    "atrem_geom",
    "atrem_model",
    "iop_opt",
    "seawater_opt",
    "polfile",
    "pol_opt",
    "ocrvc_opt",
    "rad_opt",
    "xcalfile",
    "xcal_opt",
    "xcal_wave",
    "band_shift_opt",
    "add_noise_sigma",
    "noise_scale",
    "wiggle_band",
    "wiggle_by",
    "resolution",
    "giop_aph_opt",
    "giop_aph_file",
    "giop_aph_s",
    "giop_adg_opt",
    "giop_adg_file",
    "giop_adg_s",
    "giop_bbp_opt",
    "giop_bbp_file",
    "giop_bbp_s",
    "giop_acdom_opt",
    "giop_acdom_file",
    "giop_anap_opt",
    "giop_anap_file",
    "giop_bbph_opt",
    "giop_bbph_file",
    "giop_bbnap_opt",
    "giop_bbnap_file",
    "giop_rrs_opt",
    "giop_rrs_diff",
    "giop_grd",
    "giop_wave",
    "giop_maxiter",
    "giop_fit_opt",
    "gsm_opt",
    "gsm_fit",
    "gsm_adg_s",
    "gsm_bbp_s",
    "gsm_aphw",
    "gsm_aphs",
    "qaa_adg_s",
    "qaa_wave",
    "chloc2_wave",
    "chloc2_coef",
    "chloc3_wave",
    "chloc3_coef",
    "chloc4_wave",
    "chloc4_coef",
    "kd2_wave",
    "kd2_coef",
    "flh_offset",
    "vcnnfile",
    "picfile",
    "owtfile",
    "owtchlerrfile",
    "aermodfile",
    "aermodels",
    "met1",
    "met2",
    "met3",
    "ozone1",
    "ozone2",
    "ozone3",
    "anc_cor_file",
    "land",
    "water",
    "demfile",
    "elevfile",
    "elev_auxfile",
    "icefile",
    "ice_threshold",
    "sstcoeffile",
    "sstssesfile",
    "sst4coeffile",
    "sst4ssesfile",
    "sst3coeffile",
    "sst3ssesfile",
    "sstfile",
    "sstreftype",
    "sstrefdif",
    "viirsnv7",
    "viirsnosisaf",
    "newavhrrcal",
    "ch22detcor",
    "ch23detcor",
    "no2file",
    "alphafile",
    "tauafile",
    "cldfile",
    "calfile",
    "offset",
    "gain",
    "flaguse",
    "xcalbox",
    "xcalboxcenter",
    "xcalpervalid",
    "xcalsubsmpl",
    "chlthreshold",
    "aotthreshold",
    "coccolith",
    "cirrus_thresh",
    "taua",
    "cloud_thresh",
    "cloud_wave",
    "cloud_eps",
    "glint_thresh",
    "absaer",
    "rhoamin",
    "sunzen",
    "satzen",
    "epsmin",
    "epsmax",
    "tauamax",
    "nLwmin",
    "hipol",
    "wsmax",
    "windspeed",
    "windangle",
    "pressure",
    "ozone",
    "relhumid",
    "watervapor",
    "maskland",
    "maskbath",
    "maskcloud",
    "maskglint",
    "masksunzen",
    "masksatzen",
    "maskhilt",
    "maskstlight",
    "sl_frac",
    "sl_pixl",
    "vcal_opt",
    "vcal_chl",
    "vcal_solz",
    "vcal_nLw",
    "vcal_Lw",
    "vcal_depth",
    "vcal_min_nbin",
    "vcal_min_nscene",
    "owmcfile",
    "north",
    "south",
    "east",
    "west",
    "xbox",
    "ybox",
    "subsamp",
    "prodxmlfile",
    "breflectfile",
    "deflate",
    "raman_opt",
    "viirscalparfile",
    NULL
};

static char *l1bgen_optionKeys[] = {
        "-help",
        "-version",
        "-dump_options",
        "-dump_options_paramfile",
        "-dump_options_xmlfile",
        "par",
        "pversion",
        "suite",
        "ifile",
        "ofile",
        "oformat",
        "fqfile",
        "parfile",
        "spixl",
        "epixl",
        "sline",
        "eline",
        "calfile",
        "xcalfile",
        "gain",
        "offset",
        "sl_pixl",
        "sl_frac",
        NULL
};

static char *l1mapgen_optionKeys[] = {
        "-help",
        "-version",
        "-dump_options",
        "-dump_options_paramfile",
        "-dump_options_xmlfile",
        "par",
        "pversion",
        "suite",
        "ifile",
        "geofile",
        "resolution",
        "ofile",
        "oformat",
        "north",
        "south",
        "east",
        "west",
        "width",
        "threshold",
        "rgb",
        "atmocor",
        "datamin",
        "datamax",
        "stype",
        "help_true_color",
        "cirrus_opt",
        "atrem_opt",
        "atrem_full",
        "atrem_geom",
        "atrem_model",
        NULL
};

static char *l1brsgen_optionKeys[] = {
        "-help",
        "-version",
        "-dump_options",
        "-dump_options_paramfile",
        "-dump_options_xmlfile",
        "par",
        "pversion",
        "suite",
        "ifile",
        "geofile",
        "resolution",
        "ofile",
        "oformat",
        "oformat_depth",
        "spixl",
        "epixl",
        "sline",
        "eline",
        "subsamp",
        "rgb",
        "atmocor",
        "datamin",
        "datamax",
        "stype",
        "help_true_color",
        "cirrus_opt",
        "atrem_opt",
        "atrem_full",
        "atrem_geom",
        "atrem_model",
        NULL
};

static char *l1det2det_optionKeys[] = {
        "-help",
        "-version",
        "-dump_options",
        "-dump_options_paramfile",
        "-dump_options_xmlfile",
        "par",
        "pversion",
        "suite",
        "ifile",
        "help_ifile",
        "ofile",
        "oformat",
        "geofile",
        "l2prod",
        "spixl",
        "epixl",
        "sline",
        "eline",
        "ybox",
        "chlthreshold",
        "aotthreshold",
        "glint_thresh",
        "cloud_thresh",
        "flaguse",
        NULL
};

static char *vcalmerge_optionKeys[] = {
        "-help",
        "-version",
        "-dump_options",
        "-dump_options_paramfile",
        "-dump_options_xmlfile",
        "par",
        "pversion",
        "ifile",
        "help_ifile",
        "ofile",
        "oformat",
        "spixl",
        "epixl",
        "sline",
        "eline",
        "flaguse",
        "deflate",
        NULL
};

 float* calloc_nbandsf( int32_t nbands, float *nbarray, float init_val) {
    int i;

    if ((nbarray = (float *) calloc(nbands, sizeof(float))) == NULL) {
        printf("-E- : Error allocating float memory in alloc_nbandsf\n");
        exit(FATAL_ERROR);
    }
    for (i=0; i< nbands; i++)
        nbarray[i] = init_val;

    return nbarray;
}
 int32_t* calloc_nbandsi32t( int32_t nbands, int32_t *nbarray, int32_t init_val) {
    int i;
    if ((nbarray = (int32_t *) calloc(nbands, sizeof(int32_t))) == NULL) {
        printf("-E- : Error allocating float memory in alloc_nbandsi\n");
        exit(FATAL_ERROR);
    }
    for (i=0; i< nbands; i++)
        nbarray[i] = init_val;

    return nbarray;
}
 int* calloc_nbandsi( int32_t nbands, int *nbarray, int init_val) {
    int i;
    if ((nbarray = (int *) calloc(nbands, sizeof(int))) == NULL) {
        printf("-E- : Error allocating float memory in alloc_nbandsi\n");
        exit(FATAL_ERROR);
    }
    for (i=0; i< nbands; i++)
        nbarray[i] = init_val;

    return nbarray;
}
void msl12_input_nbands_init(instr *input) {
    /* allocate and initialize dynamic arrays in input struc */

    input->gsm_aphs     = calloc_nbandsf(input->nbands, input->gsm_aphs,     -1.0);
    input->gsm_aphw     = calloc_nbandsf(input->nbands, input->gsm_aphw,     -1.0);
    input->giop_wave    = calloc_nbandsf(input->nbands, input->giop_wave,    -1.0);
    input->giop_rrs_unc = calloc_nbandsf(input->nbands, input->giop_rrs_unc, -1.0);
    input->xcal_opt     = calloc_nbandsi32t(input->nbands, input->xcal_opt,     0);
    input->xcal_wave    = calloc_nbandsf(input->nbands, input->xcal_wave,    -1.0);
    input->gain         = calloc_nbandsf(input->nbands, input->gain,          1.0);
    input->gain_unc     = calloc_nbandsf(input->nbands, input->gain_unc,      0.0);
    input->offset       = calloc_nbandsf(input->nbands, input->offset,        0.0);
    input->taua         = calloc_nbandsf(input->nbands, input->taua,          0.0);
    input->vcal_nLw     = calloc_nbandsf(input->nbands, input->vcal_nLw,      0.0);
    input->vcal_Lw      = calloc_nbandsf(input->nbands, input->vcal_Lw,       0.0);

}

/*-----------------------------------------------------------------------------
    Function:  msl12_input_init

    Returns:   int (status)
        The return code is a negative value if any error occurs, otherwise,
        returns 0.

    Description:
        Set default values for input structure.

    Parameters: (in calling order)
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        instr     input     O   structure variable for inputs

----------------------------------------------------------------------------*/

void msl12_input_init(instr *input) {
    int i;

    /* Start fresh                                                      */
    memset(input, 0, sizeof(*input));

    /*                                                                  */
    /* Set input values to defaults                                     */
    /*                                                                  */
    input->mode = FORWARD;
    input->sensorID = -1;
    input->subsensorID = -1;
    input->format = -1;
    input->modis_subset_compat = 0;

    input->spixl = 1;
    input->epixl = -1;
    input->dpixl = 1;
    input->sline = 1;
    input->eline = -1;
    input->dline = 1;
    input->ctl_pt_incr = 8;

    input->sl_pixl = -1;
    input->sl_frac = 0.25;

    input->proc_ocean = 1;
    input->proc_sst = -1;
    input->proc_land = 0;
    input->atmocor = 1;
    input->resolution = -1;

    input->glint_opt = 1;
    input->aer_iter_max = 10;
    input->brdf_opt = -1;
    input->gas_opt = 1;
    input->atrem_opt = 0;
    input->atrem_full = 0;
    input->atrem_geom = 0;
    input->atrem_model = 0;
    input->ocrvc_opt = 0;
    input->iop_opt = IOPNONE;
    input->pol_opt = -1;
    input->rad_opt = -1;
    input->absaer_opt = -1;
    input->cirrus_opt = 0;
    input->gsm_opt = 0;
    input->gsm_fit = 0;
    input->gsm_adg_s = 0.02061;
    input->gsm_bbp_s = 1.03373;

    input->giop_maxiter = -1;
    input->giop_fit_opt = -1;
    input->giop_aph_opt = -1;
    input->giop_adg_opt = -1;
    input->giop_bbp_opt = -1;
    input->giop_acdom_opt= -1;
    input->giop_anap_opt= -1;
    input->giop_bbph_opt= -1;
    input->giop_bbnap_opt= -1;
    input->giop_rrs_opt = -1;
    input->giop_rrs_diff = -1.0;
    input->giop_iterate = 0;
    input->giop_aph_s = -1000.0;
    input->giop_adg_s = -1000.0;
    input->giop_bbp_s = -1000.0;
    input->giop_aph_w = -1.0;
    input->giop_adg_w = -1.0;
    input->giop_bbp_w = -1.0;
    input->giop_grd[0] = -1000.0;
    input->giop_grd[1] = -1000.0;

    input->qaa_adg_s = 0.015;

    input->flh_offset = 0.0;

    input->seawater_opt = 0;
    input->aer_opt = 99;
    input->outband_opt = 99;
    input->oxaband_opt = 99;
    input->filter_opt = 99;

    input->aer_wave_short = 765;
    input->aer_wave_long = 865;
    input->aer_swir_short = -1;
    input->aer_swir_long = -1;
    input->aer_rrs_short = -1.0;
    input->aer_rrs_long = -1.0;
    input->aer_angstrom = -999.0;

    input->xcal_nwave = 0;
    input->band_shift_opt = 0;
    input->add_noise_sigma=0.0;
    input->noise_scale = 1.0;
    input->wiggle_band=-1;
    input->wiggle_by=0.0;

    input->vcal_chl = -1.0;
    input->vcal_solz = -1.0;
    input->vcal_opt = -1;

    input->aermodrat = 0.0;
    input->aermodmin = -1;
    input->aermodmax = -1;

    input->albedo = -1.0;
    input->cloud_wave = 865.0;
    input->cloud_eps = -1.0;
    input->absaer = 0.0;
    input->rhoamin = 0.0001;
    input->glint = 0.005;
    input->sunzen = 75.0;
    input->satzen = 60.0;
    input->epsmin = 0.85;
    input->epsmax = 1.35;
    input->tauamax = 0.30;
    input->nlwmin = 0.15;
    input->hipol = 0.50;
    input->wsmax = 8.0;

    input->mumm_alpha = 1.72;
    input->mumm_gamma = 1.00;
    input->mumm_epsilon = 1.00;

    input->windspeed = -1000;
    input->windangle = -1000;
    input->pressure = -1000;
    input->ozone = -1000;
    input->relhumid = -1000;
    input->watervapor = -1000;
    input->ice_threshold = 0.1;

    input->landmask = 1;
    input->bathmask = 0;
    input->cloudmask = 1;
    input->glintmask = 0;
    input->sunzenmask = 0;
    input->satzenmask = 0;
    input->hiltmask = 1;
    input->stlightmask = 1;

    for (i = 0; i < 8; i++) input->coccolith[0] = 0.0;
    for (i = 0; i < 2; i++) input->cirrus_thresh[i] = -1;
    for (i = 0; i < 2; i++) input->chloc2w [i] = -1;
    for (i = 0; i < 3; i++) input->chloc3w [i] = -1;
    for (i = 0; i < 4; i++) input->chloc4w [i] = -1;
    for (i = 0; i < 2; i++) input->kd2w [i] = -1;

    strcpy(input->pversion, "Unspecified");
    strcpy(input->suite, "");

    fctl_init(&(input->fctl));

    /* for inverse (calibration) modes */
    input->chlthreshold = CHL_MAX;
    input->aotthreshold = AOT_MAX;
    input->maxpointdist = 0.0;
    strcpy(default_flaguse, "ATMFAIL,LAND,HIGLINT,HILT,HISATZEN,STRAYLIGHT,CLDICE,COCCOLITH,LOWLW,CHLFAIL,NAVWARN,ABSAER,MAXAERITER,ATMWARN,HISOLZEN,NAVFAIL");
    strcpy(input->flaguse, default_flaguse);

    input->xcalbox = 0;
    input->xcalboxcenter[0] = 0;
    input->xcalboxcenter[1] = 0;
    input->xcalpervalid = 0;
    input->xcalsubsmpl = 1;

    input->vcal_depth = -1000;
    input->vcal_min_nbin = 4;
    input->vcal_min_nscene = 3;

    /* for sst */
    input->sstreftype  = 0;
    input->sstrefdif = 100.0;   /* some large number so it doesn't hurt */
    input->viirsnv7 = -1;        /* VIIRSN v7 high satz latband equation */
    input->viirsnosisaf = 0;    /* VIIRSN OSI-SAF equation */
    input->newavhrrcal = 0;
    for (i = 0; i < 10; i++) input->ch22detcor[i] = 1.0;    /* for modis: Ltir = Ltir / detcor */
    for (i = 0; i < 10; i++) input->ch23detcor[i] = 1.0;

    /* for l1mapgen and l1brsgen */
    input->datamin = 0.01;
    input->datamax = 0.9;
    input->west = -999;
    input->east = -999;
    input->north = -999;
    input->south = -999;
    input->width = 600;
    input->threshold = 0.1;
    input->subsamp = 1;
    input->stype = 0;
    for (i = 0; i < 3; i++) {
        input->rgb[i] = 1;
    }
    input->xbox = -1;
    input->ybox = -1;

    input->deflate = 0;
    
    input->raman_opt = 0;
    
    return;
}

//-----------------------------------------------------------------------

/** CLO callback function for the "par" option.  Loads the parameter file
 into the list that is stored in option->cb_data */
void par_option_cb(clo_option_t *option) {
    if (enableFileDecending)
        clo_readFile((clo_optionList_t*) option->cb_data, option->valStr);
}

//-----------------------------------------------------------------------

/** add all of the accepted command line options to list */
int l2gen_init_options(clo_optionList_t* list, const char* prog) {
    char tmpStr[2048];
    char tmpStr1[2048];
    char tmpStr2[32];
    clo_option_t* option;
    int i;

    // set the min program name
    strcpy(mainProgramName, prog);

    if (!strcmp(prog, "msl12")) {
        clo_setSelectOptionKeys(l2gen_optionKeys);
    } else if (!strcmp(prog, "l2gen")) {
        clo_setSelectOptionKeys(l2gen_optionKeys);
    } else if (!strcmp(prog, "l3gen")) {
        clo_setSelectOptionKeys(l2gen_optionKeys);
    } else if (!strcmp(prog, "l1bgen_generic")) {
        clo_setSelectOptionKeys(l1bgen_optionKeys);
    } else if (!strcmp(prog, "l1mapgen")) {
        clo_setSelectOptionKeys(l1mapgen_optionKeys);
    } else if (!strcmp(prog, "l1brsgen")) {
        clo_setSelectOptionKeys(l1brsgen_optionKeys);
    } else if (!strcmp(prog, "l1det2det")) {
        clo_setSelectOptionKeys(l1det2det_optionKeys);
    } else if (!strcmp(prog, "vcalmerge")) {
        clo_setSelectOptionKeys(vcalmerge_optionKeys);
    }

    sprintf(tmpStr, "%s %d.%d.%d-r%d (%s %s)", prog, L2GEN_VERSION_MAJOR,L2GEN_VERSION_MINOR,L2GEN_VERSION_PATCH_LEVEL,SVN_REVISION, __DATE__, __TIME__);
    clo_setVersion(tmpStr);
    clo_addXmlProgramMetadata("progressRegex", "Processing scan .+?\\((\\d+) of (\\d+)\\)");

    sprintf(tmpStr, "Usage: %s argument-list\n\n", prog);
    strcat(tmpStr, "  The argument-list is a set of keyword=value pairs. The arguments can\n");
    strcat(tmpStr, "  be specified on the commandline, or put into a parameter file, or the\n");
    strcat(tmpStr, "  two methods can be used together, with commandline over-riding.\n\n");
    strcat(tmpStr, "  return value: 0=OK, 1=error, 110=north,south,east,west does not intersect\n");
    strcat(tmpStr, "  file data.\n\n");
    strcat(tmpStr, "The list of valid keywords follows:\n");
    clo_setHelpStr(tmpStr);

    // add the par option and add the callback function
    option = clo_addOption(list, "par", CLO_TYPE_IFILE, NULL, "input parameter file");
    option->cb_data = (void*) list;
    option->cb = par_option_cb;

    strcpy(tmpStr, "processing version string");
    clo_addOption(list, "pversion", CLO_TYPE_STRING, "Unspecified", tmpStr);

    strcpy(tmpStr, "product suite string for loading\n");
    strcat(tmpStr, "        suite-specific defaults");
    clo_addOption(list, "suite", CLO_TYPE_STRING, "OC", tmpStr);

    strcpy(tmpStr, "evaluation bitmask\n");
    strcat(tmpStr, "          0: standard processing\n");
    strcat(tmpStr, "          1: init to old aerosol models\n");
    strcat(tmpStr, "          2: enables MODIS and MERIS cloud Mask for HABS\n");
    strcat(tmpStr, "         16: enables MODIS cirrus mask\n");
    strcat(tmpStr, "         32: use test sensor info file\n");
    strcat(tmpStr, "         64: use test rayleigh tables\n");
    strcat(tmpStr, "        128: use test aerosol tables\n");
    strcat(tmpStr, "        256: use test polarization tables\n");
    strcat(tmpStr, "       1024: mask modis mirror-side 1 (navfail)\n");
    strcat(tmpStr, "       2048: mask modis mirror-side 2 (navfail)\n");
    strcat(tmpStr, "       4096: don't apply 'cold-only' or equatorial aerosol tests for SST\n");
    strcat(tmpStr, "       8192: use alt sensor info file in eval\n");
    strcat(tmpStr, "      32768: enables spherical path geom for dtran");
    clo_addOption(list, "eval", CLO_TYPE_INT, "0", tmpStr);

    strcpy(tmpStr, "input L1 file name");
    option = clo_addOption(list, "ifile", CLO_TYPE_IFILE, NULL, tmpStr);
    clo_addOptionAlias(option, "ifile1");

    strcpy(tmpStr, "ifile[#]=input L1 file names (1-original and 2-vicarious) to be cross-calibrated\n");
    strcat(tmpStr, "      or input HDF file names containing cross-calibration pixels");
    clo_addOption(list, "help_ifile", CLO_TYPE_HELP, NULL, tmpStr);

    clo_addOption(list, "ilist", CLO_TYPE_IFILE, NULL, "file containing list of input files, one per line");

    strcpy(tmpStr, "input L1 geolocation file name (MODIS/VIIRS only)");
    clo_addOption(list, "geofile", CLO_TYPE_IFILE, NULL, tmpStr);

    strcpy(tmpStr, "VIIRS L1A calibration parameter file name (VIIRS only)");
    clo_addOption(list, "viirscalparfile", CLO_TYPE_IFILE, NULL, tmpStr);

    if (!strcmp(prog, "l1mapgen")) {
        option = clo_addOption(list, "ofile", CLO_TYPE_OFILE, "output", "output file name");
        clo_addOptionAlias(option, "ofile1");

        strcpy(tmpStr, "output file format\n");
        strcat(tmpStr, "        ppm:  output a netPBM PPM file\n");
        strcat(tmpStr, "        png:  output a PNG file\n");
        strcat(tmpStr, "        tiff: output a geoTIFF file");
        clo_addOption(list, "oformat", CLO_TYPE_STRING, "ppm", tmpStr);

    } else if (!strcmp(prog, "l1brsgen")) {
        option = clo_addOption(list, "ofile", CLO_TYPE_OFILE, "output", "output file name");
        clo_addOptionAlias(option, "ofile1");

        strcpy(tmpStr, "output file format\n");
        strcat(tmpStr, "        hdf4: output a HDF4 file\n");
        strcat(tmpStr, "        bin:  output a flat binary file\n");
        strcat(tmpStr, "        png:  output a PNG file\n");
        strcat(tmpStr, "        ppm:  output a netPBM PPM file");
        clo_addOption(list, "oformat", CLO_TYPE_STRING, "hdf4", tmpStr);

    } else if (!strcmp(prog, "l1bgen_generic") || !strcmp(prog,"l1det2det")) {
        option = clo_addOption(list, "ofile", CLO_TYPE_OFILE, "output", "output file name");
        clo_addOptionAlias(option, "ofile1");

        strcpy(tmpStr, "output file format\n");
        strcat(tmpStr, "        netcdf4: output a netCDF version 4 file\n");
        strcat(tmpStr, "        hdf4:    output a HDF version 4 file");
        clo_addOption(list, "oformat", CLO_TYPE_STRING, "netCDF4", tmpStr);

    } else {
        strcpy(tmpStr, "output file #1 name,\n");
        strcat(tmpStr, "        output vicarious L1B for inverse mode\n");
        strcat(tmpStr, "   ofile[#] = additional output L2 file name");
        option = clo_addOption(list, "ofile", CLO_TYPE_OFILE, "output", tmpStr);
        clo_addOptionAlias(option, "ofile1");

        strcpy(tmpStr, "output file format\n");
        strcat(tmpStr, "        netcdf4: output a netCDF version 4 file\n");
        strcat(tmpStr, "        hdf4:    output a HDF version 4 file");
        clo_addOption(list, "oformat", CLO_TYPE_STRING, "netCDF4", tmpStr);
    }

    strcpy(tmpStr, "output file color depth for HDF4 file\n");
    strcat(tmpStr, "        8bit:  output 8 bit color depth\n");
    strcat(tmpStr, "        24bit: output 24 bit color depth");
    clo_addOption(list, "oformat_depth", CLO_TYPE_STRING, "8bit", tmpStr);

    clo_addOption(list, "deflate", CLO_TYPE_INT, "0", "deflation level");

    strcpy(tmpStr, "input L2 file names for sensor to be\n");
    strcat(tmpStr, "        used as a calibrator.  Alternatively, a data point can be used as a\n");
    strcat(tmpStr, "        calibrator (e.g. MOBY)\n");
    strcat(tmpStr, "   il2file[#] = additional L2 calibration file names");
    option = clo_addOption(list, "il2file", CLO_TYPE_IFILE, NULL, tmpStr);
    clo_addOptionAlias(option, "il2file1");

    clo_addOption(list, "tgtfile", CLO_TYPE_IFILE, NULL, "vicarious calibration target file");
    clo_addOption(list, "aerfile", CLO_TYPE_IFILE, NULL, "aerosol model specification file");
    clo_addOption(list, "metafile", CLO_TYPE_IFILE, NULL, "output meta-data file");

    strcpy(tmpStr, "L2 products to be included in ofile #1\n");
    strcat(tmpStr, "   l2prod[#] = L2 products to be included in ofile[#]");
    option = clo_addOption(list, "l2prod", CLO_TYPE_STRING, NULL, tmpStr);
    clo_addOptionAlias(option, "l2prod1");

    clo_addOption(list, "spixl", CLO_TYPE_INT, "1", "start pixel number");
    clo_addOption(list, "epixl", CLO_TYPE_INT, "-1", "end pixel number (-1=the last pixel)");
    clo_addOption(list, "dpixl", CLO_TYPE_INT, "1", "pixel sub-sampling interval");
    clo_addOption(list, "sline", CLO_TYPE_INT, "1", "start line number");
    clo_addOption(list, "eline", CLO_TYPE_INT, "-1", "end line number (-1=the last line)");
    clo_addOption(list, "dline", CLO_TYPE_INT, "1", "line sub-sampling interval");

    strcpy(tmpStr, " control-point pixel increment for lon/lat\n");
    strcat(tmpStr, "        arrays");
    clo_addOption(list, "ctl_pt_incr", CLO_TYPE_INT, "8", tmpStr);

    strcpy(tmpStr, "toggle ocean processing\n");
    strcat(tmpStr, "        1: On\n");
    strcat(tmpStr, "        0: Off\n");
    strcat(tmpStr, "        2: force all pixels to be processed as ocean");
    clo_addOption(list, "proc_ocean", CLO_TYPE_INT, "1", tmpStr);

    clo_addOption(list, "proc_land", CLO_TYPE_BOOL, "off", "toggle land processing");

    strcpy(tmpStr, "toggle SST processing\n");
    strcat(tmpStr, "        (default=1 for MODIS, 0 otherwise)");
    clo_addOption(list, "proc_sst", CLO_TYPE_BOOL, NULL, tmpStr);

    strcpy(tmpStr, "processing mode\n");
    strcat(tmpStr, "        0: forward processing\n");
    strcat(tmpStr, "        1: inverse (calibration) mode, targeting to nLw=0\n");
    strcat(tmpStr, "        2: inverse (calibration) mode, given nLw target\n");
    strcat(tmpStr, "        3: inverse (calibration) mode, given Lw target (internally normalized)");
    clo_addOption(list, "mode", CLO_TYPE_INT, "0", tmpStr);

    strcpy(tmpStr, "seawater IOP options\n");
    strcat(tmpStr, "        0: static values\n");
    strcat(tmpStr, "        1: temperature & salinity-dependent seawater nw, aw, bbw\n");
    clo_addOption(list, "seawater_opt", CLO_TYPE_INT, "0", tmpStr);

    clo_addOption(list, "atmocor", CLO_TYPE_BOOL, "on", "toggle atmospheric correction");

    sprintf(tmpStr, "aerosol mode option\n");
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s      %3d: No aerosol subtraction\n", tmpStr, AERNULL);
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s      >0: Multi-scattering with fixed model (provide model number, 1-N,\n", tmpStr1);
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s           relative to aermodels list)\n", tmpStr1);
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s      %3d: White aerosol extrapolation.\n", tmpStr1, AERWHITE);
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s      %3d: Multi-scattering with 2-band model selection\n", tmpStr1, AERWANG);
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s      %3d: Multi-scattering with 2-band, RH-based model selection and\n", tmpStr1, AERRHNIR);
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s           iterative NIR correction\n", tmpStr1);
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s      %3d: Multi-scattering with 2-band model selection\n", tmpStr1, AERWANGNIR);
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s           and iterative NIR correction\n", tmpStr1);
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s      %3d: Multi-scattering with fixed model pair\n", tmpStr1, FIXMODPAIR);
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s           (requires aermodmin, aermodmax, aermodrat specification)\n", tmpStr1);
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s      %3d: Multi-scattering with fixed model pair\n", tmpStr1, FIXMODPAIRNIR);
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s           and iterative NIR correction\n", tmpStr1);
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s           (requires aermodmin, aermodmax, aermodrat specification)\n", tmpStr1);
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s      %3d: Multi-scattering with fixed angstrom\n", tmpStr1, FIXANGSTROM);
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s           (requires aer_angstrom specification)\n", tmpStr1);
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s      %3d: Multi-scattering with fixed angstrom\n", tmpStr1, FIXANGSTROMNIR);
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s           and iterative NIR correction\n", tmpStr1);
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s           (requires aer_angstrom specification)\n", tmpStr1);
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s      %3d: Multi-scattering with fixed aerosol optical thickness\n", tmpStr1, FIXAOT);
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s           (requires taua specification)\n", tmpStr1);
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s      %3d: Multi-scattering with 2-band model selection using Wang et al. 2009\n", tmpStr1, AERWANGSWIR);
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s           to switch between SWIR and NIR. (MODIS only, requires aer_swir_short,\n", tmpStr1);
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s           aer_swir_long, aer_wave_short, aer_wave_long)\n", tmpStr1);
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s      %3d: Multi-scattering with MUMM correction\n", tmpStr1, AERMUMM);
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s           and MUMM NIR calculation \n", tmpStr1);
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s      %3d: Spectral optimization via Kuchinke (SeaWiFS-only)\n", tmpStr1, AERSOA);
    strncpy(tmpStr1, tmpStr, 2048);
    sprintf(tmpStr, "%s      %3d: Spectral matching via Gordon (SeaWiFS-only)", tmpStr1, AERSMA);
    clo_addOption(list, "aer_opt", CLO_TYPE_INT, "99", tmpStr);

    clo_addOption(list, "aermodfile", CLO_TYPE_IFILE, NULL, "aerosol model filename leader");

    clo_addOption(list, "aer_wave_short", CLO_TYPE_INT, "765", "shortest sensor wavelength for aerosol\n        model selection");
    clo_addOption(list, "aer_wave_long", CLO_TYPE_INT, "865", "longest sensor wavelength for aerosol\n        model selection");
    clo_addOption(list, "aer_swir_short", CLO_TYPE_INT, "-1", "shortest sensor wavelength for\n        SWIR-based NIR Lw correction");
    clo_addOption(list, "aer_swir_long", CLO_TYPE_INT, "-1", "longest sensor wavelength for SWIR-based\n        NIR Lw correction");
    clo_addOption(list, "aer_rrs_short", CLO_TYPE_FLOAT, "-1.0", "Rrs at shortest sensor wavelength for\n        aerosol model selection");
    clo_addOption(list, "aer_rrs_long", CLO_TYPE_FLOAT, "-1.0", "Rrs at longest sensor wavelength for\n        aerosol model selection");
    clo_addOption(list, "aermodmin", CLO_TYPE_INT, "-1", "lower-bounding model to use for fixed model\n        pair aerosol option");
    clo_addOption(list, "aermodmax", CLO_TYPE_INT, "-1", "upper-bounding model to use for fixed model\n        pair aerosol option");
    clo_addOption(list, "aermodrat", CLO_TYPE_FLOAT, "0.0", "ratio to use for fixed model pair aerosol\n        option");
    clo_addOption(list, "aer_angstrom", CLO_TYPE_FLOAT, "-999.0", "aerosol angstrom exponent for model\n        selection");
    clo_addOption(list, "aer_iter_max", CLO_TYPE_INT, "10", "maximum number of iterations for NIR\n        water-leaving radiance estimation.");
    clo_addOption(list, "mumm_alpha", CLO_TYPE_FLOAT, "1.72", "water-leaving reflectance ratio for MUMM\n       turbid water atmospheric correction");
    clo_addOption(list, "mumm_gamma", CLO_TYPE_FLOAT, "1.0", "two-way Rayleigh-aerosol transmittance\n       ratio for MUMM turbid water atmospheric correction");
    clo_addOption(list, "mumm_epsilon", CLO_TYPE_FLOAT, "1.0", "aerosol reflectance ratio for MUMM\n        turbid water atmospheric correction");

    strcpy(tmpStr, "absorbing aerosol flagging option\n");
    strcat(tmpStr, "       -1: disable\n");
    strcat(tmpStr, "        0: use rhow constant\n");
    strcat(tmpStr, "        1: apply chlorophyll climatology to calculate rhow\n");
    strcat(tmpStr, "        2: 1+validate against nLw_412 climatology");
    clo_addOption(list, "absaer_opt", CLO_TYPE_INT, "-1", tmpStr);

    strcpy(tmpStr, "glint correction:\n");
    strcat(tmpStr, "        0: glint correction off\n");
    strcat(tmpStr, "        1: standard glint correction\n");
    strcat(tmpStr, "        2: simple glint correction");
    clo_addOption(list, "glint_opt", CLO_TYPE_INT, "1", tmpStr);

    strcpy(tmpStr, "out-of-band correction for water-leaving\n");
    strcat(tmpStr, "        radiances\n");
    strcat(tmpStr, "        2: On (default for MODIS, SeaWiFS, OCTS)\n");
    strcat(tmpStr, "        0: Off (default for MOS, OSMI)");
    clo_addOption(list, "outband_opt", CLO_TYPE_INT, "99", tmpStr);

    clo_addOption(list, "cirrus_opt", CLO_TYPE_BOOL, NULL, "cirrus cloud reflectance correction option");

    clo_addOption(list, "oxaband_opt", CLO_TYPE_BOOL, NULL, "oxygen a-band correction option\n        (default On for SeaWiFS, OSMI, and OCTS, Off otherwise)");
    clo_addOption(list, "filter_opt", CLO_TYPE_BOOL, NULL, "filtering input data option\n        (default On for OCTS, Off otherwise)");
    clo_addOption(list, "filter_file", CLO_TYPE_IFILE, "$OCDATAROOT/sensor/sensor_filter.dat", "\n        data file for input filtering");

    strcpy(tmpStr, "Bidirectional reflectance correction\n");
    strcat(tmpStr, "        0: no correction\n");
    strcat(tmpStr, "        1: Fresnel reflection/refraction correction for sensor path\n");
    strcat(tmpStr, "        3: Fresnel reflection/refraction correction for sensor + solar path\n");
    strcat(tmpStr, "        7: Morel f/Q + Fresnel solar + Fresnel sensor\n");
    strcat(tmpStr, "       15: Gordon DT + Morel f/Q + Fresnel solar + Fresnel sensor\n");
    strcat(tmpStr, "       19: Morel Q + Fresnel solar + Fresnel sensor");
    clo_addOption(list, "brdf_opt", CLO_TYPE_INT, "-1", tmpStr);

    strcpy(tmpStr, "gaseous transmittance bitmask selector\n");
    strcat(tmpStr, "        0: no correction\n");
    strcat(tmpStr, "        1: Ozone\n");
    strcat(tmpStr, "        2: CO2\n");
    strcat(tmpStr, "        4: NO2\n");
    strcat(tmpStr, "        8: H2O\n");
    strcat(tmpStr, "       16: Use ATREM (H2O)");
    clo_addOption(list, "gas_opt", CLO_TYPE_INT, "1", tmpStr);

    strcpy(tmpStr, "ATREM gaseous transmittance bitmask selector\n");
    strcat(tmpStr, "        0: H2O only\n");
    strcat(tmpStr, "        1: Ozone\n");
    strcat(tmpStr, "        2: CO2\n");
    strcat(tmpStr, "        4: NO2\n");
    strcat(tmpStr, "        8: CO\n");
    strcat(tmpStr, "       16: CH4\n");
    strcat(tmpStr, "       32: O2\n");
    strcat(tmpStr, "       64: N2O ");
    clo_addOption(list, "atrem_opt", CLO_TYPE_INT, "0", tmpStr);
    strcpy(tmpStr, "ATREM gaseous transmittance geometry option\n");
    strcat(tmpStr, "        0: Only recalculate geometry when error threshold reached (fast)\n");
    strcat(tmpStr, "        1: Recalculate geometry every pixel (slow)");
    clo_addOption(list, "atrem_geom", CLO_TYPE_INT, "0", tmpStr);
    strcpy(tmpStr, "ATREM gaseous transmittance calculation option\n");
    strcat(tmpStr, "        0: Calculate transmittance using k-distribution method (fast)\n");
    strcat(tmpStr, "        1: Calculate transmittance using full method (slow)");
    clo_addOption(list, "atrem_full", CLO_TYPE_INT, "0", tmpStr);
    strcpy(tmpStr, "ATREM gaseous transmittance Atm. model selection\n");
    strcat(tmpStr, "        0: Use pixel's latitude and date to determine model \n");
    strcat(tmpStr, "        1: tropical\n");
    strcat(tmpStr, "        2: mid latitude summer\n");
    strcat(tmpStr, "        3: mid latitude winter\n");
    strcat(tmpStr, "        4: subarctic summer\n");
    strcat(tmpStr, "        5: subarctic winter\n");
    strcat(tmpStr, "        6: US standard 1962");
    clo_addOption(list, "atrem_model", CLO_TYPE_INT, "0", tmpStr);

    strcpy(tmpStr, "IOP model for use in downstream products\n");
    strcat(tmpStr, "        0: None (products requiring a or bb will fail)\n");
    strcat(tmpStr, "        1: Carder\n");
    strcat(tmpStr, "        2: GSM\n");
    strcat(tmpStr, "        3: QAA\n");
    strcat(tmpStr, "        4: PML\n");
    strcat(tmpStr, "        5: NIWA\n");
    strcat(tmpStr, "        6: LAS\n");
    strcat(tmpStr, "        7: GIOP");
    clo_addOption(list, "iop_opt", CLO_TYPE_INT, "0", tmpStr);

    clo_addOption(list, "polfile", CLO_TYPE_IFILE, NULL, "polarization sensitivities filename leader");

    strcpy(tmpStr, "polarization correction (sensor-specific)\n");
    strcat(tmpStr, "        0: no correction\n");
    strcat(tmpStr, "        1: only Rayleigh component is polarized\n");
    strcat(tmpStr, "        2: all radiance polarized like Rayleigh\n");
    strcat(tmpStr, "        3: only Rayleigh and Glint are polarized (MODIS default)\n");
    strcat(tmpStr, "        4: all radiance polarized like Rayleigh + Glint");
    clo_addOption(list, "pol_opt", CLO_TYPE_INT, "-1", tmpStr);

    strcpy(tmpStr, "radiation correction option (sensor-specific)\n");
    strcat(tmpStr, "        0: no correction\n");
    strcat(tmpStr, "        1: apply MERIS Smile correction");
    clo_addOption(list, "rad_opt", CLO_TYPE_INT, "0", tmpStr);

    strcpy(tmpStr, "switch to ocean color virtual constellation\n");
    strcat(tmpStr, "        0: no\n");
    strcat(tmpStr, "        1: NN");

    clo_addOption(list, "ocrvc_opt", CLO_TYPE_INT, "0", tmpStr);

    clo_addOption(list, "xcalfile", CLO_TYPE_IFILE, NULL, "cross-calibration file");

    strcpy(tmpStr, "cross-calibration option (sensor-specific) comma separated\n");
    strcat(tmpStr, "        list of option values, 1 per band, with bands listed in xcal_wave.\n");
    strcat(tmpStr, "        3: apply cross-calibration corrections (polarization and rvs)\n");
    strcat(tmpStr, "        2: apply cross-calibration polarization corrections\n");
    strcat(tmpStr, "        1: apply cross-calibration rvs corrections\n");
    strcat(tmpStr, "        0: no correction");
    clo_addOption(list, "xcal_opt", CLO_TYPE_INT, NULL, tmpStr);
    
	strcpy(tmpStr, "bandshifting option \n");
    strcat(tmpStr, "        1: apply bio-optical bandshift\n");
    strcat(tmpStr, "        0: linear interpolation");
    clo_addOption(list, "band_shift_opt", CLO_TYPE_INT, "0", tmpStr);

	strcpy(tmpStr, "uncertainty simulation option \n");
    strcat(tmpStr, "       >0.0: add random normal noise to gain\n");
    strcat(tmpStr, "        0.0: no random noise added");
    clo_addOption(list, "add_noise_sigma", CLO_TYPE_FLOAT, "0.0", tmpStr);
    
	strcpy(tmpStr, "noise scale factor option \n");
    strcat(tmpStr, "       !=0.0: scales noise by factoring SNR\n");
    strcat(tmpStr, "       =1.0: noise not scaled");
    clo_addOption(list, "noise_scale", CLO_TYPE_FLOAT, "1.0", tmpStr);
	
    strcpy(tmpStr, "single band perturbation option,  \n");
    strcat(tmpStr, "       >=0: wiggle specified band number \n");
    strcat(tmpStr, "        -1: no perturbation");
    clo_addOption(list, "wiggle_band", CLO_TYPE_INT, "0", tmpStr);
    
	strcpy(tmpStr, "noise level for single band perturbation option \n");
    strcat(tmpStr, "       !=0.0: perturb by specified amount\n");
    strcat(tmpStr, "        0.0: no perturbation");
    clo_addOption(list, "wiggle_by", CLO_TYPE_FLOAT, "0.0", tmpStr);
    
    strcpy(tmpStr, "wavelengths at which to apply cross-calibration.  Comma\n");
    strcat(tmpStr, "        separated list of sensor wavelength values associated with xcal_opt.");
    clo_addOption(list, "xcal_wave", CLO_TYPE_FLOAT, NULL, tmpStr);

    strcpy(tmpStr, "processing resolution (MODIS only)\n");
    strcat(tmpStr, "       -1: standard ocean 1km processing\n");
    strcat(tmpStr, "     1000: 1km resolution including aggregated 250 and 500m land bands\n");
    strcat(tmpStr, "      500: 500m resolution including aggregated 250 land bands and\n");
    strcat(tmpStr, "           replication for lower resolution bands\n");
    strcat(tmpStr, "      250: 250m resolution with replication for lower resolution bands");
    clo_addOption(list, "resolution", CLO_TYPE_INT, "-1", tmpStr);

    strcpy(tmpStr, "GIOP model aph function type\n");
    strcat(tmpStr, "        0: tabulated (supplied via giop_aph_file)\n");
    strcat(tmpStr, "        2: Bricaud et al. 1995 (chlorophyll supplied via default empirical algorithm)\n");
    strcat(tmpStr, "        3: Ciotti and Bricaud 2006 (size fraction supplied via giop_aph_s)");
    clo_addOption(list, "giop_aph_opt", CLO_TYPE_INT, "2", tmpStr);
    clo_addOption(list, "giop_aph_file", CLO_TYPE_IFILE, "$OCDATAROOT/common/aph_default.txt", "\n        GIOP model, tabulated aph spectra");
    clo_addOption(list, "giop_aph_s", CLO_TYPE_FLOAT, "-1000.0", "GIOP model, spectral parameter\n        for aph");

    strcpy(tmpStr, "GIOP model adg function type\n");
    strcat(tmpStr, "        0: tabulated (supplied via giop_adg_file)\n");
    strcat(tmpStr, "        1: exponential with exponent supplied via giop_adg_s)\n");
    strcat(tmpStr, "        2: exponential with exponent derived via Lee et al. (2002)\n");
    strcat(tmpStr, "        3: exponential with exponent derived via OBPG method");
    clo_addOption(list, "giop_adg_opt", CLO_TYPE_INT, "1", tmpStr);
    clo_addOption(list, "giop_adg_file", CLO_TYPE_STRING, "$OCDATAROOT/common/adg_default.txt", "\n        GIOP model, tabulated adg spectra");
    clo_addOption(list, "giop_adg_s", CLO_TYPE_FLOAT, "0.018", "GIOP model, spectral parameter\n        for adg");

    strcpy(tmpStr, "GIOP model acdom function type\n");
    strcat(tmpStr, "        0: tabulated (supplied via giop_acdom_file)\n"); 
    strcat(tmpStr, "        1: no data");
    clo_addOption(list, "giop_acdom_opt", CLO_TYPE_INT, "1", tmpStr);
    clo_addOption(list, "giop_acdom_file", CLO_TYPE_IFILE, NULL, "\n        file of specific CDOM absorption coefficients for aLMI");
    
    strcpy(tmpStr, "GIOP model anap function type\n");
    strcat(tmpStr, "        0: tabulated (supplied via giop_anap_file)\n");
    strcat(tmpStr, "        1: no data");
    clo_addOption(list, "giop_anap_opt", CLO_TYPE_INT, "1", tmpStr);
    clo_addOption(list, "giop_anap_file", CLO_TYPE_IFILE, NULL, "\n        file of specific NAP absorption coefficients for aLMI");

    
    strcpy(tmpStr, "GIOP model bbp function type\n");
    strcat(tmpStr, "        0: tabulated (supplied via giop_bbp_file)\n");
    strcat(tmpStr, "        1: power-law with exponent supplied via giop_bbp_s)\n");
    strcat(tmpStr, "        2: power-law with exponent derived via Hoge & Lyon (1996)\n");
    strcat(tmpStr, "        3: power-law with exponent derived via Lee et al. (2002)\n");
    strcat(tmpStr, "        5: power-law with exponent derived via Ciotti et al. (1999)\n");
    strcat(tmpStr, "        6: power-law with exponent derived via Morel & Maritorena (2001)\n");
    strcat(tmpStr, "        7: power-law with exponent derived via Loisel & Stramski (2000)\n");
    strcat(tmpStr, "        8: spectrally independent vector derived via Loisel & Stramski (2000)\n");
    strcat(tmpStr, "        9: fixed vector derived via Loisel & Stramski (2000)\n");
    strcat(tmpStr, "       10: fixed vector derived via lee et al. (2002)");
    clo_addOption(list, "giop_bbp_opt", CLO_TYPE_INT, "3", tmpStr);
    clo_addOption(list, "giop_bbp_file", CLO_TYPE_IFILE, "$OCDATAROOT/common/bbp_default.txt", "\n        GIOP model, tabulated bbp spectra");
    clo_addOption(list, "giop_bbp_s", CLO_TYPE_FLOAT, "-1000.0", "GIOP model, spectral parameter\n        for bbp");

    strcpy(tmpStr, "GIOP model bbph function type\n");
    strcat(tmpStr, "        0: tabulated (supplied via giop_bbph_file)\n");
    strcat(tmpStr, "        1: no data");
    clo_addOption(list, "giop_bbph_opt", CLO_TYPE_INT, "1", tmpStr);
    clo_addOption(list, "giop_bbph_file", CLO_TYPE_IFILE, NULL, "\n        file of specific phytoplankton backscattering coefficients for aLMI");
    
    strcpy(tmpStr, "GIOP model bbnap function type\n");
    strcat(tmpStr, "        0: tabulated (supplied via giop_bbnap_file)\n");
    strcat(tmpStr, "        1: no data");
    clo_addOption(list, "giop_bbnap_opt", CLO_TYPE_INT, "1", tmpStr);
    clo_addOption(list, "giop_bbnap_file", CLO_TYPE_IFILE, NULL, "\n        file of specific nap backscattering coefficients for aLMI");
    
    strcpy(tmpStr, "GIOP model Rrs to bb/(a+bb) method\n");
    strcat(tmpStr, "        0: Gordon quadratic (specified with giop_grd)\n");
    strcat(tmpStr, "        1: Morel f/Q");
    clo_addOption(list, "giop_rrs_opt", CLO_TYPE_INT, "0", tmpStr);
    clo_addOption(list, "giop_grd", CLO_TYPE_FLOAT, "[0.0949,0.0794]", "GIOP model, Gordon\n        Rrs to bb/(a+bb) quadratic coefficients");
    clo_addOption(list, "giop_rrs_diff", CLO_TYPE_FLOAT, "0.33", "GIOP model, maximum difference between input and modeled Rrs");
    strcpy(tmpStr, "GIOP model list of sensor wavelengths for\n");
    strcat(tmpStr, "        optimization comma-separated list, default is all visible bands (400-700nm)");
    clo_addOption(list, "giop_wave", CLO_TYPE_FLOAT, "-1", tmpStr);
    clo_addOption(list, "giop_rrs_unc", CLO_TYPE_FLOAT, "-1", tmpStr);
    clo_addOption(list, "giop_maxiter", CLO_TYPE_INT, "50", "GIOP Model iteration limit");

    strcpy(tmpStr, "GIOP model optimization method\n");
    strcat(tmpStr, "        0: Amoeba optimization\n");
    strcat(tmpStr, "        1: Levenberg-Marquardt optimization\n");
    strcat(tmpStr, "        3: SVD matrix inversion\n");
    strcat(tmpStr, "        4: SIOP adaptive matrix inversion");
    clo_addOption(list, "giop_fit_opt", CLO_TYPE_INT, "1", tmpStr);

    strcpy(tmpStr, "GSM model options\n");
    strcat(tmpStr, "        0: default coefficients\n");
    strcat(tmpStr, "        1: Chesapeake regional coefficients");
    clo_addOption(list, "gsm_opt", CLO_TYPE_INT, "0", tmpStr);

    strcpy(tmpStr, "SM fit algorithm\n");
    strcat(tmpStr, "        0: Amoeba\n");
    strcat(tmpStr, "        1: Levenberg-Marquardt");
    clo_addOption(list, "gsm_fit", CLO_TYPE_INT, "0", tmpStr);

    clo_addOption(list, "gsm_adg_s", CLO_TYPE_FLOAT, "0.02061", "GSM IOP model, spectral slope for adg");
    clo_addOption(list, "gsm_bbp_s", CLO_TYPE_FLOAT, "1.03373", "GSM IOP model, spectral slope for bbp");
    clo_addOption(list, "gsm_aphw", CLO_TYPE_FLOAT, "[412.0, 443.0, 490.0, 510.0, 555.0, 670.0]", "\n        GSM IOP model, wavelengths of ap* table");
    clo_addOption(list, "gsm_aphs", CLO_TYPE_FLOAT, "[0.00665, 0.05582, 0.02055, 0.01910, 0.01015, 0.01424]", "GSM IOP model, coefficients of ap* table");

    option = clo_addOption(list, "qaa_adg_s", CLO_TYPE_FLOAT, "0.015", "QAA IOP model, spectral\n        slope for adg");
    clo_addOptionAlias(option, "qaa_S");

    clo_addOption(list, "qaa_wave", CLO_TYPE_INT, NULL, "sensor wavelengths for QAA algorithm");
    clo_addOption(list, "chloc2_wave", CLO_TYPE_INT, "[-1,-1]", "sensor wavelengths for OC2 chlorophyll\n        algorithm");
    clo_addOption(list, "chloc2_coef", CLO_TYPE_FLOAT, "[0.0,0.0,0.0,0.0,0.0]", "coefficients for OC2\n        chlorophyll algorithm");
    clo_addOption(list, "chloc3_wave", CLO_TYPE_INT, "[-1,-1,-1]", "sensor wavelengths for OC3\n         chlorophyll algorithm");
    clo_addOption(list, "chloc3_coef", CLO_TYPE_FLOAT, "[0.0,0.0,0.0,0.0,0.0]", "coefficients for OC3\n        chlorophyll algorithm");
    clo_addOption(list, "chloc4_wave", CLO_TYPE_INT, "[-1,-1,-1,-1]", "sensor wavelengths for OC4\n        chlorophyll algorithm");
    clo_addOption(list, "chloc4_coef", CLO_TYPE_FLOAT, "[0.0,0.0,0.0,0.0,0.0]", "coefficients for OC4\n        chlorophyll algorithm");
    clo_addOption(list, "kd2_wave", CLO_TYPE_INT, "[-1,-1]", "sensor wavelengths for polynomial Kd(490)\n        algorithm");
    clo_addOption(list, "kd2_coef", CLO_TYPE_FLOAT, "[0.0,0.0,0.0,0.0,0.0,0.0]", "sensor wavelengths\n        for polynomial Kd(490) algorithm");
    clo_addOption(list, "flh_offset", CLO_TYPE_FLOAT, "0.0", "bias to subtract\n        from retrieved fluorescence line height");
    clo_addOption(list, "btfile", CLO_TYPE_IFILE, NULL, "IR brightness temperature file");
    clo_addOption(list, "sstcoeffile", CLO_TYPE_IFILE, NULL, "IR sst algorithm coefficients file");
    clo_addOption(list, "sstssesfile", CLO_TYPE_IFILE, NULL, "IR sst algorithm error statistics file");
    clo_addOption(list, "sst4coeffile", CLO_TYPE_IFILE, NULL, "SWIR sst algorithm coefficients file");
    clo_addOption(list, "sst4ssesfile", CLO_TYPE_IFILE, NULL, "SWIR sst algorithm error statistics file");
    clo_addOption(list, "sst3coeffile", CLO_TYPE_IFILE, NULL, "Triple window sst algorithm coefficients file");
    clo_addOption(list, "sst3ssesfile", CLO_TYPE_IFILE, NULL, "Triple window sst algorithm error statistics file");
    clo_addOption(list, "vcnnfile", CLO_TYPE_IFILE, NULL, "virtual constellation neural net file");
    clo_addOption(list, "picfile", CLO_TYPE_IFILE, NULL, "pic table for Balch 2-band algorithm");
    clo_addOption(list, "owtfile", CLO_TYPE_IFILE, NULL, "optical water type file");
    clo_addOption(list, "owtchlerrfile", CLO_TYPE_IFILE, NULL, "chl error file associate with optical water type");

    strcpy(tmpStr, "[");
    strcat(tmpStr, defaermodels[0]);
    for (i = 1; i < defnaermodels; i++) {
        strcat(tmpStr, ",");
        strcat(tmpStr, defaermodels[i]);
    }
    strcat(tmpStr, "]");
    clo_addOption(list, "aermodels", CLO_TYPE_STRING, tmpStr, "aerosol models");

    clo_addOption(list, "met1", CLO_TYPE_IFILE, "$OCDATAROOT/common/met_climatology.hdf", "\n        1st meteorological ancillary data file");
    clo_addOption(list, "met2", CLO_TYPE_IFILE, NULL, "2nd meteorological ancillary data file");
    clo_addOption(list, "met3", CLO_TYPE_IFILE, NULL, "3rd meteorological ancillary data file");
    clo_addOption(list, "ozone1", CLO_TYPE_IFILE, "$OCDATAROOT/common/ozone_climatology.hdf", "\n        1st ozone ancillary data file");
    clo_addOption(list, "ozone2", CLO_TYPE_IFILE, NULL, "2nd ozone ancillary data file");
    clo_addOption(list, "ozone3", CLO_TYPE_IFILE, NULL, "3rd ozone ancillary data file");
    clo_addOption(list, "anc_cor_file", CLO_TYPE_IFILE, NULL, "ozone correction file");
    clo_addOption(list, "land", CLO_TYPE_IFILE, "$OCDATAROOT/common/landmask.dat", "land mask file");
    clo_addOption(list, "water", CLO_TYPE_IFILE, "$OCDATAROOT/common/watermask.dat", "\n        shallow water mask file");
    clo_addOption(list, "demfile", CLO_TYPE_IFILE, "$OCDATAROOT/common/digital_elevation_map.hdf", "\n        digital elevation map file");
    clo_addOption(list, "elevfile", CLO_TYPE_IFILE, "$OCDATAROOT/common/ETOPO1_ocssw.nc", "global elevation netCDF file");
    clo_addOption(list, "elev_auxfile", CLO_TYPE_IFILE, NULL, "auxiliary elevation netCDF file");
    clo_addOption(list, "icefile", CLO_TYPE_IFILE, "$OCDATAROOT/common/ice_mask.hdf", "sea ice file");
    clo_addOption(list, "ice_threshold", CLO_TYPE_FLOAT, "0.1", "sea ice fraction above which will be\n        flagged as sea ice");
    clo_addOption(list, "sstfile", CLO_TYPE_IFILE, "$OCDATAROOT/common/sst_climatology.hdf", "input\n        SST reference file");

    strcpy(tmpStr, "Reference SST field source\n");
    strcat(tmpStr, "        0: Reynolds OI SST reference file\n");
    strcat(tmpStr, "        1: AMSR-E daily SST reference file\n");
    strcat(tmpStr, "        2: AMSR-E 3-day SST reference file\n");
    strcat(tmpStr, "        3: ATSR monthly SST reference file\n");
    strcat(tmpStr, "        4: NTEV2 monthly SST reference file\n");
    strcat(tmpStr, "        5: AMSR-E 3-day or night SST reference file\n");
    strcat(tmpStr, "        6: WindSat daily SST reference file\n");
    strcat(tmpStr, "        7: WindSat 3-day SST reference file\n");
    strcat(tmpStr, "        8: WindSat 3-day or night SST reference file\n");
    clo_addOption(list, "sstreftype", CLO_TYPE_INT, "0", tmpStr);

    clo_addOption(list, "sssfile", CLO_TYPE_IFILE, "$OCDATAROOT/common/sss_climatology_woa2009.hdf", "input\n        SSS reference file");
    clo_addOption(list, "no2file", CLO_TYPE_IFILE, "$OCDATAROOT/common/no2_climatology.hdf", "no2\n        ancillary file");
    clo_addOption(list, "alphafile", CLO_TYPE_IFILE, "$OCDATAROOT/common/alpha510_climatology.hdf", "\n        alpha510 climatology file");
    clo_addOption(list, "tauafile", CLO_TYPE_IFILE, "$OCDATAROOT/common/taua865_climatology.hdf", "\n        taua865 climatology file");
    clo_addOption(list, "cldfile", CLO_TYPE_IFILE, NULL, "cloud mask file (MODIS only)");
    clo_addOption(list, "calfile", CLO_TYPE_IFILE, NULL, "system calibration file");
    clo_addOption(list, "fqfile", CLO_TYPE_IFILE, "$OCDATAROOT/common/morel_fq.nc", "f/Q correction file");
    clo_addOption(list, "parfile", CLO_TYPE_IFILE, NULL, "par climatology file for OPP calculation");
    clo_addOption(list, "offset", CLO_TYPE_FLOAT, NULL, "calibration offset adjustment");
    clo_addOption(list, "gain", CLO_TYPE_FLOAT, NULL, "calibration gain multiplier");
    clo_addOption(list, "gain_unc", CLO_TYPE_FLOAT, NULL, "calibration gain uncertainty");
    clo_addOption(list, "flaguse", CLO_TYPE_STRING, default_flaguse, "Flags to use");
    clo_addOption(list, "xcalbox", CLO_TYPE_INT, "0", "pixel size of the central box in the L1 scene\n        (e.g. 5 pixels around MOBY) to be extracted into xcalfile for the\n        cross-calibration, 0=whole L1");
    clo_addOption(list, "xcalboxcenter", CLO_TYPE_INT, "[0,0]", "Central [ipix, iscan] of the box in\n        the L1 scene, [0,0] = center of the L1 scene");
    clo_addOption(list, "xcalpervalid", CLO_TYPE_INT, "0", "min percent of valid cross-calibration\n        pixels within the box or the L1 scene, 0 = at least 1 pixel");
    clo_addOption(list, "xcalsubsmpl", CLO_TYPE_INT, "1", "Sub-sampling rate for the data to be used\n        for the cross-calibration");


    strcpy(tmpStr, "max distance between L1 and L2 pixels\n");
    strcat(tmpStr, "       in km\n");
    strcat(tmpStr, "       -1.0: use mean res of L1 data\n");
    strcat(tmpStr, "        0.0: max{mean(L1 res),mean(L2 res)}");
    clo_addOption(list, "maxpointdist", CLO_TYPE_FLOAT, "0.0", tmpStr);

    sprintf(tmpStr2, "%f", CHL_MAX);
    sprintf(tmpStr, "threshold on L2 data chlorophyll\n        (%f=CHL_MAX)",
            CHL_MAX);
    clo_addOption(list, "chlthreshold", CLO_TYPE_FLOAT, tmpStr2, tmpStr);

    sprintf(tmpStr2, "%f", AOT_MAX);
    sprintf(tmpStr, "threshold on L2 data AOTs\n        (%f=AOT_MAX)", AOT_MAX);
    clo_addOption(list, "aotthreshold", CLO_TYPE_FLOAT, tmpStr2, tmpStr);

    strcpy(tmpStr, "\n        coccolithophore algorithm coefs");
    clo_addOption(list, "coccolith", CLO_TYPE_FLOAT, "[1.1,0.9,0.75,1.85,1.0,1.65,0.6,1.15]", tmpStr);

    clo_addOption(list, "cirrus_thresh", CLO_TYPE_FLOAT, "[-1.0,-1.0]", "cirrus reflectance thresholds");
    clo_addOption(list, "taua", CLO_TYPE_FLOAT, NULL, "[taua_band1,...,taua_bandn] aerosol optical thickness of the\n        calibration data point");
    option = clo_addOption(list, "cloud_thresh", CLO_TYPE_FLOAT, "0.027", "cloud reflectance\n        threshold");
    clo_addOptionAlias(option, "albedo");
    clo_addOption(list, "cloud_wave", CLO_TYPE_FLOAT, "865.0", "wavelength of cloud reflectance test");
    clo_addOption(list, "cloud_eps", CLO_TYPE_FLOAT, "-1.0", "cloud reflectance ratio threshold\n        (-1.0=disabled)");
    option = clo_addOption(list, "glint_thresh", CLO_TYPE_FLOAT, "0.005", "high sun glint threshold");
    clo_addOptionAlias(option, "glint");
    clo_addOption(list, "absaer", CLO_TYPE_FLOAT, "0.0", "absorbing aerosol threshold on aerosol index");
    clo_addOption(list, "rhoamin", CLO_TYPE_FLOAT, "0.0001", "min NIR aerosol reflectance to attempt\n        model lookup");
    clo_addOption(list, "sunzen", CLO_TYPE_FLOAT, "75.0", "sun zenith angle threshold in deg.");
    clo_addOption(list, "satzen", CLO_TYPE_FLOAT, "60.0", "satellite zenith angle threshold");
    clo_addOption(list, "epsmin", CLO_TYPE_FLOAT, "0.85", "minimum epsilon to trigger atmospheric\n        correction failure flag");
    clo_addOption(list, "epsmax", CLO_TYPE_FLOAT, "1.35", "maximum epsilon to trigger atmospheric\n         correction failure flag");
    clo_addOption(list, "tauamax", CLO_TYPE_FLOAT, "0.3", "maximum 865 aerosol optical depth to trigger\n        hitau flag");
    clo_addOption(list, "nLwmin", CLO_TYPE_FLOAT, "0.15", "minimum nLw(555) to trigger low Lw flag");
    clo_addOption(list, "hipol", CLO_TYPE_FLOAT, "0.5", "threshold on degree-of-polarization to set\n        HIPOL flag");
    clo_addOption(list, "wsmax", CLO_TYPE_FLOAT, "8.0", "windspeed limit on white-cap correction in m/s");
    clo_addOption(list, "windspeed", CLO_TYPE_FLOAT, "-1000.0", "user over-ride of windspeed in m/s\n        (-1000=use ancillary files)");
    clo_addOption(list, "windangle", CLO_TYPE_FLOAT, "-1000.0", "user over-ride of wind angle in deg\n        (-1000=use ancillary files)");
    clo_addOption(list, "pressure", CLO_TYPE_FLOAT, "-1000.0", "user over-ride of atmospheric pressure\n        in mb (-1000=use ancillary files)");
    clo_addOption(list, "ozone", CLO_TYPE_FLOAT, "-1000.0", "user over-ride of ozone concentration in\n        cm (-1000=use ancillary files)");
    clo_addOption(list, "relhumid", CLO_TYPE_FLOAT, "-1000.0", "user over-ride of relative humidity in\n        percent (-1000=use ancillary files)");
    clo_addOption(list, "watervapor", CLO_TYPE_FLOAT, "-1000.0", "user over-ride of water vapor in\n        g/cm^2 (-1000=use ancillary files)");
    clo_addOption(list, "maskland", CLO_TYPE_BOOL, "on", "land mask option");
    clo_addOption(list, "maskbath", CLO_TYPE_BOOL, "off", "shallow water mask option");
    clo_addOption(list, "maskcloud", CLO_TYPE_BOOL, "on", "cloud mask option");
    clo_addOption(list, "maskglint", CLO_TYPE_BOOL, "off", "glint mask option");
    clo_addOption(list, "masksunzen", CLO_TYPE_BOOL, "off", "large sun zenith angle mask option");
    clo_addOption(list, "masksatzen", CLO_TYPE_BOOL, "off", "large satellite zenith angle mask option");
    clo_addOption(list, "maskhilt", CLO_TYPE_BOOL, "on", "high Lt masking");
    clo_addOption(list, "maskstlight", CLO_TYPE_BOOL, "on", "stray light masking");
    clo_addOption(list, "sl_frac", CLO_TYPE_FLOAT, "0.25", "SeaWiFS only, straylight fractional\n        threshold on Ltypical");
    clo_addOption(list, "sl_pixl", CLO_TYPE_INT, "-1", "SeaWiFS only, number of LAC pixels for\n        straylight flagging");
    clo_addOption(list, "vcal_opt", CLO_TYPE_INT, "-1", "Vicarious calibration option");
    clo_addOption(list, "vcal_chl", CLO_TYPE_FLOAT, "-1.0", "Vicarious calibration chl");
    clo_addOption(list, "vcal_solz", CLO_TYPE_FLOAT, "-1.0", "Vicarious calibration solz");
    clo_addOption(list, "vcal_nLw", CLO_TYPE_FLOAT, NULL, "Vicarious calibration normalized water leaving radiances");
    clo_addOption(list, "vcal_Lw", CLO_TYPE_FLOAT, NULL, "Vicarious calibration water leaving radiances");

    strcpy(tmpStr, "depth to use to exclude data from target file\n");
    strcat(tmpStr, "   e.g. -1000 excludes depths less than 1000m\n");
    clo_addOption(list, "vcal_depth", CLO_TYPE_FLOAT, "-1000.0", tmpStr);

    strcpy(tmpStr, "minimum # of samples in a bin for acceptance");
    clo_addOption(list, "vcal_min_nbin", CLO_TYPE_INT, "4", tmpStr);

    strcpy(tmpStr, "minimum # of scenes in a bin for acceptance");
    clo_addOption(list, "vcal_min_nscene", CLO_TYPE_INT, "3", tmpStr);

    clo_addOption(list, "owmcfile", CLO_TYPE_IFILE, "$OCDATAROOT/common/owmc_lut.hdf", "lut for OWMC\n        classification");

    clo_addOption(list, "stype", CLO_TYPE_INT, "0", "scaling type\n        0: log\n        1: linear");
    clo_addOption(list, "rgb", CLO_TYPE_INT, "[1,1,1]", "bands to use for red, green and blue");
    clo_addOption(list, "north", CLO_TYPE_FLOAT, "-999", "north boundary");
    clo_addOption(list, "south", CLO_TYPE_FLOAT, "-999", "south boundary");
    clo_addOption(list, "east", CLO_TYPE_FLOAT, "-999", "east boundary");
    clo_addOption(list, "west", CLO_TYPE_FLOAT, "-999", "west boundary");
    clo_addOption(list, "xbox", CLO_TYPE_INT, "-1", "number of pixels on either side of the SW point");
    if (!strcmp(prog, "l1det2det")) {
        clo_addOption(list, "ybox", CLO_TYPE_INT, "-1", "number of scan lines to require for valid detector runs");
    }else {
        clo_addOption(list, "ybox", CLO_TYPE_INT, "-1", "number of scan lines on either side of the SW point");

    }

    clo_addOption(list, "width", CLO_TYPE_INT, "600", "width of output image");
    clo_addOption(list, "threshold", CLO_TYPE_FLOAT, "0.1", "threshold for the number of good pixels\n        before an image is produced");
    clo_addOption(list, "datamin", CLO_TYPE_FLOAT, "0.01", "minimum reflectance for scaling");
    clo_addOption(list, "datamax", CLO_TYPE_FLOAT, "0.9", "maximum reflectance for scaling");
    clo_addOption(list, "subsamp", CLO_TYPE_INT, "1", "sub-sampling interval");

    clo_addOption(list, "viirsnv7", CLO_TYPE_INT, "-1", "=1 to use the VIIRSN V7 high senz latband sst and sst3 equations");
    clo_addOption(list, "viirsnosisaf", CLO_TYPE_INT, "0", "=1 to use the VIIRSN OSI-SAF sst and sst3 equations");
    clo_addOption(list, "sstrefdif", CLO_TYPE_FLOAT, "100.0", "stricter sst-ref difference threshold");
    clo_addOption(list, "newavhrrcal", CLO_TYPE_INT, "0", "=1 for new noaa-16 calibration");
    strcpy(tmpStr, "\n        Channel 22 detector corrections (MODIS only)");
    clo_addOption(list, "ch22detcor", CLO_TYPE_FLOAT, "[1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]", tmpStr);
    strcpy(tmpStr, "\n        Channel 23 detector corrections (MODIS only)");
    clo_addOption(list, "ch23detcor", CLO_TYPE_FLOAT, "[1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]", tmpStr);

    clo_addOption(list, "prodxmlfile", CLO_TYPE_OFILE, NULL, "output XML file describing all possible products");
    clo_addOption(list, "breflectfile", CLO_TYPE_IFILE, NULL, "input NetCDF file for bottom reflectances and bottom types");

    strcpy(tmpStr, "\nThis program produces a PPM-formatted output image rendered in a Plate Carree\n");
    strcat(tmpStr, "projection.\n\n");
    strcat(tmpStr, "The default band combination produces a \"true-color\" image. Other combinations\n");
    strcat(tmpStr, "may be chosen with the \"rgb=\" option.  The expected argument to this option is\n");
    strcat(tmpStr, "a comma separated string of wavelengths that specifies the desired bands in\n");
    strcat(tmpStr, "red-green-blue order.  For example, to produce a false color SeaWiFS output\n");
    strcat(tmpStr, "image using 865, 670 and 555 nm as the red, green, and blue values\n");
    strcat(tmpStr, "respectively, the option would be specified as \"rgb=865,670,555\".\n");
    clo_addOption(list, "help_true_color", CLO_TYPE_HELP, "0", tmpStr);
    
    strcpy(tmpStr, "Raman scattering Rrs corection options\n");
    strcat(tmpStr, "        0: no correction \n");
    strcat(tmpStr, "        1: Lee et al. (2013) empirical correction \n");
    strcat(tmpStr, "        2: Westberry et al. (2013) analytical correction \n");
    strcat(tmpStr, "        3: Lee et al. (1994) analytical correction \n");
    clo_addOption(list, "raman_opt", CLO_TYPE_INT, "0", tmpStr);

    return 0;
}

//-----------------------------------------------------------------------

/* 
   Read the command line option and all of the default parameter files.

   This is the order for loading the options:
    - read the command line to get the ifile and suite options
    - load the main program defaults file
    - load the sensor specific defaults file
    - load the sub-sensor defaults file
    - load the suite file
    - load the command line (including specified par files)
    - re-load the command line so they take precedence

 */
int l2gen_read_options(clo_optionList_t* list, char* progName,
        int argc, char* argv[], filehandle *l1file) {
    char *dataRoot;
    char tmpStr[FILENAME_MAX];
    char *ifile;
    clo_option_t* option;
    char localSuite[FILENAME_MAX];
    int i;
    char l2genProgName[] = "msl12";

    assert(list);

    if ((dataRoot = getenv("OCDATAROOT")) == NULL) {
        printf("-E- OCDATAROOT environment variable is not defined.\n");
        return (-1);
    }

    // disable the dump option until we have read all of the files
    clo_setEnableDumpOptions(0);

    // read command line args to get the ifile parameter
    enableFileDecending = 1;
    clo_readArgs(list, argc, argv);
    // if a file was decended, make sure the command args over-rides
    enableFileDecending = 0;
    clo_readArgs(list, argc, argv);
    enableFileDecending = 1;

    // see if suite param was set
    localSuite[0] = '\0';
    if (clo_isSet(list, "suite")) {
        strcpy(localSuite, clo_getString(list, "suite"));
    } // suite option was set

    ifile = clo_getString(list, "ifile");
    i = clo_getInt(list, "resolution");
    if (i == -1000)
        l1file->modis_subset_compat = 1;
    i = clo_getInt(list, "ocrvc_opt");
    l1file->ocrvc_opt = i;
    strcpy(l1file->name, ifile);
    if (getFormat(l1file) == -1) {
        printf("-E- %s Line %d: Could not find type for file %s.\n",
                __FILE__, __LINE__, ifile);
        return (-1);
    }

    // load l2gen program defaults
    sprintf(tmpStr, "%s/common/%s_defaults.par", dataRoot, l2genProgName);
    if (want_verbose)
        printf("Loading default parameters from %s\n", tmpStr);
    clo_readFile(list, tmpStr);

    // load non-l2gen program defaults file
    if (progName && progName[0] && strcmp(progName, l2genProgName)) {
      sprintf(tmpStr, "%s/common/%s_defaults.par", dataRoot, progName);
      if( access( tmpStr, F_OK ) != -1 ) {
        // file exists
        if (want_verbose)
	  printf("Loading default parameters from %s\n", tmpStr);
        clo_readFile(list, tmpStr);
      }
    }
    
    // load the sensor specific defaults file
    sprintf(tmpStr, "%s/%s/%s_defaults.par", dataRoot,
            sensorDir[l1file->sensorID], l2genProgName);
    if (want_verbose)
        printf("Loading default parameters for %s from %s\n",
                sensorName[l1file->sensorID], tmpStr);
    clo_readFile(list, tmpStr);

    // load the sub-sensor specific defaults file
    if (l1file->subsensorID >= 0) {
        sprintf(tmpStr, "%s/%s/%s_defaults_%s.par", dataRoot,
                sensorDir[l1file->sensorID], l2genProgName,
                sensorSub[l1file->subsensorID]);
        if (want_verbose)
            printf("Loading default sub-sensor parameters for %s from %s\n",
                    sensorName[l1file->sensorID], tmpStr);
        clo_readFile(list, tmpStr);
    } // if sub-sensor

    // if suite not set on command line or param file then use the default
    if (localSuite[0] == '\0')
        strcpy(localSuite, clo_getString(list, "suite"));

    // load the suite file
    sprintf(tmpStr, "%s/%s/%s_defaults_%s.par", dataRoot,
            sensorDir[l1file->sensorID], l2genProgName, localSuite);
    if( access( tmpStr, F_OK ) != -1 ) {
        // file exists
        if (want_verbose)
            printf("Loading parameters for suite %s from %s\n", localSuite, tmpStr);
        clo_readFile(list, tmpStr);
    }
    // load the program specific defaults file, if exists and is not "msl12"
    if (progName && progName[0] && strcmp(progName, l2genProgName)) {

        sprintf(tmpStr, "%s/%s/%s_defaults.par", dataRoot,
                sensorDir[l1file->sensorID], progName);
        if( access( tmpStr, F_OK ) != -1 ) {
            if (want_verbose)
                printf("Loading default parameters for %s from %s\n",
                        sensorName[l1file->sensorID], tmpStr);
            clo_readFile(list, tmpStr);
        }
    }

    // set the default l2prod lists before the command line or par file
    // is loaded
    option = clo_findOption(list, "l2prod");
    if (option && clo_isOptionSet(option))
        strcpy(default_l2prod[0], option->valStr);
    for (i = 0; i < MAX_OFILES; i++) {
        sprintf(tmpStr, "l2prod%d", i + 1);
        option = clo_findOption(list, tmpStr);
        if (option && clo_isOptionSet(option))
            strcpy(default_l2prod[i], option->valStr);
        else
            default_l2prod[i][0] = '\0';
    }

    // re-load the command line and par file
    if (want_verbose)
        printf("Loading command line parameters\n\n");
    enableFileDecending = 1;
    clo_readArgs(list, argc, argv);

    // enable the dump option the last time through
    clo_setEnableDumpOptions(1);

    // disable file decending so the command line options override
    // anything that is in a parameter file
    enableFileDecending = 0;
    clo_readArgs(list, argc, argv);

    return 0;
}

//-----------------------------------------------------------------------------

int l2gen_load_input(clo_optionList_t *list, instr *input) {
    char str_buf[FILENAME_MAX];
    char str_buf2[FILENAME_MAX];
    char tmp_file[FILENAME_MAX];
    char *strVal;
    clo_option_t *option;
    int numOptions;
    int optionId;
    char keyword[FILENAME_MAX];
    int count;
    char **strArray;
    float *fArray;
    int *iArray;
    int i, j;
    FILE *fp;
    const char* tmpStr;


    // first load up the default_l2prod
    for (i = 0; i < MAX_OFILES; i++) {
        if (default_l2prod[i][0]) {
            strcpy(input->def_l2prod[i], default_l2prod[i]);
        }
    }

    /* allocate and initialize dynamic arrays in input struc */

    msl12_input_nbands_init(input);

    numOptions = clo_getNumOptions(list);
    for (optionId = 0; optionId < numOptions; optionId++) {
        option = clo_getOption(list, optionId);

        // ignore options of type CLO_TYPE_HELP 
        if(option->dataType == CLO_TYPE_HELP)
            continue;

        strcpy(keyword, option->key);

        /* change keyword to lower case */
        strVal = keyword;
        while (*strVal != '\0') {
            *strVal = tolower(*strVal);
            strVal++;
        }

        if (strcmp(keyword, "-help") == 0)
            ;
        else if (strcmp(keyword, "-version") == 0)
            ;
        else if (strncmp(keyword, "-dump_options", 13) == 0)
            ;
        else if (strcmp(keyword, "par") == 0)
            ;
        else if (strcmp(keyword, "ifile") == 0) {
            strVal = clo_getOptionString(option);
            parse_file_name(strVal, tmp_file);
            strcpy(input->ifile[0], tmp_file);

        } else if (strncmp(keyword, "ifile", 5) == 0) {
            for (i = 1; i < MAX_IFILES; i++) {
                sprintf(str_buf, "ifile%d", i + 1);
                if (strcmp(keyword, str_buf) == 0) {
                    if (i == 0)
                        break;
                    strVal = clo_getOptionString(option);
                    parse_file_name(strVal, tmp_file);
                    strcpy(input->ifile[i], tmp_file);
                    break;
                }
            }
            if (i >= MAX_IFILES) {
                printf("-E- l2gen_load_input: %s bigger than MAX_IFILES(%d)\n",
                        keyword, MAX_IFILES);
                return -1;
            }

        } else if (strcmp(keyword, "ilist") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                if ((fp = fopen(tmp_file, "r")) == NULL) {
                    printf("Error: Unable to open input file list %s\n", tmp_file);
                    return -1;
                }
                tmp_file[0] = '\x0';
                for (i = 0; i < MAX_IFILES; i++) {
                    if (fscanf(fp, "%s\n", tmp_file) != EOF) {
                        strcpy((input->ifile)[i], tmp_file);
                        sprintf(str_buf, "ifile%d", i + 1);
                        sprintf(str_buf2, "ilist=%s", strVal);
                        clo_setString(list, str_buf, tmp_file, str_buf2);
                    } else
                        break;
                }
                fclose(fp);
            }

        } else if (strcmp(keyword, "calfile") == 0) {
            strVal = clo_getOptionString(option);
            parse_file_name(strVal, tmp_file);
            strcpy(input->calfile, tmp_file);

        } else if (strcmp(keyword, "fqfile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->fqfile, tmp_file);
            }
        } else if (strcmp(keyword, "parfile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->parfile, tmp_file);
            }
        } else if (strcmp(keyword, "cldfile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->cldfile, tmp_file);
            }

        } else if (strcmp(keyword, "ofile") == 0) {
            strVal = clo_getOptionString(option);
            parse_file_name(strVal, tmp_file);
            strcpy(input->ofile[0], tmp_file);

        } else if (strncmp(keyword, "ofile", 5) == 0) {
            for (i = 1; i < MAX_OFILES; i++) {
                sprintf(str_buf, "ofile%d", i + 1);
                if (strcmp(keyword, str_buf) == 0) {
                    strVal = clo_getOptionString(option);
                    parse_file_name(strVal, tmp_file);
                    strcpy(input->ofile[i], tmp_file);
                    break;
                }
            }
            if (i >= MAX_OFILES) {
                printf("-E- l2gen_load_input: %s bigger than MAX_OFILES(%d)\n",
                        keyword, MAX_OFILES);
                return -1;
            }

        } else if (strcmp(keyword, "oformat") == 0) {
            strVal = clo_getOptionString(option);
            tmpStr = getFileFormatName(strVal);
            if(tmpStr == NULL) {
                printf("-E- l2gen_load_input: oformat=%s is not a recognized file format\n",
                        strVal);
                return -1;
            }
            strcpy(input->oformat, tmpStr);

        } else if (strcmp(keyword, "oformat_depth") == 0) {
            strVal = clo_getOptionString(option);
            if(strcasecmp(strVal, "8bit") == 0 || strcasecmp(strVal, "24bit") == 0) {
                strcpy(input->oformat_depth, strVal);
            } else {
                printf("-E- l2gen_load_input: oformat_depth=%s is not a valid color depth\n",
                        strVal);
                return -1;
            }

        } else if (strcmp(keyword, "il2file") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->il2file[0], tmp_file);
            }

        } else if (strncmp(keyword, "il2file", 7) == 0) {
            for (i = 1; i < MAX_OFILES; i++) {
                sprintf(str_buf, "il2file%d", i + 1);
                if (strcmp(keyword, str_buf) == 0) {
                    if (i == 0)
                        break;
                    strVal = clo_getOptionString(option);
                    parse_file_name(strVal, tmp_file);
                    strcpy(input->il2file[i], tmp_file);
                    break;
                }
            }
            if (i >= MAX_OFILES) {
                printf("-E- l2gen_load_input: %s bigger than MAX_OFILES(%d)\n",
                        keyword, MAX_OFILES);
                return -1;
            }

        } else if (strcmp(keyword, "l2prod") == 0) {
            strArray = clo_getOptionStrings(option, &count);
            input->l2prod[0][0] = '\0';
            for (i = 0; i < count; i++) {
                if (i != 0)
                    strcat(input->l2prod[0], " ");
                strcat(input->l2prod[0], strArray[i]);
            }

        } else if (strncmp(keyword, "l2prod", 6) == 0) {
            for (i = 1; i < MAX_OFILES; i++) {
                sprintf(str_buf, "l2prod%d", i + 1);
                if (strcmp(keyword, str_buf) == 0) {
                    if (i == 0)
                        break;
                    strArray = clo_getOptionStrings(option, &count);
                    input->l2prod[i][0] = '\0';
                    for (j = 0; j < count; j++) {
                        if (j != 0)
                            strcat(input->l2prod[i], " ");
                        strcat(input->l2prod[i], strArray[j]);
                    }
                    break;
                }
            }
            if (i >= MAX_OFILES) {
                printf("-E- l2gen_load_input: %s bigger than MAX_OFILES(%d)\n",
                        keyword, MAX_OFILES);
                return -1;
            }

        } else if (strcmp(keyword, "pversion") == 0) {
            // if I copy it this way then spaces are allowed in the string
            if (option->valStr)
                strcpy(input->pversion, option->valStr);

        } else if (strcmp(keyword, "suite") == 0) {
            strcpy(input->suite, clo_getOptionString(option));

        } else if (strcmp(keyword, "eval") == 0) {
            input->evalmask = clo_getOptionInt(option);

        } else if (strcmp(keyword, "mode") == 0) {
            input->mode = clo_getOptionInt(option);

        } else if (strcmp(keyword, "deflate") == 0) {
            input->deflate = clo_getOptionInt(option);

        } else if (strcmp(keyword, "spixl") == 0) {
            input->spixl = clo_getOptionInt(option);

        } else if (strcmp(keyword, "epixl") == 0) {
            input->epixl = clo_getOptionInt(option);

        } else if (strcmp(keyword, "dpixl") == 0) {
            input->dpixl = clo_getOptionInt(option);

        } else if (strcmp(keyword, "sline") == 0) {
            input->sline = clo_getOptionInt(option);

        } else if (strcmp(keyword, "eline") == 0) {
            input->eline = clo_getOptionInt(option);

        } else if (strcmp(keyword, "dline") == 0) {
            input->dline = clo_getOptionInt(option);

        } else if (strcmp(keyword, "ctl_pt_incr") == 0) {
            input->ctl_pt_incr = clo_getOptionInt(option);

        } else if (strcmp(keyword, "proc_ocean") == 0) {
            input->proc_ocean = clo_getOptionInt(option);

        } else if (strcmp(keyword, "proc_land") == 0) {
            input->proc_land = clo_getOptionBool(option);

        } else if (strcmp(keyword, "proc_sst") == 0) {
            input->proc_sst = clo_getOptionBool(option);

        } else if (strcmp(keyword, "atmocor") == 0) {
            input->atmocor = clo_getOptionBool(option);

        } else if (strcmp(keyword, "aer_opt") == 0) {
            input->aer_opt = clo_getOptionInt(option);

        } else if (strcmp(keyword, "aer_wave_short") == 0) {
            input->aer_wave_short = clo_getOptionInt(option);

        } else if (strcmp(keyword, "aer_wave_long") == 0) {
            input->aer_wave_long = clo_getOptionInt(option);

        } else if (strcmp(keyword, "aer_swir_short") == 0) {
            input->aer_swir_short = clo_getOptionInt(option);

        } else if (strcmp(keyword, "aer_swir_long") == 0) {
            input->aer_swir_long = clo_getOptionInt(option);

        } else if (strcmp(keyword, "aer_rrs_short") == 0) {
            input->aer_rrs_short = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "aer_rrs_long") == 0) {
            input->aer_rrs_long = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "aer_angstrom") == 0) {
            input->aer_angstrom = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "aer_iter_max") == 0) {
            input->aer_iter_max = clo_getOptionInt(option);

        } else if (strcmp(keyword, "seawater_opt") == 0) {
            input->seawater_opt = clo_getOptionInt(option);

        } else if (strcmp(keyword, "brdf_opt") == 0) {
            input->brdf_opt = clo_getOptionInt(option);

        } else if (strcmp(keyword, "gas_opt") == 0) {
            input->gas_opt = clo_getOptionInt(option);

        } else if (strcmp(keyword, "atrem_opt") == 0) {
            input->atrem_opt = clo_getOptionInt(option);

        } else if (strcmp(keyword, "atrem_full") == 0) {
            input->atrem_full = clo_getOptionInt(option);
            if (option->count > 0) {
               if (input->atrem_full <=0) input->atrem_full = 0; else input->atrem_full = 1;}

        } else if (strcmp(keyword, "atrem_geom") == 0) {
            input->atrem_geom = clo_getOptionInt(option);
            if (option->count > 0) {
               if (input->atrem_geom <= 0) input->atrem_geom = 0; else input->atrem_geom = 1;}

        } else if (strcmp(keyword, "atrem_model") == 0) {
            input->atrem_model = clo_getOptionInt(option);

        } else if (strcmp(keyword, "gsm_opt") == 0) {
            input->gsm_opt = clo_getOptionInt(option);

        } else if (strcmp(keyword, "gsm_fit") == 0) {
            input->gsm_fit = clo_getOptionInt(option);

        } else if (strcmp(keyword, "gsm_adg_s") == 0) {
            input->gsm_adg_s = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "gsm_bbp_s") == 0) {
            input->gsm_bbp_s = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "gsm_aphs") == 0) {
            fArray = clo_getOptionFloats(option, &count);
            if (count > input->nbands) {
                printf("-E- number of gsm_aphs elements (%d) must be %d or less\n", count, input->nbands);
                exit(1);
            }

            for (i = 0; i < count; i++)
                input->gsm_aphs[i] = fArray[i];

        } else if (strcmp(keyword, "gsm_aphw") == 0) {
            fArray = clo_getOptionFloats(option, &count);
            if (count > input->nbands) {
                printf("-E- number of gsm_aphw elements (%d) must be %d or less\n", count, input->nbands);
                exit(1);
            }

            for (i = 0; i < count; i++)
                input->gsm_aphw[i] = fArray[i];

        } else if (strcmp(keyword, "qaa_adg_s") == 0) {
            input->qaa_adg_s = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "qaa_wave") == 0) {
            if (clo_isOptionSet(option)) {
                iArray = clo_getOptionInts(option, &count);
                if (count != 5) {
                    printf("-E- number of qaa_wave elements must be 5.\n");
                    exit(1);
                }

                for (i = 0; i < count; i++)
                    input->qaa_wave[i] = iArray[i];
            }

        } else if (strcmp(keyword, "giop_maxiter") == 0) {
            input->giop_maxiter = clo_getOptionInt(option);

        } else if (strcmp(keyword, "giop_fit_opt") == 0) {
            input->giop_fit_opt = clo_getOptionInt(option);

        } else if (strcmp(keyword, "giop_aph_opt") == 0) {
            input->giop_aph_opt = clo_getOptionInt(option);

        } else if (strcmp(keyword, "giop_adg_opt") == 0) {
            input->giop_adg_opt = clo_getOptionInt(option);

        } else if (strcmp(keyword, "giop_acdom_opt") == 0) {
            input->giop_acdom_opt = clo_getOptionInt(option);

        } else if (strcmp(keyword, "giop_anap_opt") == 0) {
            input->giop_anap_opt = clo_getOptionInt(option);
            
        } else if (strcmp(keyword, "giop_bbp_opt") == 0) {
            input->giop_bbp_opt = clo_getOptionInt(option);

        } else if (strcmp(keyword, "giop_bbph_opt") == 0) {
            input->giop_bbph_opt = clo_getOptionInt(option);
            
        } else if (strcmp(keyword, "giop_bbnap_opt") == 0) {
            input->giop_bbnap_opt = clo_getOptionInt(option);

        } else if (strcmp(keyword, "giop_rrs_opt") == 0) {
            input->giop_rrs_opt = clo_getOptionInt(option);

        } else if (strcmp(keyword, "giop_rrs_diff") == 0) {
            input->giop_rrs_diff = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "giop_aph_file") == 0) {
            strVal = clo_getOptionString(option);
            parse_file_name(strVal, tmp_file);
            strcpy(input->giop_aph_file, tmp_file);

        } else if (strcmp(keyword, "giop_aph_s") == 0) {
            input->giop_aph_s = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "giop_adg_file") == 0) {
            strVal = clo_getOptionString(option);
            parse_file_name(strVal, tmp_file);
            strcpy(input->giop_adg_file, tmp_file);

        } else if (strcmp(keyword, "giop_adg_s") == 0) {
            input->giop_adg_s = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "giop_bbp_file") == 0) {
            strVal = clo_getOptionString(option);
            parse_file_name(strVal, tmp_file);
            strcpy(input->giop_bbp_file, tmp_file);

        } else if (strcmp(keyword, "giop_bbp_s") == 0) {
            input->giop_bbp_s = clo_getOptionFloat(option);
            
        } else if (strcmp(keyword, "giop_acdom_file") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->giop_acdom_file, tmp_file);
            }
                
        } else if (strcmp(keyword, "giop_anap_file") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->giop_anap_file, tmp_file);
            }
            
        } else if (strcmp(keyword, "giop_bbph_file") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->giop_bbph_file, tmp_file);
            }   
            
        } else if (strcmp(keyword, "giop_bbnap_file") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->giop_bbnap_file, tmp_file);
            }

        } else if (strcmp(keyword, "giop_grd") == 0) {
            fArray = clo_getOptionFloats(option, &count);
            if (count > 2) {
                printf("-E- number of giop_grd elements must be 2 or less.\n");
                exit(1);
            }
            for (i = 0; i < count; i++)
                input->giop_grd[i] = fArray[i];

        } else if (strcmp(keyword, "giop_wave") == 0) {
            if (clo_isOptionSet(option)) {
                fArray = clo_getOptionFloats(option, &count);
                if (count > input->nbands) {
                    printf("-E- number of giop_wave elements (%d) must be %d or less\n", count, input->nbands);
                    exit(1);
                }

                for (i = 0; i < count; i++)
                    input->giop_wave[i] = fArray[i];
            }


        } else if (strcmp(keyword, "giop_rrs_unc") == 0) {
            if (clo_isOptionSet(option)) {
                fArray = clo_getOptionFloats(option, &count);
                if (count > input->nbands) {
                    printf("-E- number of giop_rrs_unc elements (%d) must be %d or less\n", count, input->nbands);
                    exit(1);
                }

                for (i = 0; i < count; i++)
                    input->giop_rrs_unc[i] = fArray[i];
            }

        } else if (strcmp(keyword, "iop_opt") == 0) {
            input->iop_opt = clo_getOptionInt(option);

        } else if (strcmp(keyword, "xcalfile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->xcalfile, tmp_file);
            }

        } else if (strcmp(keyword, "xcal_opt") == 0) {
            if (clo_isOptionSet(option)) {
                iArray = clo_getOptionInts(option, &count);
                input->xcal_nwave = count;
                if (count > input->nbands) {
                    printf("-E- number of xcal_opt elements (%d) must be %d or less\n", count, input->nbands);
                    exit(1);
                }

                for (i = 0; i < count; i++)
                    input->xcal_opt[i] = iArray[i];
            }

        } else if (strcmp(keyword, "xcal_wave") == 0) {
            if (clo_isOptionSet(option)) {
                fArray = clo_getOptionFloats(option, &count);
                if (count > input->nbands) {
                    printf("-E- number of xcal_wave elements (%d) must be %d or less\n", count, input->nbands);
                    exit(1);
                }

                input->xcal_nwave = count;

                for (i = 0; i < count; i++)
                    input->xcal_wave[i] = fArray[i];
            }

        } else if (strcmp(keyword, "band_shift_opt") == 0) {
            input->band_shift_opt = clo_getOptionInt(option);
	
        } else if (strcmp(keyword, "add_noise_sigma") == 0) {
            input->add_noise_sigma = clo_getOptionFloat(option);
        } else if (strcmp(keyword, "noise_scale") == 0) {
            input->noise_scale = clo_getOptionFloat(option);
        } else if (strcmp(keyword, "wiggle_band") == 0){
            input->wiggle_band = clo_getOptionInt(option);
        } else if (strcmp(keyword, "wiggle_by") == 0){
            input->wiggle_by = clo_getOptionFloat(option);
        } else if (strcmp(keyword, "aermodfile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->aermodfile, tmp_file);
            }

        } else if (strcmp(keyword, "polfile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->polfile, tmp_file);
            }

        } else if (strcmp(keyword, "pol_opt") == 0) {
            input->pol_opt = clo_getOptionInt(option);

        } else if (strcmp(keyword, "rad_opt") == 0) {
            input->rad_opt = clo_getOptionInt(option);

        } else if (strcmp(keyword, "ocrvc_opt") == 0) {
            input->ocrvc_opt = clo_getOptionInt(option);

        } else if (strcmp(keyword, "absaer_opt") == 0) {
            input->absaer_opt = clo_getOptionInt(option);

        } else if (strcmp(keyword, "sl_pixl") == 0) {
            input->sl_pixl = clo_getOptionInt(option);

        } else if (strcmp(keyword, "sl_frac") == 0) {
            input->sl_frac = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "resolution") == 0) {
            input->resolution = clo_getOptionInt(option);

        } else if (strcmp(keyword, "glint_opt") == 0) {
            input->glint_opt = clo_getOptionInt(option);

        } else if (strcmp(keyword, "outband_opt") == 0) {
            input->outband_opt = clo_getOptionInt(option);

        } else if (strcmp(keyword, "cirrus_opt") == 0) {
            input->cirrus_opt = clo_getOptionBool(option);

        } else if (strcmp(keyword, "oxaband_opt") == 0) {
            input->oxaband_opt = clo_getOptionBool(option);

        } else if (strcmp(keyword, "filter_opt") == 0) {
            if (clo_isOptionSet(option))
                input->filter_opt = clo_getOptionBool(option);

        } else if (strcmp(keyword, "filter_file") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->filter_file, tmp_file);
            }

        } else if (strcmp(keyword, "tgtfile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->tgtfile, tmp_file);
            }

        } else if (strcmp(keyword, "aerfile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->aerfile, tmp_file);
            }

        } else if (strcmp(keyword, "btfile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->btfile, tmp_file);
            }

        } else if (strcmp(keyword, "sstcoeffile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->sstcoeffile, tmp_file);
            }

        } else if (strcmp(keyword, "sstssesfile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->sstssesfile, tmp_file);
            }

        } else if (strcmp(keyword, "sst4coeffile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->sst4coeffile, tmp_file);
            }

        } else if (strcmp(keyword, "sst4ssesfile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->sst4ssesfile, tmp_file);
            }

        } else if (strcmp(keyword, "vcnnfile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->vcnnfile, tmp_file);
            }

        } else if (strcmp(keyword, "picfile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->picfile, tmp_file);
            }

        } else if (strcmp(keyword, "sst3coeffile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->sst3coeffile, tmp_file);
            }

        } else if (strcmp(keyword, "sst3ssesfile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
            strcpy(input->sst3ssesfile, tmp_file);
        }

        } else if (strcmp(keyword, "owtfile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->owtfile, tmp_file);
            }

        } else if (strcmp(keyword, "owtchlerrfile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->owtchlerrfile, tmp_file);
            }

        } else if (strcmp(keyword, "metafile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->metafile, tmp_file);
            }

        } else if (strcmp(keyword, "aerbinfile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->aerbinfile, tmp_file);
            }

        } else if (strcmp(keyword, "met1") == 0) {
            strVal = clo_getOptionString(option);
            parse_file_name(strVal, tmp_file);
            strcpy(input->met1, tmp_file);

        } else if (strcmp(keyword, "met2") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->met2, tmp_file);
            }

        } else if (strcmp(keyword, "met3") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->met3, tmp_file);
            }

        } else if (strcmp(keyword, "ozone1") == 0) {
            strVal = clo_getOptionString(option);
            parse_file_name(strVal, tmp_file);
            strcpy(input->ozone1, tmp_file);

        } else if (strcmp(keyword, "ozone2") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->ozone2, tmp_file);
            }

        } else if (strcmp(keyword, "ozone3") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->ozone3, tmp_file);
            }

        } else if (strcmp(keyword, "anc_cor_file") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->anc_cor_file, tmp_file);
            }

        } else if (strcmp(keyword, "land") == 0) {
            strVal = clo_getOptionString(option);
            parse_file_name(strVal, tmp_file);
            strcpy(input->land, tmp_file);

        } else if (strcmp(keyword, "water") == 0) {
            strVal = clo_getOptionString(option);
            parse_file_name(strVal, tmp_file);
            strcpy(input->water, tmp_file);

        } else if (strcmp(keyword, "demfile") == 0) {
            strVal = clo_getOptionString(option);
            parse_file_name(strVal, tmp_file);
            strcpy(input->demfile, tmp_file);

        } else if (strcmp(keyword, "elevfile") == 0) {
            strVal = clo_getOptionString(option);
            parse_file_name(strVal, tmp_file);
            strcpy(input->elevfile, tmp_file);

        } else if (strcmp(keyword, "elev_auxfile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->elev_auxfile, tmp_file);
            }

        } else if (strcmp(keyword, "icefile") == 0) {
            strVal = clo_getOptionString(option);
            parse_file_name(strVal, tmp_file);
            strcpy(input->icefile, tmp_file);


        } else if (strcmp(keyword, "sstfile") == 0) {
            strVal = clo_getOptionString(option);
            parse_file_name(strVal, tmp_file);
            strcpy(input->sstfile, tmp_file);
        } else if (strcmp(keyword, "sstreftype") == 0) {
            input->sstreftype = clo_getOptionInt(option);

        } else if (strcmp(keyword, "sssfile") == 0) {
            strVal = clo_getOptionString(option);
            parse_file_name(strVal, tmp_file);
            strcpy(input->sssfile, tmp_file);

        } else if (strcmp(keyword, "no2file") == 0) {
            strVal = clo_getOptionString(option);
            parse_file_name(strVal, tmp_file);
            strcpy(input->no2file, tmp_file);

        } else if (strcmp(keyword, "alphafile") == 0) {
            strVal = clo_getOptionString(option);
            parse_file_name(strVal, tmp_file);
            strcpy(input->alphafile, tmp_file);

        } else if (strcmp(keyword, "tauafile") == 0) {
            strVal = clo_getOptionString(option);
            parse_file_name(strVal, tmp_file);
            strcpy(input->tauafile, tmp_file);

        } else if (strcmp(keyword, "geofile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->geofile, tmp_file);
            }

        } else if (strcmp(keyword, "viirscalparfile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->viirscalparfile, tmp_file);
	    }

        } else if (strcmp(keyword, "owmcfile") == 0) {
            strVal = clo_getOptionString(option);
            parse_file_name(strVal, tmp_file);
            strcpy(input->owmcfile, tmp_file);

        } else if (strcmp(keyword, "prodxmlfile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->prodXMLfile, tmp_file);
            }

        } else if (strcmp(keyword, "breflectfile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->breflectfile, tmp_file);
            }

        } else if (strcmp(keyword, "gain") == 0) {
            if (clo_isOptionSet(option)) {
                fArray = clo_getOptionFloats(option, &count);
                if (count > input->nbands) {
                    printf("-E- number of gain elements (%d) must be %d or less\n", count, input->nbands);
                    exit(1);
                }

                for (i = 0; i < count; i++)
                    input->gain[i] = fArray[i];
                numGains = count;
            }

        } else if (strcmp(keyword, "gain_unc") == 0) {
            if (clo_isOptionSet(option)) {
                fArray = clo_getOptionFloats(option, &count);
                if (count > input->nbands) {
                    printf("-E- number of gain_unc elements (%d) must be %d or less\n", count, input->nbands);
                    exit(1);
                }


                for (i = 0; i < count; i++)
                    input->gain_unc[i] = fArray[i];
                numGainUncs = count;
            }

        } else if (strcmp(keyword, "offset") == 0) {
            if (clo_isOptionSet(option)) {
                fArray = clo_getOptionFloats(option, &count);
                if (count > input->nbands) {
                    printf("-E- number of offset elements (%d) must be %d or less\n", count, input->nbands);
                    exit(1);
                }


                for (i = 0; i < count; i++)
                    input->offset[i] = fArray[i];
                numOffsets = count;
            }

        } else if (strcmp(keyword, "flaguse") == 0) {
            strArray = clo_getOptionStrings(option, &count);
            input->flaguse[0] = '\0';
            for (i = 0; i < count; i++) {
                if (i != 0)
                    strcat(input->flaguse, ",");
                strcat(input->flaguse, strArray[i]);
            }

        } else if (strcmp(keyword, "xcalbox") == 0) {
            input->xcalbox = clo_getOptionInt(option);

        } else if (strcmp(keyword, "xcalboxcenter") == 0) {
            iArray = clo_getOptionInts(option, &count);
            if (count > 2) {
                printf("-E- number of xcalboxcenter elements must be 2 or less\n");
                exit(1);
            }
            for (i = 0; i < count; i++)
                input->xcalboxcenter[i] = iArray[i];

        } else if (strcmp(keyword, "xcalpervalid") == 0) {
            input->xcalpervalid = clo_getOptionInt(option);

        } else if (strcmp(keyword, "xcalsubsmpl") == 0) {
            input->xcalsubsmpl = clo_getOptionInt(option);

        } else if (strcmp(keyword, "maxpointdist") == 0) {
            input->maxpointdist = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "chlthreshold") == 0) {
            input->chlthreshold = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "aotthreshold") == 0) {
            input->aotthreshold = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "coccolith") == 0) {
            fArray = clo_getOptionFloats(option, &count);
            if (count > 8) {
                printf("-E- number of coccolith elements must be 8 or less\n");
                exit(1);
            }
            for (i = 0; i < count; i++)
                input->coccolith[i] = fArray[i];

        } else if (strcmp(keyword, "cirrus_thresh") == 0) {
            fArray = clo_getOptionFloats(option, &count);
            if (count > 2) {
                printf("-E- number of cirrus_thresh elements must be 2 or less\n");
                exit(1);
            }
            for (i = 0; i < count; i++)
                input->cirrus_thresh[i] = fArray[i];

        } else if (strcmp(keyword, "chloc2_wave") == 0) {
            iArray = clo_getOptionInts(option, &count);
            if (count > 2) {
                printf("-E- number of chloc2_wave elements must be 2 or less\n");
                exit(1);
            }
            for (i = 0; i < count; i++)
                input->chloc2w[i] = iArray[i];

        } else if (strcmp(keyword, "chloc2_coef") == 0) {
            fArray = clo_getOptionFloats(option, &count);
            if (count > 5) {
                printf("-E- number of chloc2_coef elements must be 5 or less\n");
                exit(1);
            }
            for (i = 0; i < count; i++)
                input->chloc2c[i] = fArray[i];

        } else if (strcmp(keyword, "chloc3_wave") == 0) {
            iArray = clo_getOptionInts(option, &count);
            if (count > 3) {
                printf("-E- number of chloc3_wave elements must be 3 or less\n");
                exit(1);
            }
            for (i = 0; i < count; i++)
                input->chloc3w[i] = iArray[i];

        } else if (strcmp(keyword, "chloc3_coef") == 0) {
            fArray = clo_getOptionFloats(option, &count);
            if (count > 5) {
                printf("-E- number of chloc3_coef elements must be 5 or less\n");
                exit(1);
            }
            for (i = 0; i < count; i++)
                input->chloc3c[i] = fArray[i];

        } else if (strcmp(keyword, "chloc4_wave") == 0) {
            iArray = clo_getOptionInts(option, &count);
            if (count > 4) {
                printf("-E- number of chloc4_wave elements must be 4 or less\n");
                exit(1);
            }
            for (i = 0; i < count; i++)
                input->chloc4w[i] = iArray[i];

        } else if (strcmp(keyword, "chloc4_coef") == 0) {
            fArray = clo_getOptionFloats(option, &count);
            if (count > 5) {
                printf("-E- number of chloc4_coef elements must be 5 or less\n");
                exit(1);
            }
            for (i = 0; i < count; i++)
                input->chloc4c[i] = fArray[i];

        } else if (strcmp(keyword, "kd2_wave") == 0) {
            iArray = clo_getOptionInts(option, &count);
            if (count > 2) {
                printf("-E- number of kd2_wave elements must be 2 or less\n");
                exit(1);
            }
            for (i = 0; i < count; i++)
                input->kd2w[i] = iArray[i];

        } else if (strcmp(keyword, "kd2_coef") == 0) {
            fArray = clo_getOptionFloats(option, &count);
            if (count > 6) {
                printf("-E- number of kd2_coef elements must be 6 or less\n");
                exit(1);
            }
            for (i = 0; i < count; i++)
                input->kd2c[i] = fArray[i];

        } else if (strcmp(keyword, "flh_offset") == 0) {
            input->flh_offset = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "aermodels") == 0) {
            strArray = clo_getOptionStrings(option, &count);
            if (count > MAXAERMOD) {
                printf("-E- number of aermodels must be %d or less\n",
                        MAXAERMOD);
                exit(1);
            }
            for (i = 0; i < count; i++)
                strcpy(input->aermodels[i], strArray[i]);
            for (; i < MAXAERMOD; i++)
                input->aermodels[i][0] = '\0';
            input->naermodels = count;

        } else if (strcmp(keyword, "taua") == 0) {
            if (clo_isOptionSet(option)) {
                fArray = clo_getOptionFloats(option, &count);
                if (count > input->nbands) {
                    printf("-E- number of taua elements (%d) must be %d or less\n", count, input->nbands);
                    exit(1);
                }

                for (i = 0; i < count; i++)
                    input->taua[i] = fArray[i];
                numTauas = count;
            }

        } else if (strcmp(keyword, "aermodrat") == 0) {
            input->aermodrat = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "aermodmin") == 0) {
            input->aermodmin = clo_getOptionInt(option);

        } else if (strcmp(keyword, "aermodmax") == 0) {
            input->aermodmax = clo_getOptionInt(option);

        } else if (strcmp(keyword, "absaer") == 0) {
            input->absaer = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "rhoamin") == 0) {
            input->rhoamin = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "cloud_eps") == 0) {
            input->cloud_eps = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "cloud_thresh") == 0) {
            input->albedo = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "cloud_wave") == 0) {
            input->cloud_wave = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "glint_thresh") == 0) {
            input->glint = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "sunzen") == 0) {
            input->sunzen = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "satzen") == 0) {
            input->satzen = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "epsmin") == 0) {
            input->epsmin = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "epsmax") == 0) {
            input->epsmax = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "tauamax") == 0) {
            input->tauamax = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "nlwmin") == 0) {
            input->nlwmin = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "hipol") == 0) {
            input->hipol = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "wsmax") == 0) {
            input->wsmax = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "windspeed") == 0) {
            input->windspeed = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "windangle") == 0) {
            input->windangle = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "pressure") == 0) {
            input->pressure = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "ozone") == 0) {
            input->ozone = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "relhumid") == 0) {
            input->relhumid = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "watervapor") == 0) {
            input->watervapor = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "ice_threshold") == 0) {
            input->ice_threshold = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "maskland") == 0) {
            input->landmask = clo_getOptionBool(option);

        } else if (strcmp(keyword, "maskbath") == 0) {
            input->bathmask = clo_getOptionBool(option);

        } else if (strcmp(keyword, "maskcloud") == 0) {
            input->cloudmask = clo_getOptionBool(option);

        } else if (strcmp(keyword, "maskglint") == 0) {
            input->glintmask = clo_getOptionBool(option);

        } else if (strcmp(keyword, "masksunzen") == 0) {
            input->sunzenmask = clo_getOptionBool(option);

        } else if (strcmp(keyword, "masksatzen") == 0) {
            input->satzenmask = clo_getOptionBool(option);

        } else if (strcmp(keyword, "maskhilt") == 0) {
            input->hiltmask = clo_getOptionBool(option);

        } else if (strcmp(keyword, "maskstlight") == 0) {
            input->stlightmask = clo_getOptionBool(option);

        } else if (strcmp(keyword, "mumm_alpha") == 0) {
            input->mumm_alpha = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "mumm_gamma") == 0) {
            input->mumm_gamma = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "mumm_epsilon") == 0) {
            input->mumm_epsilon = clo_getOptionFloat(option);


        } else if (strcmp(keyword, "viirsnv7") == 0) {
            input->viirsnv7 = clo_getOptionInt(option);

        } else if (strcmp(keyword, "viirsnosisaf") == 0) {
            input->viirsnosisaf = clo_getOptionInt(option);

        } else if (strcmp(keyword, "sstrefdif") == 0) {
            input->sstrefdif = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "newavhrrcal") == 0) {
            input->newavhrrcal = clo_getOptionInt(option);

        } else if (strcmp(keyword, "ch22detcor") == 0) {
            fArray = clo_getOptionFloats(option, &count);
            if (count != 10) {
                printf("-E- number of ch22detcor elements must be 10 \n");
                exit(1);
            }
            for (i = 0; i < count; i++)
                input->ch22detcor[i] = fArray[i];

        } else if (strcmp(keyword, "ch23detcor") == 0) {
            fArray = clo_getOptionFloats(option, &count);
            if (count != 10) {
                printf("-E- number of ch23detcor elements must be 10 \n");
                exit(1);
            }
            for (i = 0; i < count; i++)
                input->ch23detcor[i] = fArray[i];

        } else if (strcmp(keyword, "vcal_opt") == 0) {
            input->vcal_opt = clo_getOptionInt(option);

        } else if (strcmp(keyword, "vcal_nlw") == 0) {
            if (clo_isOptionSet(option)) {
                fArray = clo_getOptionFloats(option, &count);
                if (count > input->nbands) {
                    printf("-E- number of vcal_nlw elements (%d) must be %d or less\n", count, input->nbands);
                    exit(1);
                }

                for (i = 0; i < count; i++)
                    input->vcal_nLw[i] = fArray[i];
                if (input->vcal_opt < 0)
                    input->vcal_opt = INVERSE_NLW;
            }

        } else if (strcmp(keyword, "vcal_lw") == 0) {
            if (clo_isOptionSet(option)) {
                fArray = clo_getOptionFloats(option, &count);
                if (count > input->nbands) {
                    printf("-E- number of vcal_lw elements (%d) must be %d or less\n", count, input->nbands);
                    exit(1);
                }


                for (i = 0; i < count; i++)
                    input->vcal_Lw[i] = fArray[i];
                if (input->vcal_opt < 0)
                    input->vcal_opt = INVERSE_LW;
            }

        } else if (strcmp(keyword, "vcal_chl") == 0) {
            input->vcal_chl = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "vcal_solz") == 0) {
            input->vcal_solz = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "vcal_depth") == 0) {
            input->vcal_depth = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "vcal_min_nbin") == 0) {
            input->vcal_min_nbin = clo_getOptionInt(option);

        } else if (strcmp(keyword, "vcal_min_nscene") == 0) {
            input->vcal_min_nscene = clo_getOptionInt(option);

        }else if (strcmp(keyword, "stype") == 0) {
            input->stype = clo_getOptionInt(option);

        } else if (strcmp(keyword, "datamin") == 0) {
            input->datamin = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "datamax") == 0) {
            input->datamax = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "west") == 0) {
            input->west = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "east") == 0) {
            input->east = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "north") == 0) {
            input->north = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "south") == 0) {
            input->south = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "width") == 0) {
            input->width = clo_getOptionInt(option);

        } else if (strcmp(keyword, "threshold") == 0) {
            input->threshold = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "rgb") == 0) {
            if (clo_isOptionSet(option)) {
                iArray = clo_getOptionInts(option, &count);
                if (count != 3) {
                    printf("-E- number of rgb elements must be 3\n");
                    exit(1);
                }
                for (i = 0; i < count; i++)
                    input->rgb[i] = iArray[i];
            }

        } else if (strcmp(keyword, "subsamp") == 0) {
            input->subsamp = clo_getOptionInt(option);

        } else if (strcmp(keyword, "xbox") == 0) {
            input->xbox = clo_getOptionInt(option);

        } else if (strcmp(keyword, "ybox") == 0) {
            input->ybox = clo_getOptionInt(option);
            
        } else if (strcmp(keyword, "raman_opt") == 0) {
            input->raman_opt = clo_getOptionInt(option); 

        } else {
            printf("-E- Invalid argument \"%s\"\n", keyword);
            clo_dumpOption(option);
            exit(1);
        }

    } // for optionIDs

    return 0;
}

int msl12_option_input(int argc, char **argv, clo_optionList_t* list,
        char *progName, instr *input, filehandle *l1file) {
    int i;
    char *tmp_str;
    char str_buf[FILENAME_MAX] = "";

    /* For reading sensor table */
    int32_t numBands;

    /* read the command line options, default, sensor, sub-sensor 
       and suite par files.  set elements of l1file */
    if (l2gen_read_options(list, progName, argc, argv, l1file) != 0) {
        printf("-E- %s: Error reading program options.\n", __FILE__);
        return (-1);
    }

    // add more info into the input structure
    strcpy(input->program_name, progName);
    input->sensorID = l1file->sensorID;
    input->subsensorID = l1file->subsensorID;
    input->format = l1file->format;
    input->modis_subset_compat = l1file->modis_subset_compat;

    if (want_verbose)
        printf("Loading user parameters for %s\n\n", sensorName[l1file->sensorID]);

    /*                                                                  */
    /* Now, loop through command arguments again and update input struct*/
    /*                                                                  */
    strcpy(input->pro_control, basename(argv[0]));
    for (i = 1; i < argc; i++) {
        strcat(input->pro_control, " ");
        strcat(input->pro_control, argv[i]);
    }

    
    /* Make sure band-dependent inputs have values for all bands */
    if ((numBands = rdsensorinfo(l1file->sensorID, input->evalmask, NULL, NULL)) < 0) {
        fprintf(stderr, "-E- %s Line %d:  Error reading sensor table.\n",
                __FILE__, __LINE__);
        return (-1);
    }

    input->nbands = numBands;

    if (l2gen_load_input(list, input) != 0) {
        printf("-E- %s: Error loading options into input structure.\n", __FILE__);
        return (-1);
    }

    l1file->geofile = input->geofile;
    l1file->calfile = input->calfile;
    l1file->viirscalparfile = input->viirscalparfile;
    l1file->input = input;
    if (input->resolution == -1000)
        input->resolution = 1000;

    if ((tmp_str = getenv("OCDATAROOT")) == NULL) {
        printf("OCDATAROOT environment variable is not defined.\n");
        return (1);
    }

    if (want_verbose && (input->deflate > 0))
        printf("Internal data compression requested at compression level: %d\n",input->deflate);

    /* Load filter list, if filtering requested */
    if (input->filter_opt == 1) {
        if (input->filter_file[0] == '\0') {
            sprintf(input->filter_file, "%s/%s/%s_filter.dat", tmp_str,
                    sensorDir[l1file->sensorID], sensorDir[l1file->sensorID]);
        }
        rdfilter(input->filter_file, &(input->fctl), numBands);
    }


    if (numGains > 0 && numGains != numBands) {
        fprintf(stderr, "Parameter input error:  Number of input gains must equal %d.\n",
                numBands);
        return (-1);
    }
    if (numGainUncs > 0 && numGainUncs != numBands) {
        fprintf(stderr, "Parameter input error:  Number of input gain uncertainties must equal %d.\n",
                numBands);
        return (-1);
    }
    if (numOffsets > 0 && numOffsets != numBands) {
        fprintf(stderr, "Parameter input error:  Number of input offsets must equal %d.\n",
                numBands);
        return (-1);
    }
    if (numTauas > 0 && numTauas != numBands) {
        fprintf(stderr, "Parameter input error:  Number of input taua values must equal %d.\n",
                numBands);
        return (-1);
    }


    /* Echo calibration */
/*
    if (want_verbose) {
        printf("\nCalibration to be applied:\n");
        printf("Calfile:  %s\n", input->calfile);

        for (i = 0; i < numBands; i++) {
            printf("Vicarious Gain[%d]=%f +/-%f, Offset[%d]=%f\n", i + 1, input->gain[i], input->gain_unc[i], i + 1, input->offset[i]);
        }
        printf("\n\n");

    } // want_verbose
*/
    /*                                                                  */
    /* Build string of input parameters for meta-data documentation     */
    /*                                                                  */

    strcat(input->input_parms, "\n");
    for (i = 0; i < MAX_IFILES; i++) {

        if (input->ifile[i][0] != '\0') {
            if(i==0)
                sprintf(str_buf, "ifile = %s ", input->ifile[i]);
            else
                sprintf(str_buf, "ifile%d = %s ", i + 1, input->ifile[i]);
            strcat(input->input_parms, str_buf);
            strcat(input->input_parms, "\n");
        } else break;
    }

    for (i = 0; i < MAX_OFILES; i++) {

        if (input->ofile[i][0] != '\0') {
            if(i==0)
                sprintf(str_buf, "ofile = %s", input->ofile[i]);
            else
                sprintf(str_buf, "ofile%d = %s", i + 1, input->ofile[i]);
            strcat(input->input_parms, str_buf);
            strcat(input->input_parms, "\n");
        }

        if (input->l2prod[i][0] != '\0') {
            if(i==0)
                sprintf(str_buf, "l2prod = %s", input->l2prod[i]);
            else
                sprintf(str_buf, "l2prod%d = %s", i + 1, input->l2prod[i]);
            strcat(input->input_parms, str_buf);
            strcat(input->input_parms, "\n");
        }

        if (input->mode != FORWARD && input->il2file[i][0] != '\0') {
            sprintf(str_buf, "il2file%d = %s", i + 1, input->il2file[i]);
            strcat(input->input_parms, str_buf);
            strcat(input->input_parms, "\n");
        }

    }

    sprintf(str_buf, "oformat = %s", input->oformat);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "oformat_depth = %s", input->oformat_depth);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "calfile = %s", input->calfile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "fqfile = %s", input->fqfile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "parfile = %s", input->parfile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "cldfile = %s", input->cldfile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "geofile = %s", input->geofile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "viirscalparfile = %s", input->viirscalparfile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "metafile = %s", input->metafile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    if (input->mode != FORWARD) {

        sprintf(str_buf, "flaguse = %s", input->flaguse);
        strcat(input->input_parms, str_buf);
        strcat(input->input_parms, "\n");

        sprintf(str_buf, "maxpointdist = %8.3f", input->maxpointdist);
        strcat(input->input_parms, str_buf);
        strcat(input->input_parms, "\n");

        sprintf(str_buf, "chlthreshold = %8.3f", input->chlthreshold);
        strcat(input->input_parms, str_buf);
        strcat(input->input_parms, "\n");

        sprintf(str_buf, "aotthreshold = %8.3f", input->aotthreshold);
        strcat(input->input_parms, str_buf);
        strcat(input->input_parms, "\n");

        sprintf(str_buf, "xcalbox = %d", (int) input->xcalbox);
        strcat(input->input_parms, str_buf);
        strcat(input->input_parms, "\n");

        sprintf(str_buf, "xcalboxcenter = %d, %d",
                (int) input->xcalboxcenter[0], (int) input->xcalboxcenter[1]);
        strcat(input->input_parms, str_buf);
        strcat(input->input_parms, "\n");

        sprintf(str_buf, "xcalpervalid = %d", (int) input->xcalpervalid);
        strcat(input->input_parms, str_buf);
        strcat(input->input_parms, "\n");

        sprintf(str_buf, "xcalsubsmpl = %d", (int) input->xcalsubsmpl);
        strcat(input->input_parms, str_buf);
        strcat(input->input_parms, "\n");

    }

    sprintf(str_buf, "pversion = %s", input->pversion);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "suite = %s", input->suite);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "eval = %5d", input->evalmask);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "mode = %5d", input->mode);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "deflate = %5d", input->deflate);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "spixl = %5d", input->spixl);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "epixl = %5d", input->epixl);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "dpixl = %5d", input->dpixl);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "sline = %5d", input->sline);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "eline = %5d", input->eline);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "dline = %5d", input->dline);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "ctl_pt_incr = %5d", input->ctl_pt_incr);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "proc_ocean = %3d", input->proc_ocean);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "proc_land = %3d", input->proc_land);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "proc_sst = %3d", input->proc_sst);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "atmocor = %3d", input->atmocor);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "seawater_opt = %3d", input->seawater_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "aermodfile = %s", input->aermodfile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "aer_opt = %3d", input->aer_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "aer_wave_short = %3d", input->aer_wave_short);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "aer_wave_long = %3d", input->aer_wave_long);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "aer_swir_short = %3d", input->aer_swir_short);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "aer_swir_long = %3d", input->aer_swir_long);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "aer_rrs_short = %8.5f", input->aer_rrs_short);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "aer_rrs_long = %8.5f", input->aer_rrs_long);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "aer_angstrom = %8.5f", input->aer_angstrom);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "aer_iter_max = %3d", input->aer_iter_max);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "brdf_opt = %3d", input->brdf_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "gas_opt = %3d", input->gas_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "atrem_opt = %3d", input->atrem_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "atrem_full = %3d", input->atrem_full);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "atrem_geom = %3d", input->atrem_geom);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "atrem_model = %3d", input->atrem_model);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "iop_opt = %3d", input->iop_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "gsm_opt = %3d", input->gsm_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "gsm_fit = %3d", input->gsm_fit);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "gsm_adg_s = %8.5f", input->gsm_adg_s);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "gsm_bbp_s = %8.5f", input->gsm_bbp_s);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "gsm_aphw = %8.5f", input->gsm_aphw[0]);
    strcat(input->input_parms, str_buf);
    for (i = 1; i < numBands; i++) {
        sprintf(str_buf, ", %8.5f", input->gsm_aphw[i]);
        strcat(input->input_parms, str_buf);
    }
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "gsm_aphs = %8.5f", input->gsm_aphs[0]);
    strcat(input->input_parms, str_buf);
    for (i = 1; i < numBands; i++) {
        sprintf(str_buf, ", %8.5f", input->gsm_aphs[i]);
        strcat(input->input_parms, str_buf);
    }
    strcat(input->input_parms, "\n");
    sprintf(str_buf, "qaa_adg_s = %8.5f", input->qaa_adg_s);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "qaa_wave = %4d", input->qaa_wave[0]);
    strcat(input->input_parms, str_buf);
    for (i = 1; i < 5; i++) {
        sprintf(str_buf, ", %4d", input->qaa_wave[i]);
        strcat(input->input_parms, str_buf);
    }
    strcat(input->input_parms, "\n");
    sprintf(str_buf, "giop_maxiter = %3d", input->giop_maxiter);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "giop_fit_opt = %3d", input->giop_fit_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "giop_aph_opt = %3d", input->giop_aph_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "giop_acdom_opt = %3d", input->giop_acdom_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");
    
    sprintf(str_buf, "giop_anap_opt = %3d", input->giop_anap_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "giop_adg_opt = %3d", input->giop_adg_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "giop_bbp_opt = %3d", input->giop_bbp_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "giop_bbnap_opt = %3d", input->giop_bbnap_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");
    
    sprintf(str_buf, "giop_bbph_opt = %3d", input->giop_bbph_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "giop_rrs_opt = %3d", input->giop_rrs_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "giop_rrs_diff = %8.5f", input->giop_rrs_diff);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "giop_aph_file = %s", input->giop_aph_file);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "giop_aph_s = %8.5f", input->giop_aph_s);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "giop_adg_file = %s", input->giop_adg_file);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "giop_adg_s = %8.5f", input->giop_adg_s);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "giop_bbp_file = %s", input->giop_bbp_file);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "giop_bbp_s = %8.5f", input->giop_bbp_s);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "giop_acdom_file = %s", input->giop_acdom_file);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");
    
    sprintf(str_buf, "giop_anap_file = %s", input->giop_anap_file);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");
    
    sprintf(str_buf, "giop_bbph_file = %s", input->giop_bbph_file);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");
    
    sprintf(str_buf, "giop_bbnap_file = %s", input->giop_bbnap_file);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");
    
    sprintf(str_buf, "giop_grd = %8.5f", input->giop_grd[0]);
    strcat(input->input_parms, str_buf);
    for (i = 1; i < 2; i++) {
        sprintf(str_buf, ", %8.5f", input->giop_grd[i]);
        strcat(input->input_parms, str_buf);
    }
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "giop_wave = %6.1f", input->giop_wave[0]);
    strcat(input->input_parms, str_buf);
    for (i = 1; i < input->nbands; i++) {
        sprintf(str_buf, ", %6.1f", input->giop_wave[i]);
        strcat(input->input_parms, str_buf);
    }
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "giop_rrs_unc = %6.1f", input->giop_rrs_unc[0]);
    strcat(input->input_parms, str_buf);
    for (i = 1; i < input->nbands; i++) {
        sprintf(str_buf, ", %6.1f", input->giop_rrs_unc[i]);
        strcat(input->input_parms, str_buf);
    }
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "polfile = %s", input->polfile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "pol_opt = %3d", input->pol_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "rad_opt = %3d", input->rad_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "ocrvc_opt = %3d", input->ocrvc_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "vcnnfile = %s", input->vcnnfile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "absaer_opt = %3d", input->absaer_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "sl_pixl = %3d", input->sl_pixl);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "sl_frac = %8.4f", input->sl_frac);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "glint_opt = %3d", input->glint_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "resolution = %3d", input->resolution);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "outband_opt = %3d", input->outband_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "cirrus_opt = %3d", input->cirrus_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "oxaband_opt = %3d", input->oxaband_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "filter_opt = %3d", input->filter_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");
    if (input->filter_opt == 1) {
        sprintf(str_buf, "filter_file = %s", input->filter_file);
        strcat(input->input_parms, str_buf);
        strcat(input->input_parms, "\n");
        for (i = 0; i < input->fctl.nfilt; i++) {
            sprintf(str_buf, "# filter_%d = %d x %d (%d) %s %d ",
                    i + 1,
                    input->fctl.f[i].nx,
                    input->fctl.f[i].ny,
                    input->fctl.f[i].minfill,
                    filter_names[input->fctl.f[i].func],
                    input->fctl.f[i].band + 1);
            strcat(input->input_parms, str_buf);
            strcat(input->input_parms, "\n");
        }
    }

    sprintf(str_buf, "aerfile = %s", input->aerfile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "tgtfile = %s", input->tgtfile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "met1 = %s", input->met1);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "met2 = %s", input->met2);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "met3 = %s", input->met3);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "ozone1 = %s", input->ozone1);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "ozone2 = %s", input->ozone2);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "ozone3 = %s", input->ozone3);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "anc_cor_file = %s", input->anc_cor_file);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "land = %s", input->land);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "water = %s", input->water);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "demfile = %s", input->demfile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "elevfile = %s", input->elevfile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "elev_auxfile = %s", input->elev_auxfile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "icefile = %s", input->icefile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "sstfile = %s", input->sstfile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "sstreftype = %d", input->sstreftype);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "sssfile = %s", input->sssfile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "no2file = %s", input->no2file);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "alphafile = %s", input->alphafile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "tauafile = %s", input->tauafile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "picfile = %s", input->picfile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "owmcfile = %s", input->owmcfile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "prodxmlfile = %s", input->prodXMLfile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "breflectfile = %s", input->breflectfile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "aerbinfile = %s", input->aerbinfile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "gain = %8.4f", input->gain[0]);
    strcat(input->input_parms, str_buf);
    for (i = 1; i < numBands; i++) {
        sprintf(str_buf, ", %8.4f", input->gain[i]);
        strcat(input->input_parms, str_buf);
    }
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "gain_unc = %8.4f", input->gain_unc[0]);
    strcat(input->input_parms, str_buf);
    for (i = 1; i < numBands; i++) {
        sprintf(str_buf, ", %8.6f", input->gain_unc[i]);
        strcat(input->input_parms, str_buf);
    }
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "offset = %8.5f", input->offset[0]);
    strcat(input->input_parms, str_buf);
    for (i = 1; i < numBands; i++) {
        sprintf(str_buf, ", %8.5f", input->offset[i]);
        strcat(input->input_parms, str_buf);
    }
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "naermodels = %d", input->naermodels);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "aermodels = %3s", input->aermodels[0]);
    strcat(input->input_parms, str_buf);
    for (i = 1; i < input->naermodels; i++) {
        sprintf(str_buf, ", %3s", input->aermodels[i]);
        strcat(input->input_parms, str_buf);
    }
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "taua = %8.4f", input->taua[0]);
    strcat(input->input_parms, str_buf);
    for (i = 1; i < numBands; i++) {
        sprintf(str_buf, ", %8.4f", input->taua[i]);
        strcat(input->input_parms, str_buf);
    }
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "aermodrat = %8.5f", input->aermodrat);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "aermodmin = %3d", input->aermodmin);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "aermodmax = %3d", input->aermodmax);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "cloud_wave = %8.3f", input->cloud_wave);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "cloud_thresh = %8.3f", input->albedo);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "cloud_eps = %8.3f", input->cloud_eps);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "cirrus_thresh = %8.5f, %8.5f", input->cirrus_thresh[0], input->cirrus_thresh[1]);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "absaer = %8.3f", input->absaer);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "rhoamin = %8.5f", input->rhoamin);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "nlwmin = %8.3f", input->nlwmin);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "hipol = %8.3f", input->hipol);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "wsmax = %8.3f", input->wsmax);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "coccolith = %8.4f", input->coccolith[0]);
    strcat(input->input_parms, str_buf);
    for (i = 1; i < 8; i++) {
        sprintf(str_buf, ", %8.4f", input->coccolith[i]);
        strcat(input->input_parms, str_buf);
    }
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "tauamax = %8.3f", input->tauamax);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "epsmin = %8.3f", input->epsmin);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "epsmax = %8.3f", input->epsmax);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "glint_thresh = %8.3f", input->glint);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "windspeed = %8.3f", input->windspeed);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "windangle = %8.3f", input->windangle);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "pressure = %8.3f", input->pressure);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "ozone = %8.3f", input->ozone);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "watervapor = %8.3f", input->watervapor);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "relhumid = %8.3f", input->relhumid);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "ice_threshold = %8.3f", input->ice_threshold);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "sunzen = %8.3f", input->sunzen);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "satzen = %8.3f", input->satzen);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "maskland = %2d", input->landmask);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "maskbath = %2d", input->bathmask);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "maskcloud = %2d", input->cloudmask);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "maskglint = %2d", input->glintmask);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "masksunzen = %2d", input->sunzenmask);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "masksatzen = %2d", input->satzenmask);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "maskhilt = %2d", input->hiltmask);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "maskstlight = %2d", input->stlightmask);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "mumm_alpha = %8.3f", input->mumm_alpha);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "mumm_gamma = %8.3f", input->mumm_gamma);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "mumm_epsilon = %8.3f", input->mumm_epsilon);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "chloc2_wave = %4d", input->chloc2w[0]);
    strcat(input->input_parms, str_buf);
    for (i = 1; i < 2; i++) {
        sprintf(str_buf, ", %4d", input->chloc2w[i]);
        strcat(input->input_parms, str_buf);
    }
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "chloc2_coef = %8.5f", input->chloc2c[0]);
    strcat(input->input_parms, str_buf);
    for (i = 1; i < 5; i++) {
        sprintf(str_buf, ", %8.5f", input->chloc2c[i]);
        strcat(input->input_parms, str_buf);
    }
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "chloc3_wave = %4d", input->chloc3w[0]);
    strcat(input->input_parms, str_buf);
    for (i = 1; i < 3; i++) {
        sprintf(str_buf, ", %4d", input->chloc3w[i]);
        strcat(input->input_parms, str_buf);
    }
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "chloc3_coef = %8.5f", input->chloc3c[0]);
    strcat(input->input_parms, str_buf);
    for (i = 1; i < 5; i++) {
        sprintf(str_buf, ", %8.5f", input->chloc3c[i]);
        strcat(input->input_parms, str_buf);
    }
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "chloc4_wave = %4d", input->chloc4w[0]);
    strcat(input->input_parms, str_buf);
    for (i = 1; i < 4; i++) {
        sprintf(str_buf, ", %4d", input->chloc4w[i]);
        strcat(input->input_parms, str_buf);
    }
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "chloc4_coef = %8.5f", input->chloc4c[0]);
    strcat(input->input_parms, str_buf);
    for (i = 1; i < 5; i++) {
        sprintf(str_buf, ", %8.5f", input->chloc4c[i]);
        strcat(input->input_parms, str_buf);
    }
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "kd2_wave = %4d", input->kd2w[0]);
    strcat(input->input_parms, str_buf);
    for (i = 1; i < 2; i++) {
        sprintf(str_buf, ", %4d", input->kd2w[i]);
        strcat(input->input_parms, str_buf);
    }
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "kd2_coef = %8.5f", input->kd2c[0]);
    strcat(input->input_parms, str_buf);
    for (i = 1; i < 6; i++) {
        sprintf(str_buf, ", %8.5f", input->kd2c[i]);
        strcat(input->input_parms, str_buf);
    }
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "flh_offset = %8.5f", input->flh_offset);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "btfile = %s", input->btfile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "sstcoeffile = %s", input->sstcoeffile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "sstssesfile = %s", input->sstssesfile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "sst4coeffile = %s", input->sst4coeffile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "sst4ssesfile = %s", input->sst4ssesfile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "sst3coeffile = %s", input->sst3coeffile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "sst3ssesfile = %s", input->sst3ssesfile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "vcal_opt = %3d", input->vcal_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "vcal_nlw = %8.4f", input->vcal_nLw[0]);
    strcat(input->input_parms, str_buf);
    for (i = 1; i < numBands; i++) {
        sprintf(str_buf, ", %8.4f", input->vcal_nLw[i]);
        strcat(input->input_parms, str_buf);
    }
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "vcal_lw = %8.4f", input->vcal_Lw[0]);
    strcat(input->input_parms, str_buf);
    for (i = 1; i < numBands; i++) {
        sprintf(str_buf, ", %8.4f", input->vcal_Lw[i]);
        strcat(input->input_parms, str_buf);
    }
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "vcal_chl = %8.4f", input->vcal_chl);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "vcal_solz = %8.4f", input->vcal_solz);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "vcal_depth = %8.4f", input->vcal_depth);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "vcal_min_nbin = %d", input->vcal_min_nbin);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "vcal_min_nscene = %d", input->vcal_min_nscene);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "xcalfile = %s", input->xcalfile);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "xcal_opt = %3d", input->xcal_opt[0]);
    strcat(input->input_parms, str_buf);
    for (i = 1; i < input->xcal_nwave; i++) {
        sprintf(str_buf, ", %3d", input->xcal_opt[i]);
        strcat(input->input_parms, str_buf);
    }

    strcat(input->input_parms, "\n");
    sprintf(str_buf, "xcal_wave = %8.4f", input->xcal_wave[0]);
    strcat(input->input_parms, str_buf);
    for (i = 1; i < input->xcal_nwave; i++) {
        sprintf(str_buf, ", %8.4f", input->xcal_wave[i]);
        strcat(input->input_parms, str_buf);
    }
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "band_shift_opt = %3d", input->band_shift_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "add_noise_sigma = %2.4f",input->add_noise_sigma);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "noise_scale = %2.4f",input->noise_scale);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");
    
    sprintf(str_buf, "wiggle_band = %3d",input->wiggle_band);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "wiggle_by = %2.4f",input->wiggle_by);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "stype = %d", input->stype);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "datamin = %8.4f", input->datamin);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "datamax = %8.4f", input->datamax);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "north = %8.4f", input->north);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "south = %8.4f", input->south);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "east = %8.4f", input->east);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "west = %8.4f", input->west);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "xbox = %d", input->xbox);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "ybox = %d", input->ybox);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");
    
    //new Raman test
    sprintf(str_buf, "raman_opt = %d", input->raman_opt);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "width = %d", input->width);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "threshold = %8.4f", input->threshold);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "rgb = %d", input->rgb[0]);
    strcat(input->input_parms, str_buf);
    for (i = 1; i < 3; i++) {
        sprintf(str_buf, ", %d", input->rgb[i]);
        strcat(input->input_parms, str_buf);
    }
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "subsamp = %d", input->subsamp);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "viirsnv7 = %d", input->viirsnv7);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "viirsnosisaf = %d", input->viirsnosisaf);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "sstrefdif = %8.4f", input->sstrefdif);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "newavhrrcal = %d", input->newavhrrcal);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "ch22detcor = %9.6f", input->ch22detcor[0]);
    strcat(input->input_parms, str_buf);
    for (i = 1; i < 10; i++) {
        sprintf(str_buf, ", %9.6f", input->ch22detcor[i]);
        strcat(input->input_parms, str_buf);
    }
    strcat(input->input_parms, "\n");

    sprintf(str_buf, "ch23detcor = %9.6f", input->ch23detcor[0]);
    strcat(input->input_parms, str_buf);
    for (i = 1; i < 10; i++) {
        sprintf(str_buf, ", %9.6f", input->ch23detcor[i]);
        strcat(input->input_parms, str_buf);
    }
    strcat(input->input_parms, "\n");

    /*                                                                  */
    /* Build string of input files for meta-data documentation          */
    /*                                                                  */
    for (i = 0; i < MAX_IFILES; i++) {
        if (input->ifile[i][0] != '\0') {
            tmp_str = strrchr(input->ifile[i], '/');
            tmp_str = (tmp_str == 0x0) ? input->ifile[i] : tmp_str + 1;
            if (i == 0) sprintf(input->input_files, "%s", tmp_str);
            else {
                sprintf(str_buf, ",%s", tmp_str);
                strcat(input->input_files, str_buf);
            }
        } else break;
    }
    if (strlen(input->geofile)){
        tmp_str = strrchr(input->geofile, '/');
        tmp_str = (tmp_str == 0x0) ? input->geofile : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    if (strlen(input->viirscalparfile)){
        tmp_str = strrchr(input->viirscalparfile, '/');
        tmp_str = (tmp_str == 0x0) ? input->viirscalparfile : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    if (strlen(input->met1)){
        tmp_str = strrchr(input->met1, '/');
        tmp_str = (tmp_str == 0x0) ? input->met1 : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    if (strlen(input->met2)){
        tmp_str = strrchr(input->met2, '/');
        tmp_str = (tmp_str == 0x0) ? input->met2 : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    if (strlen(input->met3)){
        tmp_str = strrchr(input->met3, '/');
        tmp_str = (tmp_str == 0x0) ? input->met3 : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    if (strlen(input->ozone1)){
        tmp_str = strrchr(input->ozone1, '/');
        tmp_str = (tmp_str == 0x0) ? input->ozone1 : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    if (strlen(input->ozone2)){
        tmp_str = strrchr(input->ozone2, '/');
        tmp_str = (tmp_str == 0x0) ? input->ozone2 : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    if (strlen(input->ozone3)){
        tmp_str = strrchr(input->ozone3, '/');
        tmp_str = (tmp_str == 0x0) ? input->ozone3 : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    if (strlen(input->anc_cor_file)){
        tmp_str = strrchr(input->anc_cor_file, '/');
        tmp_str = (tmp_str == 0x0) ? input->anc_cor_file : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    if (strlen(input->calfile)){
        tmp_str = strrchr(input->calfile, '/');
        tmp_str = (tmp_str == 0x0) ? input->calfile : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    if (strlen(input->xcalfile)){
        tmp_str = strrchr(input->xcalfile, '/');
        tmp_str = (tmp_str == 0x0) ? input->xcalfile : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    if (strlen(input->fqfile)){
        tmp_str = strrchr(input->fqfile, '/');
        tmp_str = (tmp_str == 0x0) ? input->fqfile : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    if (strlen(input->parfile)){
        tmp_str = strrchr(input->parfile, '/');
        tmp_str = (tmp_str == 0x0) ? input->parfile : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    if (strlen(input->polfile)){
        tmp_str = strrchr(input->polfile, '/');
        tmp_str = (tmp_str == 0x0) ? input->polfile : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    if (strlen(input->aermodfile)){
        tmp_str = strrchr(input->aermodfile, '/');
        tmp_str = (tmp_str == 0x0) ? input->aermodfile : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    if (strlen(input->vcnnfile)){
        tmp_str = strrchr(input->vcnnfile, '/');
        tmp_str = (tmp_str == 0x0) ? input->vcnnfile : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    if (strlen(input->cldfile)){
        tmp_str = strrchr(input->cldfile, '/');
        tmp_str = (tmp_str == 0x0) ? input->cldfile : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    if (strlen(input->land)){
        tmp_str = strrchr(input->land, '/');
        tmp_str = (tmp_str == 0x0) ? input->land : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    if (strlen(input->water)){
        tmp_str = strrchr(input->water, '/');
        tmp_str = (tmp_str == 0x0) ? input->water : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    if (strlen(input->demfile)){
        tmp_str = strrchr(input->demfile, '/');
        tmp_str = (tmp_str == 0x0) ? input->demfile : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    if (strlen(input->elevfile)){
        tmp_str = strrchr(input->elevfile, '/');
        tmp_str = (tmp_str == 0x0) ? input->elevfile : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    if (strlen(input->elev_auxfile)){
        tmp_str = strrchr(input->elev_auxfile, '/');
        tmp_str = (tmp_str == 0x0) ? input->elev_auxfile : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    if (strlen(input->icefile)){
        tmp_str = strrchr(input->icefile, '/');
        tmp_str = (tmp_str == 0x0) ? input->icefile : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    if (strlen(input->sstfile)){
        tmp_str = strrchr(input->sstfile, '/');
        tmp_str = (tmp_str == 0x0) ? input->sstfile : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }

    sprintf(str_buf, "sstreftype = %d", input->sstreftype);
    strcat(input->input_parms, str_buf);
    strcat(input->input_parms, "\n");

    if (strlen(input->sssfile)){
        tmp_str = strrchr(input->sssfile, '/');
        tmp_str = (tmp_str == 0x0) ? input->sssfile : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    if (strlen(input->no2file)){
        tmp_str = strrchr(input->no2file, '/');
        tmp_str = (tmp_str == 0x0) ? input->no2file : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    if (strlen(input->alphafile)){
        tmp_str = strrchr(input->alphafile, '/');
        tmp_str = (tmp_str == 0x0) ? input->alphafile : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    if (strlen(input->tauafile)){
        tmp_str = strrchr(input->tauafile, '/');
        tmp_str = (tmp_str == 0x0) ? input->tauafile : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }

    if (strlen(input->picfile)){
        tmp_str = strrchr(input->picfile, '/');
        tmp_str = (tmp_str == 0x0) ? input->picfile : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    if (strlen(input->owmcfile)){
        tmp_str = strrchr(input->owmcfile, '/');
        tmp_str = (tmp_str == 0x0) ? input->owmcfile : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input->input_files, str_buf);
    }
    /*                                                                  */
    /* Build string of mask names for meta-data documentation           */
    /*                                                                  */
    strcpy(input->mask_names, "ATMFAIL");
    if (input->landmask == 1) strcat(input->mask_names, ",LAND");
    if (input->bathmask == 1) strcat(input->mask_names, ",COASTZ");
    if (input->cloudmask == 1) strcat(input->mask_names, ",CLDICE");
    if (input->glintmask == 1) strcat(input->mask_names, ",HIGLINT");
    if (input->sunzenmask == 1) strcat(input->mask_names, ",HISOLZEN");
    if (input->satzenmask == 1) strcat(input->mask_names, ",HISATZEN");
    if (input->hiltmask == 1) strcat(input->mask_names, ",HILT");
    if (input->stlightmask == 1) strcat(input->mask_names, ",STRAYLIGHT");

    return 0;
}

/*-----------------------------------------------------------------------------
    Function:  msl12_input

    Returns:   int (status)
        The return code is a negative value if any error occurs, otherwise,
        returns 0.

    Description:
        Convert the arguments from the command line into a structure input
        variable.

    Parameters: (in calling order)
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int       argc          I   number of arguments
        char      **argv        I   list of arguments
        instr     input     O   structure variable for inputs

----------------------------------------------------------------------------*/
int msl12_input(int argc, char *argv[], const char* progName, instr *input, filehandle *l1file) {
    /* hold all of the command line options */
    clo_optionList_t* list;

    list = clo_createList();
    clo_setEnableExtraOptions(1);

    /* initialize the option list with descriptions and default values */
    l2gen_init_options(list, progName);

    return msl12_option_input(argc, argv, list, (char*) progName, input, l1file);
}


/* This function takes a predefined L1 file handle and loads all defaults    */

/* It's used by MSl1info, MSll2snpx, etc which don't use standard par files. */
int msl12_input_defaults(filehandle *l1file, instr *input) {
    int argc = 0;
    char *argv[10];
    static char resolutionStr[FILENAME_MAX];
    static char ifileStr[FILENAME_MAX];
    static char geofileStr[FILENAME_MAX];

    argv[argc++] = "msl12";

    sprintf(ifileStr, "ifile=%s", l1file->name);
    argv[argc++] = ifileStr;

    if (l1file->geofile) {
        sprintf(geofileStr, "geofile=%s", l1file->geofile);
        argv[argc++] = geofileStr;
        argc = 3;
    }

    if (input->resolution != -1) {
        sprintf(resolutionStr, "resolution=%d", input->resolution);
        argv[argc++] = resolutionStr;
    }

    if (msl12_input(argc, argv, "msl12", input, l1file))
        return (-1);

    return 0;
}

int l2gen_usage(const char *prog) {
    int i;
    clo_optionList_t* list;
    clo_option_t* option;

    list = clo_createList();
    l2gen_init_options(list, prog);
    clo_printUsage(list);

    return 0;
}
