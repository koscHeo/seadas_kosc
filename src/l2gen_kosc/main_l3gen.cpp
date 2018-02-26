// =====================================================================
// l3gen - level-3 to level-3 ocean algorithm processor                 
// B. Franz, NASA/OBPG, August 2008            
//                         
// Modification History
// --------------------
// Use hdf_bin class rather than hdf4_bin class
// J. Gales     Futuretech     07/18/11
// =====================================================================

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <list>
#include "netcdf.h"
#include "hdf_bin.h"
#include <seaproto.h>
#include "version.h"

#define NBINREAD MAXPIX
#define FLAGMASK PRODFAIL

extern "C" {
#include "l12_proto.h"
}

using namespace std;

typedef Hdf::binListStruct blstr;

void l3gen_usage(char *prog) {
    l2gen_usage("l3gen");
}

void load_l3file_handle(instr *input, int filenum, filehandle *file) {
    strcpy(file->name, (const char *) input->ofile[filenum]);
    strcpy(file->l2prod, (const char *) input->l2prod[filenum]);
    strcpy(file->def_l2prod, input->def_l2prod[filenum]);

    file->format = FMT_L3BIN;
    file->mode = WRITE;
    file->sensorID = input->sensorID;
    file->pro_control = input->pro_control;
    file->ocrvc_opt = input->ocrvc_opt;

    file->tot_prod = prodlist(input->sensorID, input->evalmask, file->l2prod,
            file->def_l2prod, file->l2_prod_names);
}

void load_l12(int32 nscans, int32 npix, instr *input, l1str *l1rec,
        l2str *l2rec) {
    int32 sensorID = input->sensorID;
    int32 evalmask = input->evalmask;
    int32 nbands;
    int32 nbandsir;
    int32 *Lambda;
    float *Fobar;
    float *Tau_r;
    float *k_oz;
    float *aw;
    float *bbw;
    int32 iw;

    nbands = rdsensorinfo(sensorID, evalmask, NULL, NULL);
    nbandsir = rdsensorinfo(sensorID, evalmask, "NbandsIR", NULL);

    rdsensorinfo(sensorID, evalmask, "Lambda", (void **) &Lambda);
    rdsensorinfo(sensorID, evalmask, "Fobar", (void **) &Fobar);
    rdsensorinfo(sensorID, evalmask, "Tau_r", (void **) &Tau_r);
    rdsensorinfo(sensorID, evalmask, "k_oz", (void **) &k_oz);
    rdsensorinfo(sensorID, evalmask, "aw", (void **) &aw);
    rdsensorinfo(sensorID, evalmask, "bbw", (void **) &bbw);

    memset(l1rec->data, 0, l1rec->length);
    memset(l2rec->data, 0, l2rec->length);

    for (iw = 0; iw < nbands; iw++) {

        l1rec->iwave[iw] = Lambda[iw];
        l1rec->fwave[iw] = (float) Lambda[iw];
        l1rec->Tau_r[iw] = Tau_r[iw];
        l1rec->k_oz[iw] = k_oz[iw];
        l1rec->aw[iw] = aw[iw];
        l1rec->bbw[iw] = bbw[iw];
        l1rec->Fobar[iw] = Fobar[iw];
        l1rec->Fo[iw] = Fobar[iw];  // for the moment

        get_f0_thuillier_ext(l1rec->iwave[iw], BANDW, &l1rec->Fonom[iw]);
    }

    l2rec->sensorID = l1rec->sensorID = sensorID;
    l2rec->nbands = l1rec->nbands = nbands;
    l2rec->nbandsir = l1rec->nbandsir = nbandsir;
    l2rec->tilt = l1rec->tilt = 0.0;
    l2rec->tilt = l1rec->tilt = 0.0;
    l2rec->mside = l1rec->mside = 0;
    l2rec->detnum = l1rec->detnum = 0;
    l2rec->ndets = l1rec->ndets = 1;
    l2rec->nscans = l1rec->nscans = nscans;
    l2rec->input = l1rec->input = input;

    l1rec->n_inprods = 0;
    l1rec->spix = 0;
    l1rec->epix = npix - 1;
    l1rec->dpix = 1;
    l1rec->sscan = 0;
    l1rec->escan = nscans - 1;
    l1rec->dscan = 1;

    int32 ip;
    for (ip = 0; ip < npix; ip++) {
        l1rec->pixnum[ip] = 0;
        l1rec->alpha[ip] = 0.0;
        l1rec->flags[ip] = 0;
        l1rec->mask[ip] = 0;
    }

    // Set pointers to shared data areas

    l2rec->ws = l1rec->ws;
    l2rec->wd = l1rec->wd;
    l2rec->mw = l1rec->mw;
    l2rec->zw = l1rec->zw;
    l2rec->pr = l1rec->pr;
    l2rec->oz = l1rec->oz;
    l2rec->wv = l1rec->wv;
    l2rec->rh = l1rec->rh;
    l2rec->pixnum = l1rec->pixnum;
    l2rec->alpha = l1rec->alpha;
    l2rec->no2_tropo = l1rec->no2_tropo;
    l2rec->no2_strat = l1rec->no2_strat;
    l2rec->no2_frac = l1rec->no2_frac;
    l2rec->height = l1rec->height;
    l2rec->sstref = l1rec->sstref;
    l2rec->iwave = l1rec->iwave;
    l2rec->fwave = l1rec->fwave;
    l2rec->Fo = l1rec->Fo;
    l2rec->Fobar = l1rec->Fobar;
    l2rec->Fonom = l1rec->Fonom;
    l2rec->Tau_r = l1rec->Tau_r;
    l2rec->k_oz = l1rec->k_oz;
    l2rec->aw = l1rec->aw;
    l2rec->bbw = l1rec->bbw;
    l2rec->delphi = l1rec->delphi;
    l2rec->alpha = l1rec->alpha;
    l2rec->flags = l1rec->flags;
    l2rec->mask = l1rec->mask;
    l2rec->rhos = l1rec->rhos;
    l2rec->sw_n = l1rec->sw_n;
    l2rec->sw_a = l1rec->sw_a;
    l2rec->sw_bb = l1rec->sw_bb;

    bindex_set((int32_t*) Lambda, nbands + nbandsir, BANDW);
}

// -------------------------------------------------------------------- 
//                            main                                      
// -------------------------------------------------------------------- 
int main(int argc, char* argv[]) {
    static l1str *l1rec;               // generic level-1b scan structure
    static l2str *l2rec;               // generic level-2  scan structure
    static instr *input;               // input parameters structure
    static filehandle ifile;            // input file handle 
    static filehandle ofile;            // output file handle 
    int iprod;

    char *ptime = ydhmsf(now(), 'G');
    double start_time = now();

    int16 year, day;
    int32_t syear;
    int32_t sday;
    int msec;
    double dsec;
    double mtime;

    char soft_id[200];
    float gmt;
    float lon;
    float lat;
    float solz;
    float sola;
    int status;
    char buf[FILENAME_MAX];

    int i;
    int ip;
    int offset;
    int bin_num;
    int nobs;
    int nscenes;
    int nwrite;
    static bool atLeastOne = false;

    if (argc == 1) {
        l3gen_usage(argv[0]);
        return 0;
    }

    // see if help on command line
    for (i = 0; i < argc; i++) {
        if ((strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "-help") == 0)) {
            l3gen_usage(argv[0]);
            return 1;
        }
    }

    setvbuf(stdout, NULL, _IOLBF, 0);
    setvbuf(stderr, NULL, _IOLBF, 0);

    // Initialize file handles
    input = (instr*) malloc(sizeof(instr));

    cdata_();
    filehandle_init(&ifile);
    ifile.input = input;
    filehandle_init(&ofile);
    ofile.input = input;

    // Parse input parameters
    if (msl12_input(argc, argv, "l3gen", input, &ifile) != 0) {
        printf("-E- %s: Error parsing input parameters.\n", argv[0]);
        exit(1);
    }
    if (input->ocrvc_opt)
        input->sensorID = OCRVC;

    char     proc_con[2048];
    strcpy(proc_con, basename(argv[0]));
    for (i=1; i < argc; i++) {
      strcat(proc_con, " ");
      strcat(proc_con, argv[i]);
    }

    // Transfer input info to file handles
    load_l3file_handle(input, 0, &ofile);

    // Build output product list

    cout << endl << "The following products will be included in " << ofile.name
            << endl;
    int32 nprods_out = ofile.tot_prod;
    char prodlist_out[PRODSTRLEN] = "";
    for (iprod = 0; iprod < ofile.tot_prod; iprod++) {
        cout << iprod + 1 << " " << ofile.l2_prod_names[iprod] << endl;
        strcat(prodlist_out, ofile.l2_prod_names[iprod]);
        if (iprod < ofile.tot_prod - 1)
            strcat(prodlist_out, ",");
    }

    // Open input binfile and get full input product list

    static Hdf::hdf_bin *input_binfile;
    const char *inputfile;

    inputfile = input->ifile[0];
    input_binfile = Hdf::openBinObject(inputfile);

    int len = input_binfile->query();
    char *fullprodlist = (char *) malloc(len);
    input_binfile->query(fullprodlist);
    int nrows = input_binfile->nrows;

    // Get input bin dimension

    int32 nbins = input_binfile->n_data_records;
    int32 npix = NBINREAD;

    // Get mean time
    mtime = (input_binfile->meta_l3b.startTime + input_binfile->meta_l3b.endTime)/2.0;

    unix2yds(mtime, &year, &day, &dsec);
    msec = long(dsec * 1000.0);
    //yes, silly, but easier than fixing the various routines that need either a 16 or 32bit int for year/day
    syear = year;
    sday = day;
    // Allocate memory for L1 and L2 scan data (the L2 record shares
    // space with the L1 record, hence we need to allocate both)
    l1rec = (l1str*) malloc(sizeof(l1str));
    l2rec = (l2str*) malloc(sizeof(l2str));

    if (alloc_l1(npix, input->nbands, NBANDSIR, 0, l1rec) == 0) {
        printf("-E- %s: Unable to allocate L1 record.\n", argv[0]);
        exit(1);
    }
    if (alloc_l2(npix, input->nbands, l2rec) == 0) {
        printf("-E- %s: Unable to allocate L2 record.\n", argv[0]);
        exit(1);
    }

    // Add meta-data to L1 and L2 structures - use a single pixel just to get things
    // started - load_l12 will be called for each bin row with an appropriate # of pixels
    load_l12(1, npix, input, l1rec, l2rec);

    // Build input product request list

    int32 n_vvv = MIN(windex(700., l1rec->fwave, l1rec->nbands) + 1,
            l1rec->nbands);
    char prodstr[PRODSTRLEN] = "";

    //   what we have

    char **prodlist_in;
    int32 nprods_in = input_binfile->query(&prodlist_in);

    //   what we want

    char prodlist_wt[MAXPROD][32];
    int32 nprods_wt = prodlist(input->sensorID, input->evalmask,
            "Rrs_vvv nLw_vvv Rrs_vc_vvv solz sola senz sena relaz", "",
            prodlist_wt);
    enum geo_order {
        SOLZ, SOLA, SENZ, SENA, RELAZ
    };

    //   what we have that we want

    int32 *rrs_bindex;
    int32 *rrs_windex;
    int32 n_rrs = 0;
    int32 *nlw_bindex;
    int32 *nlw_windex;
    int32 n_nlw = 0;
    int32 geo_bindex[5];
    int32 geo_windex[5];
    int32 n_geo = 0;
    int32 have_solz = 0;
    int32 have_sola = 0;
    int32 have_sena = 0;
    int32 have_relaz = 0;

    if ((rrs_bindex = (int32 *) malloc(input->nbands * sizeof(int32))) == NULL) {
        printf("-E- %s: Error allocating memory to the input data array.\n",
                argv[0]);
        exit(FATAL_ERROR);
    }
    if ((rrs_windex = (int32 *) malloc(input->nbands * sizeof(int32))) == NULL) {
        printf("-E- %s: Error allocating memory to the input data array.\n",
                argv[0]);
        exit(FATAL_ERROR);
    }
    if ((nlw_bindex = (int32 *) malloc(input->nbands * sizeof(int32))) == NULL) {
        printf("-E- %s: Error allocating memory to the input data array.\n",
                argv[0]);
        exit(FATAL_ERROR);
    }
    if ((nlw_windex = (int32 *) malloc(input->nbands * sizeof(int32))) == NULL) {
        printf("-E- %s: Error allocating memory to the input data array.\n",
                argv[0]);
        exit(FATAL_ERROR);
    }

    cout << endl << "Checking input for desired products: " << endl;
    iprod = 0;
    for (int32 iwant = 0; iwant < nprods_wt; iwant++) {
        for (int32 ihave = 0; ihave < nprods_in; ihave++) {
            if (strcmp(prodlist_in[ihave], prodlist_wt[iwant]) == 0) {
                strcat(prodstr, prodlist_in[ihave]);
                strcat(prodstr, ",");
                cout << "found " << prodlist_wt[iwant] << endl;
                if (iwant < n_vvv) {
                    rrs_bindex[n_rrs] = iprod;        // posiiton in bin read
                    rrs_windex[n_rrs] = iwant; // position in sensor wavelengths
                    n_rrs++;
                    iprod++;
                } else if (iwant < n_vvv * 2) {
                    nlw_bindex[n_nlw] = iprod;        // positon in bin read
                    nlw_windex[n_nlw] = iwant - n_vvv; // position in sensor wavelengths
                    n_nlw++;
                    iprod++;
                } else if (iwant < n_vvv * 3) { // Rrs_vc products overwriting Rrs
                    rrs_bindex[n_rrs] = iprod;        // posiiton in bin read
                    rrs_windex[n_rrs] = iwant - 2 * n_vvv; // position in sensor wavelengths
                    n_rrs++;
                    iprod++;
                } else {
                    geo_bindex[n_geo] = iprod;
                    geo_windex[n_geo] = iwant - 3 * n_vvv;
                    switch (geo_windex[n_geo]) {
                    case SOLZ:
                        have_solz = 1;
                        break;
                    case SOLA:
                        have_sola = 1;
                        break;
                    case SENA:
                        have_sena = 1;
                        break;
                    case RELAZ:
                        have_relaz = 1;
                        break;
                    }
                    n_geo++;
                    iprod++;
                }
                break;
            }
        }
    }
    nprods_in = iprod;

    // remove trailing comma
    len = strlen(prodstr);
    prodstr[len - 1] = '\0';

    // Allocate data arrays and bin list and masking flag
    float **outData;
    float **tmpData;
    float **inData;
    if ((inData = (float **) malloc(nprods_in * sizeof(float *))) == NULL) {
        printf("-E- %s: Error allocating memory to the input data array.\n",
                argv[0]);
        exit(FATAL_ERROR);
    }
    if ((outData = (float **) malloc(nprods_out * sizeof(float *))) == NULL) {
        printf("-E- %s: Error allocating memory to the output data array.\n",
                argv[0]);
        exit(FATAL_ERROR);
    }
    if ((tmpData = (float **) malloc(nprods_out * sizeof(float *))) == NULL) {
        printf("-E- %s: Error allocating memory to the temp data array.\n",
                argv[0]);
        exit(FATAL_ERROR);
    }
    for (i = 0; i < nprods_in; i++) {
        if ((inData[i] = (float *) calloc(2*npix, sizeof(float))) == NULL) {
            printf("-E- %s: Error allocating memory to the output data array.\n",
                    argv[0]);
            exit(FATAL_ERROR);
        }
    }
    for (i = 0; i < nprods_out; i++) {
        if ((outData[i] = (float *) calloc(2*npix, sizeof(float))) == NULL) {
            printf(
                    "-E- %s: Error allocating memory to the output data array.\n",
                    argv[0]);
            exit(FATAL_ERROR);
        }
        if ((tmpData[i] = (float *) calloc(npix, sizeof(float))) == NULL) {
            printf(
                    "-E- %s: Error allocating memory to the output data array.\n",
                    argv[0]);
            exit(FATAL_ERROR);
        }
    }
    char *outMask = (char *) calloc(npix, sizeof(char));

    // set up the bin file reading for the products we need
    input_binfile->read(prodstr);

    /* Create output file */
    /* ------------------ */
    Hdf::hdf_bin *output_binfile;

    /*
     * If the input structure does not set the output format,
     * make the output format match the input L3 format
     */
    if (getFileFormatName(input->oformat) == NULL) {
        if (input_binfile->isHDF5) {
            strcpy(input->oformat, "HDF5");
        } else if (input_binfile->isCDF4) {
            strcpy(input->oformat, "netCDF4");
        } else {
            strcpy(input->oformat, "HDF4");
        }
    }

    if (strcmp(input->oformat, "HDF4") == 0) {
        output_binfile = new Hdf::hdf4_bin;
        output_binfile->hasNoext = true;
    }
    if (strcmp(input->oformat, "HDF5") == 0)
        output_binfile = new Hdf::hdf5_bin;
    if (strcmp(input->oformat, "netCDF4") == 0)
        output_binfile = new Hdf::cdf4_bin;

    output_binfile->create(ofile.name, input_binfile->nrows);

    if (input_binfile->isCDF4)
        output_binfile->deflate = input->deflate;

    strcpy(output_binfile->meta_l3b.product_name, ofile.name);

    strncpy(output_binfile->meta_l3b.ptime, ptime, 16);

    strcpy( output_binfile->meta_l3b.proc_con, proc_con);

    strcpy( output_binfile->meta_l3b.input_parms,input->input_parms);
    // Begin processing
    for (int32 iscan = 0; iscan < nrows; iscan++) {
        if (iscan % 50 == 0){
            int secs = (int) now()-start_time;
            cout << "Reading row " << iscan  << " of " << nrows << " after "<< secs <<" seconds"<< endl;
        }
        // Read the inputs
        // Get basebin and numbin for this input row
        int basebin = input_binfile->get_basebin( iscan);
        int numbin = input_binfile->get_numbin( iscan);
        input_binfile->readBinIndex( iscan);
        int ext = input_binfile->get_ext();
        int beg = input_binfile->get_beg();
        // if the row has no filled bins, skip it
        if (beg == 0)
            continue;

        int binList_ptr = input_binfile->get_list_ptr();
        input_binfile->readBinList(ext);

        // fill the input data array with the necessary input data
        i=0;
        for (int j = 0; j < input_binfile->nprod(); j++) {

            if (input_binfile->active_data_prod[j] == true) {
                input_binfile->get_prodname(j, buf);
                if (strcmp(prodlist_in[j], buf) == 0) {
                    input_binfile->readSums(inData[i], ext, j);
                    i++;
                }
            }
        }
        if ( input_binfile->isHDF5 || input_binfile->isCDF4)
            input_binfile->setDataPtr( ext);


        // Stuff the L2 record
        // only process as many pixels as the row as filled bins
        memset(outMask, '\0', npix);
        init_l1(l1rec, ext,input->nbands);
        init_l2(l2rec,input->nbands);



        l1rec->npix = ext;
        l1rec->epix = ext - 1;
        l2rec->npix = ext;
        ifile.npix = ext;
        ifile.epix = ext-1;
        ifile.terrain_corrected = 1;

        load_l12(1, ext, input, l1rec, l2rec);

        l2rec->iscan = l1rec->iscan = iscan;
        *(l1rec->year) = *(l2rec->year) = year;
        *(l1rec->day) = *(l2rec->day) = day;
        *(l1rec->msec) = *(l2rec->msec) = msec;

        // Fill in the l1/2 structures with data from input bin file
        // nLw, Rrs, geometries - if available.
        for (ip = 0; ip < ext; ip++) {
            float weight = input_binfile->get_weights(ip);
            int32 ibw, ipw;
            float lat;
            float lon;
            bin_num = input_binfile->get_bin_num(ip);
            bin2ll( bin_num, &lat, &lon);
            l1rec->lon[ip] = l2rec->lon[ip] = lon;
            l1rec->lat[ip] = l2rec->lat[ip] = lat;
            l1rec->nobs[ip] = input_binfile->get_nobs(ip);
            for (int32 iw = 0; iw < n_nlw; iw++) {
                ipw = ip * input->nbands + nlw_windex[iw];
                l2rec->nLw[ipw] = inData[nlw_bindex[iw]][2 * ip] / weight;
                l2rec->Rrs[ipw] = MAX(l2rec->nLw[ipw] / l2rec->Fonom[iw], BAD_FLT);
            }
            for (int32 iw = 0; iw < n_rrs; iw++) {
                ipw = ip * input->nbands + rrs_windex[iw];
                l2rec->Rrs[ipw] = inData[rrs_bindex[iw]][2 * ip] / weight;
                l2rec->nLw[ipw] = MAX(l2rec->Rrs[ipw] * l2rec->Fonom[iw], BAD_FLT);
            }
            l1rec->solz[ip] = l2rec->solz[ip] = 0.0;
            l1rec->sola[ip] = l2rec->sola[ip] = 0.0;
            l1rec->senz[ip] = l2rec->senz[ip] = 0.0;
            l1rec->sena[ip] = l2rec->sena[ip] = 90.0;
            l1rec->delphi[ip] = l2rec->delphi[ip] = 90.0;

            for (int32 iw = 0; iw < n_geo; iw++) {
                switch (geo_windex[iw]) {
                case SOLZ:
                    l1rec->solz[ip] = l2rec->solz[ip] = inData[geo_bindex[iw]][2  * ip] / weight;
                    break;
                case SOLA:
                    l1rec->sola[ip] = l2rec->sola[ip] = inData[geo_bindex[iw]][2  * ip] / weight;
                    break;
                case SENZ:
                    l1rec->senz[ip] = l2rec->senz[ip] = inData[geo_bindex[iw]][2  * ip] / weight;
                    break;
                case SENA:
                    l1rec->sena[ip] = l2rec->sena[ip] = inData[geo_bindex[iw]][2  * ip] / weight;
                    break;
                case RELAZ:
                    l2rec->delphi[ip] = inData[geo_bindex[iw]][2  * ip] / weight;
                    break;
                }
            }
            if (!have_solz) {
                // assume noon orbit
                lon = l1rec->lon[ip];
                lat = l1rec->lat[ip];
                gmt = 12 - lon / 15;
                sunangs_(&syear, &sday, &gmt, &lon, &lat, &solz, &sola);
                l1rec->solz[ip] = l2rec->solz[ip] = solz;
                l1rec->sola[ip] = l2rec->sola[ip] = sola;
            }
            if (have_sola && have_sena && !have_relaz) {
                l2rec->delphi[ip] = l2rec->sena[ip] - 180.0 - l2rec->sola[ip];
                if (l2rec->delphi[ip] < -180)
                    l2rec->delphi[ip] += 360.0;
            } else if (!have_sola && !have_sena && have_relaz) {
                l1rec->sena[ip] = l2rec->sena[ip] = 90.0;
                l1rec->sola[ip] = l2rec->sola[ip] = l2rec->sena[ip] - 180.0
                        - l2rec->delphi[ip];
                if (l2rec->sola[ip] < -180)
                    l2rec->sola[ip] += 360.0;
            }
        }

        // Add ancillary, Rayleigh, etc.
//        if (!input->ocrvc_opt)
            loadl1(&ifile, input, l1rec);
        // clear out the masks
        for (ip = 0; ip < ext; ip++)
            l2rec->mask[ip] = 0;
        cpl1l2(l1rec, l2rec);

        // Add a default chl

        for (ip = 0; ip < ext; ip++)
            l2rec->chl[ip] = get_default_chl(l2rec, &l2rec->Rrs[ip * input->nbands]);

        // Add default inherent optical properties

        if (input->iop_opt > 0 && (input->proc_ocean != 0))
            get_iops(l2rec, input->iop_opt);

        // Compute output products and store in temporary buffer

        for (iprod = 0; iprod < nprods_out; iprod++) {

            // get the product index record
            for (i=0;i<ext;i++)
                tmpData[iprod][i] = BAD_FLT;

            l2prodstr *p;

            if ((p = get_l2prod_index(ofile.l2_prod_names[iprod],
                    l1rec->sensorID, l1rec->nbands + l1rec->nbandsir, ext,
                    l1rec->nscans, l1rec->iwave)) == NULL) {
                printf("-E- %s line %d: product index failure.\n", __FILE__,
                        __LINE__);
                exit(1);
            };

            // Compute or extract the product & copy to output buffer

            VOIDP pbuf = prodgen(p, l2rec);

            memcpy(tmpData[iprod], (float *) pbuf,
                    ext * sizeof(float));

            // Check flags and set masking

            for (ip = 0; ip < ext; ip++) {
                if ((tmpData[iprod][ip] == BAD_FLT)
                        || (l2rec->flags[ip] & FLAGMASK) != 0) {
                    outMask[ip] = 1;
                }
            }

        }

        /* Write output */
        /* ------------ */
        output_binfile->clear_binlist();
        nwrite = 0;

        for (ip = 0; ip < ext; ip++) {
            if (!outMask[ip]) {
                atLeastOne = true;
                bin_num = input_binfile->get_bin_num( ip);
                offset = bin_num - basebin;
                output_binfile->set_bin_num( offset, bin_num);
                nobs = input_binfile->get_nobs( ip);
                output_binfile->inc_nobs( offset, nobs);
                nscenes = input_binfile->get_nscenes( ip);
                output_binfile->inc_nscenes( offset, nscenes);
                output_binfile->set_weights( offset, 1);

                /* Loop over data products */
                /* ----------------------- */
                for (iprod = 0; iprod < nprods_out; iprod++) {
                    outData[iprod][2*nwrite] = tmpData[iprod][ip];

                } /* iprod loop */
                if (nwrite != offset)
                    output_binfile->copy_binlist( offset, nwrite);

                nwrite++;
            }
        } /* ip loop */


        /* Write BinList & Data Products */
        /* ----------------------------- */
        if (nwrite > 0) {
            output_binfile->writeBinList( nwrite);
            for (iprod = 0; iprod < nprods_out; iprod++) {
                strcpy(buf, ofile.l2_prod_names[iprod]);
                output_binfile->writeSums(outData[iprod], nwrite, buf);
            }

            if (strcmp(input->oformat, "HDF5") == 0
                    || strcmp(input->oformat, "netCDF4") == 0)
                output_binfile->incNumRec(nwrite);
        }
    }

    output_binfile->copymeta( 1, &input_binfile);
    if(input->ocrvc_opt) {
        strcpy(output_binfile->meta_l3b.sensor_name, "OCRVC");
        strcpy(output_binfile->meta_l3b.sensor, "OCRVC");
        strcpy(output_binfile->meta_l3b.mission, "OCRVC");
    } else {
        if (strcmp(input->oformat, "netCDF4")){
            strcpy(output_binfile->meta_l3b.sensor, instrumentName[ifile.sensorID]);
        } else {
            strcpy(output_binfile->meta_l3b.sensor_name, sensorName[ifile.sensorID]);
        }
        strcpy(output_binfile->meta_l3b.mission, platformName[ifile.sensorID]);
    }
    strcpy( output_binfile->meta_l3b.infiles, basename((char *)input->ifile[0]));
    strcpy( output_binfile->meta_l3b.soft_name, "l3gen");
    sprintf(soft_id, "%d.%d.%d-r%d",L2GEN_VERSION_MAJOR,L2GEN_VERSION_MINOR,L2GEN_VERSION_PATCH_LEVEL,SVN_REVISION);
    strcpy( output_binfile->meta_l3b.soft_ver, soft_id);
    if (strcmp(input->oformat, "netCDF4") == 0){
        ptime  = unix2isodate(now(),'G');
        strcpy( output_binfile->meta_l3b.ptime, ptime);
    } else {
        strcpy( output_binfile->meta_l3b.ptime, ptime);

    }
    strcpy( output_binfile->meta_l3b.proc_con, proc_con);

    // Close files and free memory
    if (atLeastOne == false){
            cout << "No valid bins to output!" << endl << "...seems a selected product may have resulted in 100% PRODFAIL..."<< endl;
        }
    output_binfile->close();
    input_binfile->close();

    free(fullprodlist);

    for (i = 0; i < nprods_in; i++)
        free(inData[i]);

    for (i = 0; i < nprods_out; i++) {
        free(outData[i]);
        free(tmpData[i]);
    }

    free_l1(l1rec);
    free_l2(l2rec);
    free(l1rec);
    free(l2rec);
    free(input);
    free(rrs_bindex);
    free(rrs_windex);
    free(nlw_bindex);
    free(nlw_windex);

    cout << "Processing Complete at " << ptime << endl << endl;

    return(EXIT_SUCCESS);
}

