#include "smeta_api.h"       /* Product metadata functionality */
#define MAXBAND 7
/**
 * Get the Landsat Nominal OLI angles from a MTL file
 * @param in input MTL file
 * @param in input number of pixels per line
 * @param in input number of lines
 * @param in input scan line to work on
 * @param out output satellite zenith angle
 * @param out output satellite azimuth angle
 * @param out output solar zenith angle
 * @param out output solar azimuth angle
 */
int get_oli_nom_angles( char *meta_filename, int32_t npix, int32_t nscan, int32_t iscan,
        float *solz, float *sola, float *senz, float *sena)
{
    static int  firstCall = 1;
    static int  nband;                     /* Number of bands                */
    static META_FRAME   frame[MAXBAND];    /* Output image frame info        */
    char 		*root_filename;			    /* Scene root file name           */
    int         band;	  			        /* User band number               */
    int			band_index;			        /* Band index                     */
    int         is;                         /* Sample (pixel) index           */
    int			status;				        /* Status return flag             */
    short		*sat_zn;			        /* Satellite zenith angle         */
    short		*sat_az;			        /* Satellite azimuth angle        */
    short		*sun_zn;			        /* Solar zenith angle             */
    short		*sun_az;			        /* Solar azimuth angle            */

    if (firstCall == 1) {
        /* Initialize the metadata package */
        if ( (nband = smeta_init( meta_filename )) < 1 )
        {
            IAS_LOG_ERROR("Initializing the metadata file %s.", meta_filename);
            fprintf(stderr,
                    "-E- %s line %d: Initializing the metadata file %s.\n",
                    __FILE__, __LINE__, meta_filename);
            exit (EXIT_FAILURE);
        }
    }
    /* Allocate the intermediary buffers for calculating average over MAXBAND*/
    if ( (sat_zn = (short *)malloc( npix*sizeof(short) )) == NULL )
    {
        smeta_close();
        fprintf(stderr,
                "-E- %s line %d: Unable to allocate satellite elevation angle array.\n",
                __FILE__, __LINE__);
        exit (EXIT_FAILURE);
    }
    if ( (sat_az = (short *)malloc( npix*sizeof(short) )) == NULL )
    {
        free(sat_zn);
        smeta_close();
        fprintf(stderr,
                "-E- %s line %d: Unable to allocate satellite azimuth angle array.\n",
                __FILE__, __LINE__);
        exit (EXIT_FAILURE);
    }
    if ( (sun_zn = (short *)malloc( npix*sizeof(short) )) == NULL )
    {
        free(sat_zn);
        free(sat_az);
        smeta_close();
        fprintf(stderr,
                "-E- %s line %d: Unable to allocate solar elevation angle array.\n",
                __FILE__, __LINE__);
        exit (EXIT_FAILURE);
    }
    if ( (sun_az = (short *)malloc( npix*sizeof(short) )) == NULL )
    {
        free(sat_zn);
        free(sat_az);
        free(sun_zn);
        smeta_close();
        fprintf(stderr,
                "-E- %s line %d: Unable to allocate solar azimuth angle array.\n",
                __FILE__, __LINE__);
        exit (EXIT_FAILURE);
    }

    /* Initialize  */
    for ( is=0; is<npix; is++ )
    {
    senz[is] = 0;
    sena[is] = 0;
    solz[is] = 0;
    sola[is] = 0;
    }

    /* Process the angles for each band */
    for ( band_index=0; band_index<MAXBAND; band_index++ )
    {
	/* Get framing information for this band */
        if ( firstCall == 1 && smeta_get_frame( band_index, &frame[band_index] ) != SUCCESS )
        {
            fprintf(stderr,
                    "-E- %s line %d: Unable to get smeta frame.\n",
                    __FILE__, __LINE__);
            exit (EXIT_FAILURE);
        }
        band = frame[band_index].band;
//        if ( (iscan) % 500 == 0 ) printf("GET_OLI_NOM_ANGLES: Band/Line %d/%d\n", band, iscan );

        /* Loop through the L1T lines and samples */
             for ( is=0; is<npix; is++ )
            {
                status = smeta_calc( (int) iscan, is, band, &sat_zn[is], &sat_az[is], &sun_zn[is], &sun_az[is] );
                if ( status != SUCCESS )
                {
                    fprintf(stderr,
                            "-E- %s line %d: Error Evaluating angles in band %d.\n",
                            __FILE__, __LINE__,band);
                    exit (EXIT_FAILURE);
                }
                /* Average the bands and convert to decimal degrees */
             //   printf("band=%d sat_az[%d]=%d\n",band_index, is,sat_az[is]);
                senz[is] += (float)sat_zn[is]/(100*MAXBAND);
                sena[is] += (float)sat_az[is]/(100*MAXBAND);
                solz[is] += (float)sun_zn[is]/(100*MAXBAND);
                sola[is] += (float)sun_az[is]/(100*MAXBAND);
                //printf("sat_az[%d]=%d, sena[%d]=%f\n",is,sat_az[is],is,sena[is]);
            }
         }

/*    for ( is=0; is<npix; is+=10) {
        printf("sena[%d]=%f\n",is,sena[is]);
        printf("sat_az[%d]=%d\n",is,sat_az[is]);
    }*/
    /* Free the angle arrays */
    free( sat_zn );
    free( sat_az );
    free( sun_zn );
    free( sun_az );

    /* Release the metadata */
    //smeta_close();  /* Probably should close this properly somewhere after processing */
    firstCall = 0;
    return (SUCCESS);
}
