#include "smeta_api.h"        		/* SMETA structure and routines   */

static SMETA	   	 smeta;		/* Enhanced metadata structure    */
static WRS2		  wrsparms;	    /* Nominal WRS-2 orbit parameters */
static double		dellon;		/* WRS-2 longitude offset         */
static double		delrow;		/* Framed scene size in WRS rows  */
static VECTOR		ecfsun;		/* Scene center ECEF solar vector */
static MTL_GRID_BAND	bgrid[IAS_MAX_NBANDS];	/* Geometric grid for each band */

int smeta_init( 			/* Returns num bands in SMETA     */
    char	*smeta_filename	)	/* Product metadata file name     */
{
    int band_index;			/* Loop index for bands           */

    /* Read the input MTL file */
    if ( smeta_read( &smeta, smeta_filename ) != SUCCESS )
    {
        IAS_LOG_ERROR("Reading the input metadata file %s.", smeta_filename);
        return(0);
    }

    /* Load WRS parameters  - these values are available in the CPF */
    wrsparms.semimajor = 6378137.0;
    wrsparms.semiminor = 6356752.3142;
    wrsparms.inertialvel = 7.2921158553E-05;
    wrsparms.orbitincl = 98.2;		/* Nominal orbital inclination */
    wrsparms.orbitrad = 7083137.0;	/* Nominal 705 km orbital altitude */
    wrsparms.path1long = -64.6;
    wrsparms.numpath = 233;
    wrsparms.numrow = 248;
    wrsparms.cycle = 16;
    wrsparms.row0lat = 60;

    /* Initialize the scene map projection */
    if ( smeta_init_projection( smeta.projection ) != SUCCESS )
    {
        IAS_LOG_ERROR("Initializing the map projection transformation.");
        return(0);
    }

    /* Adjust nominal scene frame to fit metadata */
    dellon = 0.0;
    delrow = 1.2;
    if ( frame_scene( wrsparms, smeta, &dellon, &delrow ) != SUCCESS )
    {
        IAS_LOG_ERROR("Fitting scene frame to product metadata.\n");
        smeta_release_projection();
        return(0);
    }
    /* Iterate to refine fit */
    if ( frame_scene( wrsparms, smeta, &dellon, &delrow ) != SUCCESS )
    {
        IAS_LOG_ERROR("Fitting scene frame to product metadata.\n");
        smeta_release_projection();
        return(0);
    }

    /* Compute ECEF solar vector */
    if ( get_ecfsun( wrsparms, smeta, 0.0, &ecfsun ) != SUCCESS )
    {
        IAS_LOG_ERROR("Calculating ECEF solar vector.\n");
        smeta_release_projection();
        return(0);
    }

    for ( band_index=0; band_index<smeta.num_band; band_index++ )
    {
        /* Generate grid for this band */
        if ( calc_band_grid( wrsparms, smeta, dellon, ecfsun, band_index, &bgrid[band_index] ) != SUCCESS )
        {
          IAS_LOG_ERROR("Creating grid for band %d.", smeta.band_smeta[band_index].band);
          smeta_release_projection();
          return(0);
        }
    }

    return(smeta.num_band);
}

char *smeta_file_name()
{
  return( smeta.scene_id );
}

int smeta_get_frame(
    int         band_index,		/* Input band index        */
    META_FRAME *frame )			/* Output image frame info */
{

    if ( band_index < 0 || band_index > smeta.num_band-1 )
    {
        IAS_LOG_ERROR("Invalid band requested");
        return(ERROR);
    }

    /* Establish the output image file frame */
    frame->band = smeta.band_smeta[band_index].band;
    frame->nlines = smeta.band_smeta[band_index].l1t_lines;
    frame->nsamps = smeta.band_smeta[band_index].l1t_samps;
    frame->code = smeta.projection.code;
    frame->zone = smeta.projection.zone;
    frame->pixsize = smeta.band_smeta[band_index].pixsize;
    frame->ul_x = smeta.projection.corners.upleft.x;
    frame->ul_y = smeta.projection.corners.upleft.y;

    return(SUCCESS);
}

int smeta_calc(
    int		in_line,		/* L1T line coordinate              */
    int		in_samp,		/* L1T sample coordinate            */
    int		in_band,		/* Spectral band number             */
    short	*sat_zn,		/* Satellite zenith angle           */
    short	*sat_az,		/* Satellite azimuth angle          */
    short	*sun_zn,		/* Solar zenith angle               */
    short	*sun_az )		/* Solar azimuth angle              */
{
    int		band_index;		/* Band index                       */
    int		status;			/* Status return flag               */

 //   printf("SMETA_CALC: %d %d %d\n",in_line, in_samp, in_band);
    /* Find the index for the current band */
    band_index = smeta_band_number_to_index( &smeta, in_band );
    if ( band_index < 0 )
    {
        IAS_LOG_ERROR("Invalid band requested: %d.", in_band);
        smeta_close();
        return(ERROR);
    }

    /* Generate angles for this pixel */
    if ( (status=smeta_angles( (double)in_line, (double)in_samp, bgrid[band_index],
                               sat_zn, sat_az, sun_zn, sun_az )) != SUCCESS )
    {
        IAS_LOG_ERROR("Generating angles for band %d.", smeta.band_smeta[band_index].band);
        smeta_close();
    }
//    printf("sat_az[%d]=%d\n",is,sat_az[is]);

    /* Return status */
    return( status );
}

void smeta_close()
{
    int band_index;			/* Index of current band */

    /* Release the projection transformation */
    smeta_release_projection();

    /* Release the band grids */
    for ( band_index=0; band_index<smeta.num_band; band_index++ )
        free( bgrid[band_index].sgrid );

    return;
}

