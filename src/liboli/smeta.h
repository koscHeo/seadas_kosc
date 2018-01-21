/*
 *	Name:	smeta.h
 *
 *	Purpose:
 *	Include file defining the product metadata structure that is used to generate the
 *	scene viewing geometry information needed to calculate per-pixel view angles.
 */
#ifndef	_SMETA_H_
#define	_SMETA_H_

/* Defines the maximum number of normal bands present on any of the supported
   satellites.  Note that this should only be used if absolutely necessary.
   It is preferred that arrays be sized dynamically. */
#define IAS_MAX_NBANDS 11

/* Defines the maximum number of SCAs present on any of the supported 
   sensors.   Note that this should only be used if absolutely necessary.
   It is preferred that arrays be sized dynamically. */
#define IAS_MAX_NSCAS 14

#include "ias_const.h"
#include "ias_structures.h"			/* IAS_VECTOR and IAS_CORNERS                */

typedef struct SMETA_SCENE_PROJ
{
    /* ==== WGS84 ellipsoid parameters ====                                                  */
    double	wgs84_major_axis;		/* WGS 84 ellipsoid semi-major axis (meters) */
    double	wgs84_minor_axis;		/* WGS 84 ellipsoid semi-minor axis (meters) */
    /* ==== Projection parameters and info. ====                                             */
    char	units[IAS_UNITS_SIZE];		/* Projection units string                   */
    int		code;				/* Projection code for the output space image.
						   Values for this field are defined in the
						   "gctp.h" include file.                    */
    char	datum[IAS_DATUM_SIZE];		/* Projection datum string                   */
    int		spheroid;			/* Projection spheroid code                  */
    int		zone;				/* Projection zone code for UTM or
						   State Plane projections.                  */
    double	projprms[IAS_PROJ_PARAM_SIZE];
		/* Array of 15 projection coefficients as required by the projection 
		   transformation package.  Refer to the projection package documentation 
		   for a description of each field for a given projection.                   */
    /* ==== Grid corners. ====                                                               */
    struct	IAS_CORNERS corners;		/* Projection coordinates of the resulting 
						   output image's four corners.              */
} SMETA_SCENE_PROJ;

typedef struct SMETA_BAND
{
    int			band;			/* User band number                          */
    int			l1t_lines;		/* Number of lines in the L1T image.         */
    int			l1t_samps;		/* Number of samples in the L1T image.       */
    int			nsca;			/* Number of SCAs in this band               */
    double		pixsize;		/* Projection distance per pixel in meters.  */
    double		align[3][3];		/* Instrument to spacecraft alignment        */
    double		legendre[IAS_MAX_NSCAS][2][4];   /* Per SCA Legendre coefficients    */
} SMETA_BAND;

typedef struct SMETA
{
    char		scene_id[22];		/* Scene ID as Lx8ppprrryyyydddLGNnn         */
    int			wrs_path;		/* WRS path number (target)                  */
    int			wrs_row;		/* WRS row number (target)                   */
    double		roll_angle;		/* Off-nadir roll angle in degrees           */
    double		sun_azimuth;		/* Scene center sun azimuth in degrees       */
    double		sun_elevation;		/* Scene center sun elevation in degrees     */
    SMETA_SCENE_PROJ	projection;		/* Ground reference for this L1T scene       */
    int			num_band;		/* Number of bands in the metadata structure */
    SMETA_BAND		band_smeta[IAS_MAX_NBANDS];	/* Metadata for each band            */
} SMETA;

#endif
