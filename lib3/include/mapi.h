/************************************************************************
*			mapi.h utilities header file			*
*				PROTOTYPE				*    	
*************************************************************************

*************************************************************************
*!C-INC
*
*!Description: Utilities header file containing M-API utility prototypes
*              and constants for MODIS science software.
*              The Header file mapi.h is part of a larger software
*              system called the MODIS Applications Programming Interface (API)
*              Utility, abbreviated M-API.  The M-API Utility consists of
*              subroutines which allow MODIS Science Team-supplied software
*              to read in Level 1B radiance bands and write out output
*              products and metadata to HDF files.  The functionality of the
*              M-API is defined in the MODIS API User's Guide, Version 2.3,
*              dated 6/18/98.
*
*	       The mapi.h file contains the prototypes, macros, and data
*              types needed to execute the MODIS API utility functions.
*
*!Input Parameters:      none
*
*!Output Parameters:     none
*
*!Revision History: 
* $Log: mapi.h,v $
* Revision 6.3  2010/07/13 14:02:47  kuyper
* Removed const qualifiers from multi-dimensional arrays, where they don't work.
*
* Revision 6.2  2010/05/12 20:04:05  kuyper
* Removed two lines that should not have been present.
*
* Revision 6.1  2010/05/04 18:51:16  kuyper
* Added ProductionHistory.
*
* Revision 5.1  2005/04/11 15:04:51  vlin
* pointer arguments to a M-API function that points to a
* read-only array of data are declared 'const'.
* vlin@saicmodis.com
*
*
*!Team-unique Header: 
*
*      This software is developed by the MODIS Science Data Support
*      Team for the National Aeronautics and Space Administration,
*      Goddard Space Flight Center, under contract NAS5-32373.
*
*!END
**************************************************************************/


#ifndef mapi_H
#define mapi_H

#include "PGS_MET.h"
#include "hdf.h"
#include "mfhdf.h"

#undef   VOIDP
/*#undef   ATTRIBUTE*/

#define MAX_ECS_NAME_L     H4_MAX_NC_NAME
#define INVENTORY_METADATA 1
#define ARCHIVED_METADATA  2

/********************************TYPEDEFs*************************************/
typedef struct dataid
{	int32		id;		/* object's identifier */
	char		*name;		/* object's name */
	char		*group;		/* Vgroup name containing the object */
	void		*info;		/* structure containing info of the object */
	struct		dataid	*next;	/* pointer to next DATAID structure */
	struct		dataid	*prev;	/* pointer to previous DATAID structure */
}DATAID;

typedef struct datainfo
{	int		nsds;		/* number of SDSs in the ring super structure */
	int		nvd;		/* number of Vdatas in the ring super structure */
	struct		dataid	*sds;	/* pointer to SDS DATAID structure */
	struct		dataid	*vd;	/* pointer to Vdata DATAID structure */
}DATAINFO;

typedef struct	modfile  /* File handle structure for HDF files*/
{
	char		*filename;	/*string name of the file*/
	int32		access;		/*Type of access to HDF file (e.g.*/
					/*DFACC_READ, RDWR, or CREATE*/
	int32		sd_id;		/*SD interface file id.	*/
	int32		hdf_id;		/*HDF file identifier	*/
	DATAINFO	*dinfo;		/* structure containing data info */
} MODFILE;

typedef char ECSattr_names_for_all_handles[PGSd_MET_NUM_OF_GROUPS][MAX_ECS_NAME_L];

/*********************************M-API DEFINITIONS***************************/

#define MFAIL -1
#define MAPIOK 0
#define DATATYPELENMAX		8   /* Maximum number of characters \
	required for a MODIS data type string, including comma or \0. */

/*********************************DATA TYPE CONSTANTS*************************/
#define	I8	"int8"
#define	I16	"int16" 
#define	I32	"int32"
#define	I64	"int64"
#define	R32	"float32"
#define	R64	"float64"
#define	TXT	"char *"
#define UI8     "uint8"
#define UI16    "uint16"
#define UI32    "uint32"
#define UI64    "uint64"

#define MODIS_ALL_TYPES 0
#define MODIS_ARRAY	720
#define MODIS_TABLE	1962

/********************************METADATA*************************************/
/****************** HDF Standard Attribute Names *****************************/

/* HDF standard SDS array structure label string (SDgetdatastrs)
	 and dimension label string (SDgetdimstrs) attribute name */
#define MLONG_NAME	"long_name"

/* HDF standard SDS array structure units string (SDgetdatastrs) 
	and dimension units string (SDgetdimstrs) attribute name */
#define MUNITS		"units"

/* HDF standard SDS array structure format string (SDgetdatastrs) 
	and dimension format string (SDgetdimstrs) attribute name */
#define MFORMAT		"format"

/* HDF standard SDS array structure coordinate system string (SDgetdatastrs) 
	attribute name */
#define MCOORD_SYS	"cordsys"

/* HDF Standard SDS array structure Calibration factor (SDgetcal)
	attribute name */
#define MSLOPE		"scale_factor"

/* HDF Standard SDS array structure Calibration factor error (SDgetcal)
	attribute name */
#define MSLOPE_ERROR	"scale_factor_err"

/* HDF Standard SDS array structure uncalibrated offset (SDgetcal)
	attribute name */
#define MOFFSET		"add_offset"

/* HDF Standard SDS array structure uncalibrated offset error (SDgetcal)
	attribute name */
#define MOFFSET_ERROR	"add_offset_err"

/* HDF Standard SDS array structure uncalibrated data HDF number type (SDgetcal)
	attribute name */
#define MNUM_TYPE	"calibrated_nt"

/* HDF standard data valid range (SDgetrange)[minimum,maximum] attribute name */
#define MDATA_RANGE	"valid_range"

/* HDF Standard SDS array structure Fill Value (SDgetfillvalue)  
	attribute name */
#define MFILL_VALUE	"_FillValue"

/****************** ECS core metadata names **********************************/

/* ECS core metadata global attribute name */
#define MECS_CORE	"CoreMetadata.0"

/* The following metadata are encoded in the "CoreMetadata.0" attribute string
 * in the ECS-defined metadata language.				     
 */

/* OBJECT = AdditionalAttributeName 
 * DEFINITION = Product Specific Attributes.
 */
#define MCORE_ADDATTRIBUTENAME         	"ADDITIONALATTRIBUTENAME"

/* OBJECT = ANCILLARYINPUTPOINTER
 * DEFINITION = References to all ancillary input files, i.e. all input files
 * other than MODIS products.
 */
#define MCORE_ANCIL_POINTER         	"ANCILLARYINPUTPOINTER"

/* OBJECT = ANCILLARYINPUTTYPE
 * DEFINITION = 
 */
#define MCORE_ANCIL_INPUT_TYPE         	"ANCILLARYINPUTTYPE"

/* OBJECT = AssociatedPlatformShortName 
 * DEFINITION = obsolete 
 */
#define MCORE_APSHORTNAME    "ASSOCIATEDPLATFORMSHORTNAME"

/* OBJECT = AssociatedInstrumentShortName 
 * DEFINITION = obsolete
 */
#define MCORE_AISHORTNAME        "ASSOCIATEDINSTRUMENTSHORTNAME"

/* OBJECT = ASSOCIATEDSENSORSHORTNAME 
 * DEFINITION = obsolete
 */
#define MCORE_ASSHORTNAME    "ASSOCIATEDSENSORSHORTNAME"

/* OBJECT = AUTOMATICQUALITYFLAG
 * DEFINITION = Indicates the results of QA performed during product
 * generation.
 */
#define MCORE_AUTO_QUALITY               "AUTOMATICQUALITYFLAG"

/* OBJECT = AUTOMATICQUALITYFLAGEXPLANATION
 * DEFINITION = 
 */
#define MCORE_AUTO_QUALITY_FLG            "AUTOMATICQUALITYFLAGEXPLANATION"

/* OBJECT = DAYNIGHTFLAG 
 * DEFINITION = 
 */
#define MCORE_DAYNIGHTFLAG        	"DAYNIGHTFLAG"
 
/* OBJECT = EASTBOUNDINGCOORDINATE
 * DEFINITION = Easternmost longitude of the granule spatial coverage.
 */
#define MCORE_EAST_BOUND         	"EASTBOUNDINGCOORDINATE"

/* OBJECT = EQUATORCROSSINGLONGITUDE 
 * DEFINITION = Easternmost longitude of the granule spatial coverage.
 */
#define MCORE_EQUATCROSSINGLONG              "EQUATORCROSSINGLONGITUDE"

/* OBJECT = EQUATORCROSSINGDATE 
 * DEFINITION = Easternmost longitude of the granule spatial coverage.
 */
#define MCORE_EQUATCROSSINGDATE        "EQUATORCROSSINGDATE"

/* OBJECT = EQUATORCROSSINGTIME 
 * DEFINITION = Easternmost longitude of the granule spatial coverage.
 */
#define MCORE_EQUATCROSSINGTIME          "EQUATORCROSSINGTIME"

/* OBJECT = EXCLUSIONGRINGFLAG.1
 * DEFINITION = Flag indicating whether points are on an inner (exclusion)
 * G-ring.
 */
#define MCORE_EXCLUS_GRING_FLG             "EXCLUSIONGRINGFLAG.1"

/* OBJECT = GRANULEPOINTER
 * DEFINITION = 
 *  obsolete 
 */
#define MCORE_GRAN_POINTER                  "GRANULEPOINTER"
 
/* OBJECT = GRINGPOINTLATITUDE.1
 * DEFINITION = Latitudes of a series of points representing the perimeter of
 * the granule spatial coverage (i.e., corners).
 */
#define MCORE_GRING_POINT_LAT            	"GRINGPOINTLATITUDE.1"

/* OBJECT = GRINGPOINTLONGITUDE.1
 * DEFINITION = Longitudes of a series of points representing the perimeter of
 * the granule spatial coverage.
 */
#define MCORE_GRING_POINT_LON                  	"GRINGPOINTLONGITUDE.1"

/* OBJECT = GRINGPOINTSEQUENCENO.1
 * DEFINITION = Sequence numbers corresponding to perimeter latitudes and
 * longitudes.
 */
#define MCORE_GRING_POINT_NUM                	"GRINGPOINTSEQUENCENO.1"
 
/* OBJECT = INPUTPOINTER
 * DEFINITION = References to other MODIS product granules used as input for
 * this product.
 */
#define MCORE_INPUT_POINTER	                "INPUTPOINTER"

/* OBJECT = INSTRUMENTSHORTNAME
 * DEFINITION =
 */
#define MCORE_INSTRUMENTSHORTNAME         "INSTRUMENTSHORTNAME"
 
/* OBJECT = LocalGranuleID 
 * DEFINITION = 
 */
#define MCORE_LOCALGRANULEID          	"LOCALGRANULEID"
 
/* OBJECT = LocalVersionID 
 * DEFINITION = 
 */
#define MCORE_LOCALVERSIONID         	"LOCALVERSIONID"

/* OBJECT = NORTHBOUNDINGCOORDINATE
 * DEFINITION = Northernmost latitude of the granule spatial coverage.
 */
#define MCORE_NORTH_BOUND	"NORTHBOUNDINGCOORDINATE"

/* OBJECT = OPERATIONMODE
 * DEFINITION = MODIS mode of operation.
 * obsolete   
 */
#define MPROD_OPERATIONMODE    		  "OPERATIONMODE"

/* OBJECT = OPERATIONALQUALITYFLAG
 * DEFINITION = The granule level flag applying both generally to the granule
 * and specifically to the parameters at the granule level. When applied to a
 * parameter, the flag refers to the quality of that parameter in the granule.
 */
#define MCORE_OPER_QUAL_FLAG   		 "OPERATIONALQUALITYFLAG"
#define MCORE_OPERATIONALQUALITYFLAG      MCORE_OPER_QUAL_FLAG 

/* OBJECT = OPERATIONALQUALITYFLAGEXPLANATION
 * DEFINITION = 
 */
#define MCORE_OQFLG_EXPL	 "OPERATIONALQUALITYFLAGEXPLANATION"

/* OBJECT = ORBITNUMBER
 * DEFINITION = Number of satellite orbit during which the granule data were
 * collected.
 */
#define MCORE_ORBIT_NUM			"ORBITNUMBER"

/* OBJECT = PARAMETERNAME 
 * DEFINITION = 
 */
#define MCORE_PARAMETERNAME		"PARAMETERNAME"

/* OBJECT = PARAMETERVALUE 
 * DEFINITION = 
 */
#define MCORE_PARAMETERVALUE		"PARAMETERVALUE"

/* OBJECT = PGEVERSION 
 * DEFINITION = 
 */
#define MCORE_PGEVERSION		"PGEVERSION"

#define MCORE_PROCESSING_ENV            "PROCESSINGENVIRONMENT"
/* OBJECT = PLATFORMSHORTNAME 
 * DEFINITION = 
 */
#define MCORE_PLATFORMSHORTNAME		"PLATFORMSHORTNAME"

/* OBJECT = PRODUCTIONDATETIME 
 * DEFINITION = Reference to processing history file.
 */
#define MCORE_PRODUCTIONDATETIME	"PRODUCTIONDATETIME"

/* OBJECT = PROCESSINGHISTORYPOINTER
 * DEFINITION = Reference to processing history file.
 * obsolete
 */
#define MCORE_HISTORY_POINTER		"PROCESSINGHISTORYPOINTER"

/* OBJECT = PRODUCTIONHISTORY
 * DEFINITION = Concatenation of PGE versions used to create this file and
 *	it's input files.
 */
#define MCORE_PRODUCTIONHISTORY		"PRODUCTIONHISTORY"

/* OBJECT = QAPERCENTCLOUDCOVER 
 * DEFINITION = Value indicating the percent of interpolated data in the
 * granule
 */
#define MCORE_QAPERCENTCLOUDCOVER	"QAPERCENTCLOUDCOVER"

/* OBJECT = QAPERCENTINTERPOLATEDDATA
 * DEFINITION = Value indicating the percent of interpolated data in the
 * granule
 */
#define MCORE_PERCENT_INTERP		"QAPERCENTINTERPOLATEDDATA"

/* OBJECT = QAPERCENTMISSINGDATA
 * DEFINITION = Value indicating the percent of missing data in the granule.
 */
#define MCORE_PERCENT_MISSING		"QAPERCENTMISSINGDATA"

/* OBJECT = QAPERCENTOUTOFBOUNDSDATA
 * DEFINITION = Value indicating the percent of data in the granule outside of
 * acceptable limits.
 */
#define MCORE_PERCENT_OUT		"QAPERCENTOUTOFBOUNDSDATA"

/* OBJECT = QUALITYFLAGEXPLANATION
 * DEFINITION = A text explanation of the criteria used to set each quality
 * lag; including thresholds or other criteria.
 * obsolete
 */
#define MCORE_QUAL_EXPL        	        "QUALITYFLAGEXPLANATION"
#define MCORE_QUALITYFLAGEXPLANATION    MCORE_QUAL_EXPL 

/* OBJECT = RANGEBEGINNINGDATETIME
 * DEFINITION = The date and time when the temporal coverage period of this
 * granule began. obsolete
 */
#define MCORE_RANGE_START		"RANGEBEGINNINGDATETIME"

/* OBJECT = RANGEBEGINNINGDATE
 * DEFINITION = The time when the temporal coverage period of this
 * granule began.
 */
#define MCORE_RANGE_BEG_DATE		"RANGEBEGINNINGDATE"

/* OBJECT = RANGEBEGINNINGTIME
 * DEFINITION = The date when the temporal coverage period of this
 * granule began.
 */
#define MCORE_RANGE_BEG_TIME		"RANGEBEGINNINGTIME"

/* OBJECT = RANGEENDINGDATETIME
 * DEFINITION = The date and time when the temporal coverage period of this
 * granule ended.  obsolete
 */
#define MCORE_RANGE_END			"RANGEENDINGDATETIME"

/* OBJECT = RANGEENDINGDATE
 * DEFINITION = The date and time when the temporal coverage period of this
 * granule ended.
 */
#define MCORE_RANGE_ENDING_DATE		"RANGEENDINGDATE"

/* OBJECT = RANGEENDINGTIME
 * DEFINITION = The date and time when the temporal coverage period of this
 * granule ended.
 */
#define MCORE_RANGE_ENDING_TIME		"RANGEENDINGTIME"

/* OBJECT = TimeofDay
 * DEFINITION = Specifies a single point in time covered by a data collection,
 * granule or event.
 */
#define MCORE_TIME_OF_DAY		"TimeofDay"

/* OBJECT = CalendarDate
 * DEFINITION = Specifies a single date covered by a data collection, granule,
 * or event.
 */
#define MCORE_CALENDAR_DATE		"CalendarDate"

/* OBJECT = REPROCESSINGPLANNED
  * DEFINITION = Indicator of what reprocessing is planned for the granule. 
  */
#define MCORE_TO_BE_REDONE		"REPROCESSINGPLANNED"

/* OBJECT = REPROCESSINGACTUAL
 * DEFINITION = Indicator of the reprocessing status of the granule.
 */
#define MCORE_ACTUALLY_REDONE		"REPROCESSINGACTUAL"

/* OBJECT = SCIENCEQUALITYFLAG
 * DEFINITION = The granule level flag applying both generally to the granule
 * and specifically to the parameters at the granule level. When applied to a
 * parameter, the flag refers to the quality of that parameter in the granule.
 */
#define MCORE_SCIENCE_QUAL_FLG 		"SCIENCEQUALITYFLAG"
#define MCORE_SCIENCEQUALITYFLAG 	MCORE_SCIENCE_QUAL_FLG 

/* OBJECT = SCIENCEQUALITYFLAGEXPLANATION
 * DEFINITION = The granule level flag applying both generally to the granule
 * and specifically to the parameters at the granule level. When applied to a
 * parameter, the flag refers to the quality of that parameter in the granule.
 */
#define MCORE_SCIENCE_QUAL_FLG_EXPL 	       "SCIENCEQUALITYFLAGEXPLANATION"

/* OBJECT = SENSORSHORTNAME
 * DEFINITION = 
 * obsolete
 */
#define MCORE_SENSOR_SHORT_NAME                  "SENSORSHORTNAME"

/* OBJECT = SensorCharacteristicName 
 * DEFINITION = 
   obsolete
 */
#define MCORE_SENSORCHARACTERISTICNAME        "SENSORCHARACTERISTICNAME"

/* OBJECT = SensorCharacteristicValue 
 * DEFINITION =
   obsolete
 */
#define MCORE_SENSORCHARACTERISTICVALUE        "SENSORCHARACTERISTICVALUE"

/* OBJECT = SHORTNAME
 * DEFINITION = The identifier for the data collection.
 */
#define MCORE_SHORT_NAME      		  "SHORTNAME"
#define MCORE_SHORTNAME      		  MCORE_SHORT_NAME 

/* OBJECT = SIZEMBECSDATAGRANULE
 * DEFINITION = The size of the data granule in megabytes.
 */
#define MCORE_SIZE_OF_GRANULE            "SIZEMBECSDATAGRANULE"

/* OBJECT = SOUTHBOUNDINGCOORDINATE
 * DEFINITION = Southernmost latitude of the granule spatial coverage.
 */
#define MCORE_SOUTH_BOUND	"SOUTHBOUNDINGCOORDINATE"

/* OBJECT = VERSIONID 
 * DEFINITION = 
 */
#define MCORE_VERSIONID			"VERSIONID"

/* OBJECT = WESTBOUNDINGCOORDINATE
 * DEFINITION = Westernmost longitude of the granule spatial coverage.
 */
#define MCORE_WEST_BOUND	"WESTBOUNDINGCOORDINATE"

/****************** ECS product specific metadata ****************************/

/* ECS product specific  metadata global attribute name */
#define MECS_PRODUCT	MECS_ARCHIVE
#define MECS_ARCHIVE	"ArchiveMetadata.0"

/****************** ECS product inventory metadata ****************************/

/* OBJECT = MODISPRODUCTFILENAME
 * DEFINITION = The MODIS filename for this granule.
 * obsolete 
 */
#define MPROD_FILENAME 	          "MODISPRODUCTFILENAME"
#define MPROD_PRODUCTFILENAME   MPROD_FILENAME 

/* OBJECT = PROCESSINGDATETIME
 * DEFINITION = This field contains the date and time the process that
 * created this file was started.
 */
#define MPROD_PROC_DATE_TIME    "PROCESSINGDATETIME"
#define MPROD_PROCESSINGDATETIME MPROD_PROC_DATE_TIME

/* OBJECT = SPSOPARAMETERS
 * DEFINITION = The SPSO parameters for all data contained in this file, as
 * listed in the SPSO database.
 */
#define MPROD_SPSO_PARAM        "SPSOPARAMETERS"
#define MPROD_SPSOPARAMETERS    MPROD_SPSO_PARAM 

/****************** ECS product inventory-a metadata ****************************/

/* OBJECT = GRANULENUMBER
DEFINITION = The number of this MODIS granule.  This is a "PSA".
*/
#define MPROD_GRANULE_NUM       "GRANULENUMBER"
#define MPROD_GRANULENUMBER     MPROD_GRANULE_NUM 

/****************** ECS product inventory-b metadata ****************************/
/* OBJECT = GRIDTYPE
 * DEFINITION = The type of MODLAND global grid used for L2G, L3 and L4
 * products.
 */
#define MPROD_GRID_TYPE                "GRIDTYPE"
#define MPROD_GRIDTYPE                 MPROD_GRID_TYPE 

/* OBJECT = HORIZONTALTILENUMBER
 * DEFINITION = The horizontal tile number of this tile in the MODLAND
 * integerized sinusoidal grid.
 */
#define MPROD_HORIZ_TILE_NUM         "HORIZONTALTILENUMBER"
#define MPROD_HORIZONTALTILENUMBER    MPROD_HORIZ_TILE_NUM 

/* OBJECT = VERTICALTILENUMBER
 * DEFINITION = The vertical tile number of this tile in the MODLAND
 * integerized sinusoidal grid.
 */
#define MPROD_VERT_TILE_NUM   	      "VERTICALTILENUMBER"
#define MPROD_VERTICALTILENUMBER       MPROD_VERT_TILE_NUM 

/****************** ECS product-archive metadata ****************************/

/* OBJECT = ALGORITHMPACKAGEACCEPTANCEDATE
 * DEFINITION = The date this algorithm package version successfully passed
 * AI&T procedures and was accepted as an ECS standard algorithm.
 "*/
#define MPROD_ALGO_PCK_ACPT_DATE         "ALGOPACKAGEACCEPTANCEDATE"
#define MPROD_ALGOPACKAGEACCEPTANCEDATE  MPROD_ALGO_PCK_ACPT_DATE 

/* OBJECT = ALGORITHMPACKAGEMATURITYCODE
 * DEFINITION = This specifies the maturity of the algorithm package.
 */
#define MPROD_ALGO_PACK_MAT_CODE         "ALGOPACKAGEMATURITYCODE"
#define MPROD_ALGOPACKAGEMATURITYCODE     MPROD_ALGO_PACK_MAT_CODE 

/* OBJECT = ALGORITHMPACKAGENAME
 * DEFINITION = The name given to the complete delivered algorithm package
 * submitted for algorithm integration and test.
 */
#define MPROD_ALGO_PACK_NAME         "ALGORITHMPACKAGENAME"
#define MPROD_ALGORITHMPACKAGENAME    MPROD_ALGO_PACK_NAME 

/* OBJECT = ALGORITHMPACKAGEVERSION
 * DEFINITION = The version of the algorithm package.
 */ 
#define MPROD_ALGO_PACK_VER              "ALGORITHMPACKAGEVERSION" 
#define MPROD_ALGORITHMPACKAGEVERSION     MPROD_ALGO_PACK_VER 

/* OBJECT = INSTRUMENTNAME
 * DEFINITION = The long name by which the instrument is known.
 */
#define MPROD_INSTR_NAME                  "INSTRUMENTNAME"
#define MPROD_INSTRUMENTNAME               MPROD_INSTR_NAME                  
 
/* OBJECT = LocalInputGranuleID 
 * DEFINITION = 
 */
#define MCORE_LOCALINPUTGRANULEID          "LOCALINPUTGRANULEID"

/* OBJECT = LONGNAME
 * DEFINITION = A descriptive name for the data collection.
 */
#define MCORE_LONG_NAME		    "LONGNAME"

/* BJECT = PLATFORMSHORTNAME
 * DEFINITION = The short name assigned to the platform carrying the
 "* instrument.
 */ 
#define MPROD_PLATFORM_SHORT_NAM    "PLATFORMSHORTNAME"
#define MPROD_PLATFORMSHORTNAME      MPROD_PLATFORM_SHORT_NAM 
 
/* OBJECT = PROCESSINGCENTER
 * DEFINITION = DAAC where product is processed.
 */
#define MPROD_PROC_CENTER            "PROCESSINGCENTER" 
#define MPROD_PROCESSINGCENTER        MPROD_PROC_CENTER 

/********************************UTILITY PROTOTYPES***************************/

			int closeMODISfile(	MODFILE 	**file);

			int getMODISardims(	MODFILE		*file,
						char const	*arrayname,
						char const 	*groupname,
						char 		*data_type,
						int32		*rank,
						int32		dimsizes[ ]);

			int getMODISarray(	MODFILE		*file,
						char const	*arrayname,
						char const	*groupname,
						int32	 	start[ ],
						int32		dimsizes[ ],
						void 		*data);

			int32 MODISsizeof(	char const	*data_type);

	                MODFILE *openMODISfile(	char const	*filename,
						char const 	*access);

			int putMODISarray(	MODFILE		*file,
						char const	*arrayname,
						char const	*groupname,
						int32 const	start[ ],
						int32 const	dimsizes[ ],
						void const      *data);

  			int createMODISarray( MODFILE		*file,
						char const 	*arrayname,
						char const 	*groupname,
						char const 	*data_type,
						int32		rank,
						int32 const  dimsizes[] );

			int createMODISgroup(	MODFILE		*file,
						char const 	*groupname,
						char const 	*classname);

			int createMODIStable(	MODFILE		*file,
						char const	*tablename,
						char const 	*classname,
						char const 	*groupname,
					        char const      *fieldname,
						char const	*data_type);

			int getMODISfields(	MODFILE		*file,
						char const	*tablename,
						char const 	*groupname,
 						int32		*stringlen,
 						int32		*recno,
						int32		*fieldno,
						char		*fieldname,
						char		*data_type,
						char		*classname);

			int getMODISfileinfo(	MODFILE		*file,
						char const	*attribute,
                                                char            *data_type,
                                                int32	        *n_elements,
                                                void            *value);

			int getMODIStable(	MODFILE		*file,
						char const	*tablename,
						char const 	*groupname,
						char const      *fieldname,
						int32		start,
						int32		recno,
						int32		*buffsize,
						unsigned char 	*data);

			int putMODISfileinfo(	MODFILE		*file,
						char const 	*attribute,
						char const 	*data_type,
						int32       n_elements,
						void const  	*value);

                	int putMODIStable(	MODFILE		*file,
						char const	*tablename,
						char const 	*groupname,
						int32		start,
						int32		recno,
						unsigned char const *data);  

                        int getMODISarinfo(     MODFILE         *file,
						char const	*arrayname,
						char const 	*groupname,
                                                char const      *attribute,
                                                char            *data_type,
                                                int32	        *n_elements,
                                                void            *value);

                        int putMODISarinfo(     MODFILE         *file,
						char const	*arrayname,
						char const 	*groupname,
                                                char const      *attribute,
                                                char const      *data_type,
                                                int32	        n_elements,
                                                void const      *value);

                        int putMODISdiminfo(    MODFILE         *file,
						char const	*arrayname,
						char const 	*groupname,
                                                int32	        dimension,
                                                char const      *attribute, 
                                                char const      *data_type,
                                                int32	        n_elements,
                                                void const      *value);
                       int getMODISdiminfo(     MODFILE         *file,
						char const	*arrayname,
						char const 	*groupname,
                                                int32	        dimension,
                                                char const      *attribute,
                                                char            *data_type,
                                                int32	        *n_elements,
                                                void            *value);

                       int getMODISECSinfo(     MODFILE         *file, 
                                                char const      *PVLAttrName, 
                                                char const      *parmName,
                                                char            *data_type,
                                                int32	        *n_elements, 
                                                void            *value); 

                       int completeMODISfile(  
			   MODFILE              **file,
			   PGSt_MET_all_handles mdHandles,
			   ECSattr_names_for_all_handles HDFattrNames,
			   int32	            NumHandles); 

                       int putMODISdimname(     MODFILE *file, 
						char const *arrayname,
                                                char const *groupname,
                                                int32 dimension,  
                                                char const *dimname);

                       int substrMODISECSinfo(  char const	*char_value, 
                                                int32		n_elements,
                                                int32		*n_strings, 
                                                char		*substr[]); 

                       int getMODISdimname(     MODFILE *file, 
						char const *arrayname,
                                                char const *groupname,
                                                int32 dimension,
                                                char *dimname);

                       int endMODISobjaccess(   MODFILE *file,
                                                char const *name,
                                                char const *group,
                                                int32 type);

                       int getMODIShobjid(      MODFILE *file, 
                                                char const *name, 
                                                char const *group,
                                                int32 type, 
                                                char const *access);

		       MODFILE *createMAPIfilehandle(int32 fid);

                       int     releaseMAPIfilehandle(MODFILE **file);

#endif
