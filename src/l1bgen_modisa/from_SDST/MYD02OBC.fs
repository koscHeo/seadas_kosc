+=============================================================================+
|                                                                             |
|                L1B OBC File Specification -- MODIS/AQUA                     |
|                           V6.0.1, 05/22/2008                                |
|       Effective for all higher numbered versions until superceded           |
|                                                                             |
+=============================================================================+

   This document describes the format and content for the MODIS/AQUA Level 1B
   on-board calibrator (OBC) product. The product is implemented in an HDF
   file.

   +=====================================================================+
   |                                                                     |
   |  Contents:                                                          |
   |                                                                     |
   |       i)  Band and Detector Order Conventions                       |
   |                                                                     |
   |       I)  Global Metadata and Attributes                            |
   |                                                                     |
   |      II)  Scientific Data Sets (SDSs)                               |
   |                                                                     |
   |     III)  HDF Vdatas                                                |
   |                                                                     |
   +=====================================================================+

--------------------------------------
i) Band and Detector Order Conventions
--------------------------------------

This product contains data using various conventions in terms of band
groupings and detector order.  The purpose of this section is to define
the conventions used within this product.   In many cases, the name of the
data item will indicate the band convention implied.  The detector order
conventions of individual data objects is explicitly defined in a table
in subsection 1.2.

Band Conventions
----------------

Attributes dimensioned 38 (number of bands) refer to the full set of MODIS
bands, which has the following implied order:

       1, 2, 3, ... 12, 13lo, 13hi, 14lo, 14hi, 15, ... , 36

The set of 22 reflective solar bands (RefSB) are:

       1, 2, 3, ... 12, 13lo, 13hi, 14lo, 14hi, 15, ... , 19, 26.

The set of 16 thermal emissive bands (TEB or Emiss) are:

       20-25, 27-36

There are three product resolutions of the digital numbers: 250m, 500m, 1km.
The Level 1A band groupings by resolution are:

       250m bands:         1,2
       500m bands:         3,4,5,6,7
       1km_day bands:      8-12,13lo,13hi,14lo,14hi,15-19
       1km_night bands:    20-36

The Level 1B band groupings by resolution are:

       250m bands:         1,2
       500m bands:         3,4,5,6,7
       1km_RefSB bands:    8-12,13lo,13hi,14lo,14hi,15-19,26
       1km_Emissive bands: 20-25,27-36

Some of the global metadata and a few SDSs utilize the Level 1B band
groupings.  ALL the digital number (DN) SDSs in this product use the
Level 1A band groupings.  This is different than the Level 1B EV products,
where all data objects use the Level 1B band groupings.

Detector Order Convention for Detector-Dependent Data
-----------------------------------------------------

The "SBRS" detector order convention refers to the numbering of detectors
on each focal plane assembly as defined by the manufacturer.  In this
convention, the numbering of detectors is opposite in direction from the
satellite track direction.  The "product" detector order convention is the
reverse of the SBRS convention.  The product order was defined to allow
consecutive scans to be laid down on an image with the last detector of one
scan abutting the first detector of the next scan. The following table
summarizes the detector ordering convention implied for those data items in
this product that have a dimension that is detector dependent.

   Data Name                                          Detector Order
 -------------                                        --------------
 SD_250m                                                 product
 SD_500m                                                 product
 SD_1km_day                                              product
 SD_1km_night                                            product
 SRCA_250m                                               product
 SRCA_500m                                               product
 SRCA_1km_day                                            product
 SRCA_1km_night                                          product
 BB_250m                                                 product
 BB_500m                                                 product
 BB_1km_day                                              product
 BB_1km_night                                            product
 SV_250m                                                 product
 SV_500m                                                 product
 SV_1km_day                                              product
 SV_1km_night                                            product
 DN_obc_avg_250m                                         product
 DN_obc_var_250m                                         product
 DN_obc_outlier_mask_250m                                product
 DN_obc_avg_500m                                         product
 DN_obc_var_500m                                         product
 DN_obc_outlier_mask_500m                                product 
 DN_obc_avg_1km_day                                      product
 DN_obc_var_1km_day                                      product
 DN_obc_outlier_mask_1km_day                             product
 DN_obc_avg_1km_night                                    product
 DN_obc_var_1km_night                                    product
 DN_obc_outlier_mask_1km_night                           product 
 Noise in Thermal Detectors                              product
 Change in relative responses of thermal detectors       product
 DC Restore Change for Thermal Bands                     product
 DC Restore Change for Reflective 250m Bands             product
 DC Restore Change for Reflective 500m Bands             product
 DC Restore Change for Reflective 1km Bands              product
 Dead Detector List                                      product
 Noisy Detector List                                     product
 Dead Subframe List                                      product
 Noisy Subframe List                                     product
 Detector Quality Flag                                   product
 Detector Quality Flag2                                  product
 fpa_dcr_offset                                          SBRS
 raw_pv_gains                                            SBRS

+=============================================================================+

I)  Global Metadata and Attributes

This section describes granule-level metadata (i.e. ECS metadata and product
granule metadata).  Many of these data are the same as in the Level 1B 1km EV
product and thus represent the metatadata for the earth view sector.  These
are included here in the OBC file for convenience only.

Contents:
  1.1) ECS Standard Core Granule Metadata
  1.2) ECS Standard Archive Granule Metadata
  1.3) Product Granule Metadata Stored as Global Attributes

+=============================================================================+
|                                                                             |
|  1.1) ECS Standard Core Granule Metadata                                    |
|                                                                             |
|  Stored as one ECS PVL string in global attribute: "CoreMetadata.0".        |
|  The example indicates the number of values and data type of                |
|  the particular item.  If there is any question, see the                    |
|  corresponding MCF file, MOD02OBC.mcf.                                      |
|                                                                             |
+=============================================================================+
|                                                                             |
| Name                                        Example                         |
+-----------------------------------------------------------------------------+
| LOCALGRANULEID                              "MYD02OBC.A2001103.1555.002.    |
|                                             2001108115306.hdf"              |
+-----------------------------------------------------------------------------+
| PRODUCTIONDATETIME                          "1998-06-01T15:53:08.000Z"      |
+-----------------------------------------------------------------------------+
| DAYNIGHTFLAG                                "Day" or "Night" or "Both"      |
+-----------------------------------------------------------------------------+
| REPROCESSINGACTUAL                          "processed once"                |
+-----------------------------------------------------------------------------+
| REPROCESSINGPLANNED                         "no further update anticipated" |
+-----------------------------------------------------------------------------+
| SIZEMBECSDATAGRANULE                        (set by "DSS")                  |
+-----------------------------------------------------------------------------+
| AUTOMATICQUALITYFLAGEXPLANATION.1           "not being investigated"        |
+-----------------------------------------------------------------------------+
| AUTOMATICQUALITYFLAG.1                      "Suspect"                       |
+-----------------------------------------------------------------------------+
| OPERATIONALQUALITYFLAGEXPLANATION.1         (set by "DAAC")                 |
+-----------------------------------------------------------------------------+
| OPERATIONALQUALITYFLAG.1                    (set by "DAAC")                 |
+-----------------------------------------------------------------------------+
| SCIENCEQUALITYFLAGEXPLANATION.1             (set by "DP")                   |
+-----------------------------------------------------------------------------+
| SCIENCEQUALITYFLAG.1                        "Not Investigated"              |
+-----------------------------------------------------------------------------+
| QAPERCENTMISSINGDATA.1                      0                               |
+-----------------------------------------------------------------------------+
| QAPERCENTINTERPOLATEDDATA.1                 0                               |
+-----------------------------------------------------------------------------+
| QAPERCENTOUTOFBOUNDSDATA.1                  0                               |
+-----------------------------------------------------------------------------+
| PARAMETERNAME.1                             "EV_1KM_RefSB"                  |
+-----------------------------------------------------------------------------+
| AUTOMATICQUALITYFLAGEXPLANATION.2           "not being investigated"        |
+-----------------------------------------------------------------------------+
| AUTOMATICQUALITYFLAG.2                      "Suspect"                       |
+-----------------------------------------------------------------------------+
| OPERATIONALQUALITYFLAGEXPLANATION.2         (set by "DAAC")                 |
+-----------------------------------------------------------------------------+
| OPERATIONALQUALITYFLAG.2                    (set by "DAAC")                 |
+-----------------------------------------------------------------------------+
| SCIENCEQUALITYFLAGEXPLANATION.2             (set by "DP")                   |
+-----------------------------------------------------------------------------+
| SCIENCEQUALITYFLAG.2                        "Not Investigated"              |
+-----------------------------------------------------------------------------+
| QAPERCENTMISSINGDATA.2                      0                               |
+-----------------------------------------------------------------------------+
| QAPERCENTINTERPOLATEDDATA.2                 0                               |
+-----------------------------------------------------------------------------+
| QAPERCENTOUTOFBOUNDSDATA.2                  0                               |
+-----------------------------------------------------------------------------+
| PARAMETERNAME.2                             "EV_1KM_Emissive"               |
+-----------------------------------------------------------------------------+
| EQUATORCROSSINGDATE.1                       "1996-08-05"                    |
+-----------------------------------------------------------------------------+
| EQUATORCROSSINGTIME.1                       "15:55:45.854788"               |
+-----------------------------------------------------------------------------+
| ORBITNUMBER.1                               88                              |
+-----------------------------------------------------------------------------+
| EQUATORCROSSINGLONGITUDE.1                  -73.021282                      |
+-----------------------------------------------------------------------------+
| VERSIONID                                   1                               |
+-----------------------------------------------------------------------------+
| SHORTNAME                                   "MYD02OBC"   (Aqua)             |
+-----------------------------------------------------------------------------+
| INPUTPOINTER                                "L1A_v2_10scans_0.L1A",         |
|                                             "L1A_v2_10scans_1.L1A",         |
|                                             "L1A_v2_10scans_2.L1A",         |
|                                             "Reflective_Lookup_Tables_file",|
|                                             "Emissive_Lookup_Tables_file",  |
|                                             "QA_Lookup_Tables_file"         |
+-----------------------------------------------------------------------------+
| GRINGPOINTLONGITUDE.1                       (-86.567558, -59.389000,        |
|                                              -59.836601, -86.653564)        |
+-----------------------------------------------------------------------------+
| GRINGPOINTLATITUDE.1                        (41.644032, 37.729759,          |
|                                              36.824055, 40.686016)          |
+-----------------------------------------------------------------------------+
| GRINGPOINTSEQUENCENO.1                      (1, 2, 3, 4)                    |
+-----------------------------------------------------------------------------+
| EXCLUSIONGRINGFLAG.1                        "N"                             |
+-----------------------------------------------------------------------------+
| RANGEENDINGDATE                             "1996-08-05"                    |
+-----------------------------------------------------------------------------+
| RANGEENDINGTIME                             "15:56:00.626488"               |
+-----------------------------------------------------------------------------+
| RANGEBEGINNINGDATE                          "1996-08-05"                    |
+-----------------------------------------------------------------------------+
| RANGEBEGINNINGTIME                          "15:55:45.854788"               |
+-----------------------------------------------------------------------------+
| PGEVERSION                                  "2.3.0"                         |
+-----------------------------------------------------------------------------+
| ANCILLARYINPUTPOINTER.1                     "Geoloc_v2_10scans.GEO"         |
+-----------------------------------------------------------------------------+
| ANCILLARYINPUTTYPE.1                        "Geolocation"                   |
+-----------------------------------------------------------------------------+
| ADDITIONALATTRIBUTENAME.1                   "AveragedBlackBodyTemperature"  |
+-----------------------------------------------------------------------------+
| PARAMETERVALUE.1                            "268.90"                        |
+-----------------------------------------------------------------------------+
| ADDITIONALATTRIBUTENAME.2                   "AveragedMirrorTemperature"     |
+-----------------------------------------------------------------------------+
| PARAMETERVALUE.2                            "250.27"                        |
+-----------------------------------------------------------------------------+
| ADDITIONALATTRIBUTENAME.3                   "AveragedFocalPlane1Temperature"|
+-----------------------------------------------------------------------------+
| PARAMETERVALUE.3                            "361.42"                        |
+-----------------------------------------------------------------------------+
| ADDITIONALATTRIBUTENAME.4                   "AveragedFocalPlane2Temperature"|
+-----------------------------------------------------------------------------+
| PARAMETERVALUE.4                            "363.21"                        |
+-----------------------------------------------------------------------------+
| ADDITIONALATTRIBUTENAME.5                   "AveragedFocalPlane3Temperature"|
+-----------------------------------------------------------------------------+
| PARAMETERVALUE.5                            "118.25"                        |
+-----------------------------------------------------------------------------+
| ADDITIONALATTRIBUTENAME.6                   "AveragedFocalPlane4Temperature"|
+-----------------------------------------------------------------------------+
| PARAMETERVALUE.6                            "85.94"                         |
+-----------------------------------------------------------------------------+
| ADDITIONALATTRIBUTENAME.7                   "CalibrationQuality"            |
+-----------------------------------------------------------------------------+
| PARAMETERVALUE.7                            "marginal"                      |
+-----------------------------------------------------------------------------+
| ADDITIONALATTRIBUTENAME.8                   "MissionPhase"                  |
+-----------------------------------------------------------------------------+
| PARAMETERVALUE.8                            "A&E"                           |
+-----------------------------------------------------------------------------+
| ADDITIONALATTRIBUTENAME.9                   "NadirPointing"                 |
+-----------------------------------------------------------------------------+
| PARAMETERVALUE.9                            "Y"                             |
+-----------------------------------------------------------------------------+
| ASSOCIATEDPLATFORMSHORTNAME.1                "Aqua"                         |
+-----------------------------------------------------------------------------+
| ASSOCIATEDINSTRUMENTSHORTNAME.1             "MODIS"                         |
+-----------------------------------------------------------------------------+
| ASSOCIATEDSENSORSHORTNAME.1                 "MODIS"                         |
+=============================================================================+

+=============================================================================+
|                                                                             |
|  1.2) ECS Standard Archive Granule Metadata                                 |
|                                                                             |
|  Stored (as HDF ECS PVL string) in global attribute: "ArchiveMetadata.0"    |
|                                                                             |
+=============================================================================+
|                                                                             |
| Name                                        Example                         |
+-----------------------------------------------------------------------------+
| ALGORITHMPACKAGEACCEPTANCEDATE              "1999-11-01"                    |
+-----------------------------------------------------------------------------+
| ALGORITHMPACKAGEMATURITYCODE                "pre-launch"                    |
+-----------------------------------------------------------------------------+
| ALGORITHMPACKAGENAME                        "MODIS Level 1B Algorithm       |
|                                              Package"                       |
+-----------------------------------------------------------------------------+
| ALGORITHMPACKAGEVERSION                     "2.3.0.0"                       |
+-----------------------------------------------------------------------------+
| INSTRUMENTNAME                              "Moderate-Resolution Imaging    |
|                                              SpectroRadiometer"             |
+-----------------------------------------------------------------------------+
| PROCESSINGCENTER                            "GSFC"                          |
+-----------------------------------------------------------------------------+
| EASTBOUNDINGCOORDINATE                      40.000000                       |
+-----------------------------------------------------------------------------+
| WESTBOUNDINGCOORDINATE                      15.000000                       |
+-----------------------------------------------------------------------------+
| NORTHBOUNDINGCOORDINATE                     25.000000                       |
+-----------------------------------------------------------------------------+
| SOUTHBOUNDINGCOORDINATE                     10.000000                       |
+-----------------------------------------------------------------------------+
| DESCRREVISION                               "1.0"                           |
+-----------------------------------------------------------------------------+
| PRODUCTIONHISTORY                           "PGE02:3.4.5.6;PGE01:3.5.8"     |
+-----------------------------------------------------------------------------+
| LONGNAME                                    "MODIS/Aqua Calibrated          |
|                                             Radiances 5-Min L1B Swath 1km"  |
+-----------------------------------------------------------------------------+
| PROCESSINGENVIRONMENT                       "IRIX64"                        |
+=============================================================================+

+=============================================================================+
|                                                                             |
|  1.3) Product Granule Metadata Stored as Global Attributes                  |
|                                                                             |
+=============================================================================+
| Attr. Name:  "DN_obc_avg_first_frame_to_use"                                |
| Data Type:   INT16                                                          |
| Count:       1                                                              |
| Description: Number of frames to be used in computing the DN electronic     |
|              background numbers for reflective solar bands.                 |
+-----------------------------------------------------------------------------+
| Attr. Name:  "DN_obc_avg_number_of_frames_to_use"                           |
| Data Type:   INT16                                                          |
| Count:       1                                                              |
| Description: First frame of the set to be used in computing DN electronic   |
|              background numbers for reflective solar bands.                 |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Number of Scans"                                              |
| Data Type:   INT32                                                          |
| Count:       1                                                              |
| Description: The number of scans of data in the granule. This value is used |
|              to dimension many of the SDS and Vdata objects within this     |
|              file.  A dimension name of "nscans" or "num_scans" refers to   |
|              this value, which is typically 203 or 204.                     |
|              This value is copied from the input Level 1A middle granule.   |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Number of Day mode scans"                                     |
| Data Type:   INT32                                                          |
| Count:       1                                                              |
| Description: The subset of the Number of Scans that are classified by the   |
|              Level 1A code as day scans.  This value is copied from the     |
|              input Level 1A middle granule.                                 |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Number of Night mode scans"                                   |
| Data Type:   INT32                                                          |
| Count:       1                                                              |
| Description: The subset of the Number of Scans that are classified by the   |
|              Level 1A code as night scans.  This value is copied from the   |
|              input Level 1A middle granule.                                 |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Incomplete Scans"                                             |
| Data Type:   INT32                                                          |
| Count:       1                                                              |
| Description: The subset of the Number of Scans that are classified by the   |
|              Level 1A code as incomplete scans.  This value is copied from  |
|              the input Level 1A middle granule.                             |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Max Earth View Frames"                                        |
| Data Type:   INT32                                                          |
| Count:       1                                                              |
| Description: The number of 1km frames within the earth-view sector.  This   |
|              value is copied from the L1A attribute "Max Earth Frames"      |
|              of the input Level 1A middle granule.                          |
+-----------------------------------------------------------------------------+
| Attr. Name:  "%Valid EV Observations"                                       |
| Data Type:   FLOAT32                                                        |
| Count:       38                                                             |
| Description: The percentage of valid EV pixels in each band.  These values  |
|              are determined by subtracting from the total pixel count any   |
|              pixel that could not be calibrated (e.g., missing, saturated,  |
|              could not compute zero point, etc.).                           |
+-----------------------------------------------------------------------------+
| Attr. Name:  "%Saturated EV Observations"                                   |
| Data Type:   FLOAT32                                                        |
| Count:       38                                                             |
| Description: The percentage of EV pixels in each band that are determined   |
|              to be saturated.                                               |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Electronics Redundancy Vector"                                |
| Data Type:   UINT32                                                         |
| Count:       2                                                              |
| Description: Individual bits denote the status of one of the MODIS          |
|              electronics systems.  See the Level 1B EV file specifications  |
|              for a complete explanation.                                    |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Electronics Configuration Change"                             |
| Data Type:   UINT32                                                         |
| Count:       2                                                              |
| Description: Individual bits identify if any change occurred within the     |
|              granule for the corresponding telemetry point.  See the Level  |
|              1B EV file specifications for a complete explanation.          |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Reflective LUT Serial Number and Date of Last Change"         |
| Data Type:   CHAR8                                                          |
| Count:       30                                                             |
| Description: Serial number for the file "Reflective_Lookup_Tables_file",    |
|              which contains LUTs used principally in the reflective         |
|              solar bands calibration algorithm.                             |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Emissive LUT Serial Number and Date of Last Change"           |
| Data Type:   CHAR8                                                          |
| Count:       30                                                             |
| Description: Serial number for the file "Emissive_Lookup_Tables_file",      |
|              which contains LUTs used principally in the thermal            |
|              emissive bands calibration algorithm.                          |
+-----------------------------------------------------------------------------+
| Attr. Name:  "QA LUT Serial Number and Date of Last Change"                 |
| Data Type:   CHAR8                                                          |
| Count:       30                                                             |
| Description: Serial number for the file "QA_Lookup_Tables_file",            |
|              which contains LUTs used principally for QA analysis and       |
|              product identification.                                        |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Focal Plane Set Point State"                                  |
| Data Type:   INT8                                                           |
| Count:       1                                                              |
| Description: Identifies the temperature-control state of the focal planes.  |
|              0 = running uncontrolled,                                      |
|              1 = Set point is 83 degrees Celsius.                           |
|              2 = Set point is 85 degrees Celsius.                           |
|              3 = Set point is 88 degrees Celsius.                           |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Doors and Screens Configuration"                              |
| Data Type:   INT8                                                           |
| Count:       1                                                              |
| Description: Bit-wise fields denoting the state of the nadir-aperture door  |
|              (NAD), the space-view door (SVD), the solar diffuser door      |
|              (SDD) and the solar diffuser screen (SDSCRN) on the first scan |
|              of the granule.  The 4 most-significant bits of the word are   |
|              used in the order indicated above (NAD is the most-significant |
|              bit).  For the NAD, SVD and SDD, a value of 1 means open and 0 |
|              means closed.  For the SDSCRN, 1 means not screened and 0      |
|              means that the SD screen is in place.                          |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Reflective Bands With Bad Data"                               |
| Data Type:   INT8                                                           |
| Count:       22                                                             |
| Description: Flag indicating the presence (1) or absence (0) of bad data in |
|              each of the reflective solar bands.  The reflective solar      |
|              bands are MODIS bands 1-12,13lo,13hi,14lo,14hi,15-19, 26.      |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Emissive Bands With Bad Data"                                 |
| Data Type:   INT8                                                           |
| Count:       16                                                             |
| Description: Flag indicating the presence (1) or absence (0) of bad data in |
|              each of the thermal emissive bands.  The thermal emissive      |
|              bands are MODIS bands 20-25, 27-36.                            |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Noise in Black Body Thermistors"                              |
| Data Type:   UINT8                                                          |
| Count:       12                                                             |
| Description: A measure of the noise in each of the 12 black-body (BB)       |
|              thermistors.  The variance of the BB temperature over the      |
|              middle granule is divided by the pre-launch variance (from a   |
|              LUT).  Values of this ratio in the range of [0,9] are then     |
|              scaled to the range of [0,255] (values of 8-bit integer).      |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Noise in Average BB Temperature"                              |
| Data Type:   UINT8                                                          |
| Count:       1                                                              |
| Description: A measure of the noise of the average BB temperature.          |
|              (similar definition as in "Noise in Black Body Thermistors").  |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Noise in LWIR FPA Temperature"                                |
| Data Type:   UINT8                                                          |
| Count:       1                                                              |
| Description: A measure of the noise in the LWIR FPA temperature.            |
|              (similar definition as in "Noise in Black Body Thermistors").  |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Noise in MWIR FPA Temperature"                                |
| Data Type:   UINT8                                                          |
| Count:       1                                                              |
| Description: A measure of the noise in the MWIR FPA temperature.            |
|              (similar definition as in "Noise in Black Body Thermistors").  |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Noise in Scan Mirror Thermistor #1"                           |
| Data Type:   UINT8                                                          |
| Count:       1                                                              |
| Description: A measure of the noise in the Scan Mirror Thermistor #1.       |
|              (similar definition as in "Noise in Black Body Thermistors").  |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Noise in Scan Mirror Thermistor #2"                           |
| Data Type:   UINT8                                                          |
| Count:       1                                                              |
| Description: A measure of the noise in Scan Mirror Thermistor #2.           |
|              (similar definition as in "Noise in Black Body Thermistors").  |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Noise in Scan Mirror Thermistor Average"                      |
| Data Type:   UINT8                                                          |
| Count:       1                                                              |
| Description: A measure of the noise in the average SM temperature.          |
|              (similar definition as in "Noise in Black Body Thermistors").  |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Noise in Instrument Temperature"                              |
| Data Type:   UINT8                                                          |
| Count:       1                                                              |
| Description: A measure of the noise in the instrument temperature.          |
|              (similar definition as in "Noise in Black Body Thermistors").  |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Noise in Cavity Temperature"                                  |
| Data Type:   UINT8                                                          |
| Count:       1                                                              |
| Description: A measure of the noise in the cavity temperature.              |
|              (similar definition as in "Noise in Black Body Thermistors").  |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Noise in Temperature of NIR FPA"                              |
| Data Type:   UINT8                                                          |
| Count:       1                                                              |
| Description: A measure of the noise in the NIR FPA temperature.             |
|              (similar definition as in "Noise in Black Body Thermistors").  |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Noise in Temperature of Vis FPA"                              |
| Data Type:   UINT8                                                          |
| Count:       1                                                              |
| Description: A measure of the noise in the VIS FPA temperature.             |
|              (similar definition as in "Noise in Black Body Thermistors").  |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Dead Detector List"                                           |
| Data Type:   INT8                                                           |
| Count:       490                                                            |
| Description: List of detectors identifying those which do not provide data  |
|              of usable quality (= 1). Usable detectors are flagged as 0.    |
|              EV pixels are not calibrated if the associated flag for dead   |
|              detector is set to 1.  Values of this come from the atttribute |
|              "Detector Quality Flag", described below.                      |
|              The order of detectors follows the order of the 38 MODIS       |
|              bands, with the appropriate number of detectors in each band,  |
|              and with the "product" detector order convention.              |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Noisy Detector List"                                          |
| Data Type:   INT8                                                           |
| Count:       490                                                            |
| Description: List of detectors flagged as noisy, but otherwise produce      |
|              usable data.  These values come from the atttribute            |
|              "Detector Quality Flag", described below.                      |
|              The ordering is similar to "Dead Detector List", above.        |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Dead Subframe List"                                           |
| Data Type:   INT8                                                           |
| Count:       520                                                            |
| Description: List of subframes identifying those which do not provide data  |
|              of usable quality (= 1). Usable subframes are flagged as 0.    |
|              EV pixels are not calibrated if the associated flag for dead   |
|              subframe is set to 1.  Values of this come from the atttribute |
|              "Detector Quality Flag2", described below.                     |
|              The order of subframes follows the order of the 7 MODIS 250m   |
|              and 500m bands, with the appropriate number of detectors in    |
|              each band, and with the "product" detector order convention.   |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Noisy Subframe List"                                          |
| Data Type:   INT8                                                           |
| Count:       520                                                            |
| Description: List of subframes flagged as noisy, but otherwise produce      |
|              usable data.  These values come from the atttribute            |
|              "Detector Quality Flag2", described below.                     |
|              The ordering is similar to "Dead Subframe List", above.        |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Detector Quality Flag"                                        |
| Data Type:   UINT8                                                          |
| Count:       490                                                            |
| Description: Individual bits identify noisy, dead and other anomalous       |
|              detector quality conditions (see EV file specifications)       |
|              Values come from a LUT.                                        |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Detector Quality Flag2"                                       |
| Data Type:   UINT8                                                          |
| Count:       180                                                            |
| Description: Individual bits identify noisy and dead subframes.             |
|              (see EV file specifications)                                   |
|              Values come from a LUT.                                        |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Earth-Sun Distance"                                           |
| Data Type:   FLOAT32                                                        |
| Count:       1                                                              |
| Description: Identifies the earth-sun distance, relative to 1 AU, for the   |
|              granule as a whole.                                            |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Solar Irradiance on RSB Detectors over pi"                    |
| Data Type:   FLOAT32                                                        |
| Count:       330                                                            |
| Description: Contains the RSR-weighted solar irradiance at 1 AU divided by  |
|              pi for each reflective detector.                               |
+-----------------------------------------------------------------------------+
| Attr. Name:  "% L1A EV All Scan Data are Missing"                           |
| Data Type:   FLOAT32                                                        |
| Count:       1                                                              |
| Description: Identifies the percentage of completely missing scans, defined |
|              as having no packets (see the L1A SDS Scan quality array)      |
+-----------------------------------------------------------------------------+
| Attr. Name:  "% L1A EV RSB DN Not in Day Mode"                              |
| Data Type:   FLOAT32                                                        |
| Count:       490                                                            |
| Description: Identifies, for each detector, the percentage of pixels that   |
|              are set to an unusable data value because MODIS was not in     |
|              day mode.                                                      |
+-----------------------------------------------------------------------------+
| Attr. Name:  "% L1A EV DN Missing Within Scan"                              |
| Data Type:   FLOAT32                                                        |
| Count:       490                                                            |
| Description: Identifies, for each detector, the percentage of pixels that   |
|              have DN missing within a scan (that is not a completely        |
|              missing scan and MODIS was in day mode).                       |
+-----------------------------------------------------------------------------+
| Attr. Name:  "% Dead Detector EV Data"                                      |
| Data Type:   FLOAT32                                                        |
| Count:       490                                                            |
| Description: Identifies, for each detector, the percentage of pixels that   |
|              have scaled integers set to the dead detector unusable data    |
|              value.                                                         |
+-----------------------------------------------------------------------------+
| Attr. Name:  "% Dead Subframe EV Data"                                      |
| Data Type:   FLOAT32                                                        |
| Count:       180                                                            |
| Description: Identifies, for each detector in 250m and 500m bands, the      |
|              percentage of pixels that have scaled integers set to the dead |
|              subframe unusable data value.                                  |
+-----------------------------------------------------------------------------+
| Attr. Name:  "% Sector Rotation EV Data"                                    |
| Data Type:   FLOAT32                                                        |
| Count:       490                                                            |
| Description: Identifies, for each detector, the percentage of pixels that   |
|              have scaled integers set to the sector rotation unusable data  |
|              value.                                                         |
+-----------------------------------------------------------------------------+
| Attr. Name:  "% Saturated EV Data"                                          |
| Data Type:   FLOAT32                                                        |
| Count:       490                                                            |
| Description: Identifies, for each detector, the percentage of pixels that   |
|              have scaled integers set to the saturated unusable data value. |
+-----------------------------------------------------------------------------+
| Attr. Name:  "% TEB EV Data With Moon in SVP"                               |
| Data Type:   FLOAT32                                                        |
| Count:       490                                                            |
| Description: Identifies, for each detector, the percentage of pixels that   |
|              have scaled integers set to the moon in SVP unusable data      |
|              value. (RSB detectors have this value set to zero)             |
+-----------------------------------------------------------------------------+
| Attr. Name:  "% EV Data Where Cannot Compute BG DN"                         |
| Data Type:   FLOAT32                                                        |
| Count:       490                                                            |
| Description: Identifies, for each detector, the percentage of pixels that   |
|              have scaled integers set to the cannot complute zero point     |
|              DN unusable data value.                                        |
+-----------------------------------------------------------------------------+
| Attr. Name:  "% RSB EV Data With dn** Below Scale"                          |
| Data Type:   FLOAT32                                                        |
| Count:       490                                                            |
| Description: Identifies, for each detector, the percentage of pixels that   |
|              have scaled integers set to the dn** below scale unusable data |
|              value. (TEB detectors have this value set to zero)             |
+-----------------------------------------------------------------------------+
| Attr. Name:  "% EV Data Where Nadir Door Closed"                            |
| Data Type:   FLOAT32                                                        |
| Count:       490                                                            |
| Description: Identifies, for each detector, the percentage of pixels that   |
|              have scaled integers set to the NAD closed unusable data value.|
+-----------------------------------------------------------------------------+
| Attr. Name:  "% EV Data Not Calibrated"                                     |
| Data Type:   FLOAT32                                                        |
| Count:       490                                                            |
| Description: Identifies, for each detector, the percentage of pixels that   |
|              have scaled integers set to an unusable data value.  (This is  |
|              a roll-up of preceding attributes)                             |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Bit QA Flags Last Value"                                      |
| Data Type:   UINT32                                                         |
| Count:       1                                                              |
| Description: The value of the Bit QA Flags on the last scan in the granule. |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Bit QA Flags Change"                                          |
| Data Type:   UINT32                                                         |
| Count:       1                                                              |
| Description: Bit-by-bit flag identifying that a change in Bit QA Flags      |
|              occurred within granule.                                       |
+-----------------------------------------------------------------------------+
| Attr. Name:  "Granule Average QA Values"                                    |
| Data Type:   FLOAT32                                                        |
| Count:       50                                                             |
| Description: Granule average temperatures and voltages:                     |
|              TP_BB_TEMP01          BB Thermister 1                          |
|              TP_BB_TEMP02          BB Thermister 2                          |
|              TP_BB_TEMP03          BB Thermister 3                          |
|              TP_BB_TEMP04          BB Thermister 4                          |
|              TP_BB_TEMP05          BB Thermister 5                          |
|              TP_BB_TEMP06          BB Thermister 6                          |
|              TP_BB_TEMP07          BB Thermister 7                          |
|              TP_BB_TEMP08          BB Thermister 8                          |
|              TP_BB_TEMP09          BB Thermister 9                          |
|              TP_BB_TEMP10          BB Thermister 10                         |
|              TP_BB_TEMP11          BB Thermister 11                         |
|              TP_BB_TEMP12          BB Thermister 12                         |
|              TA_AO_VIS_FPAE        (FPA1) VIS Focal Plane                   |
|              TA_AO_NIR_FPAE        (FPA2) NIR Focal Plane                   |
|              TA_RC_SMIR_CFPAE      (FPA3) SMIR Focal plane Temp             |
|              TA_RC_LWIR_CFPAE      (FPA4) LWIR Focal plane Temp             |
|              TP_SA_RCT1_MIRE       Scan Mirror Temp 1                       |
|              TP_SA_RCT2_MIRE       Scan Mirror Temp 2                       |
|              TP_SA_A_MTR           Scan motor temperature A                 |
|              TP_MF_CALBKHD_SR      (CAV1) Cal bulkhead below SRCA mount     |
|              TP_SR_SNOUT           (CAV2) SRCA Bulkhead Temp                |
|              TP_MF_Z_BKHD_BB       (CAV3) Mid zenith bulkhead near BB       |
|              TP_MF_CVR_OP_SR       (CAV4) cover opposite SRCA               |
|              TP_AO_SMIR_OBJ        (INS1) SMIR Objective Lens Temp          |
|              TP_AO_LWIR_OBJ        (INS2) LWIR Objective Lens Temp          |
|              TP_AO_SMIR_LENS       (INS3) SMIR Eye Assy Temp                |
|              TP_AO_LWIR_LENS       (INS4) LWIR Eye Assy Temp                |
|              TA_RC_CS              (RC1) Cold stage Temp                    |
|              TA_RC_CS_OG           (RC2) Cold stage outgas Temp             |
|              TA_RC_IS              (RC3) Intermediate Stage Temp            |
|              TA_RC_IS_OG           (RC4) Intermediate Stage outgas Temp     |
|              TA_RC_OS_OG           (RC5) Outerstage outgas Temp             |
|              VR_RC_LW_FPA_HTR      LWIR heater voltage                      |
|              (elements 34-50 are for future use)                            |
+=============================================================================+

II)  Scientific Data Sets (SDSs)

This section contains HDF Scientific Data Sets (SDSs).  Many of these are
copied directly from the middle L1A input granule.  The order of dimensions
for a multidimensional array follows the C language convention:

      [least-rapidly varying, ... , most-rapidly varying]

Acronyms used throughout this section:

    SD    Solar Diffuser
    SRCA  Spectro-Radiometric Calibration Assembly
    SV    Space View
    BB    Black-body
    EV    Earth-view
    DN    Digital number

Dimension names for SDSs:

     name           value                  meaning
 ---------------   ---------    ---------------------------------------
 "nscans"          (variable)    value of attribute "Number of Scans"
 "vecdim"               3        number of components in cartesian vector
 "num_bands"           38        number of MODIS bands
 "Band_250m"            2        number of 250m resolution bands
 "Band_500m"            5        number of 500m resolution bands
 "Band_1km_day"        14        number of 1km_day resolution bands
 "Band_1km_night"      17        number of 1km_night resolution bands
 "250m_subsamples"      4        number of sub-samples per 1km frame for
                                 250m-resolution bands
 "500m_subsamples"      2        number of sub-samples per 1km frame for
                                 500m-resolution bands
 "1km_subsamples"       1        number of sub-samples per 1km frame for
                                 1km-resolution bands
 "SD_frames"           50        number of 1km frames in SD sector
 "SRCA_frames"         10        number of 1km frames in SRCA sector
 "SV_frames"           50        number of 1km frames in SV sector
 "BB_frames"           50        number of 1km frames in BB sector

 Note:  If a number appears as a multiplicative factor in a dimension name,
        The value of the dimension is the product indicated.  For example,
        the dimension "2*SD_frames" has the value 100.  Specifically, a
        number multiplying "nscans" is the number of detectors for that
        resolution and a number multiplying a "frames" dimension indicates
        the number of subsamples (or subframes) for that resolution.
        Additionally, some dimension names consist only of a number,
        e.g. "550" also has the value 550.

+=============================================================================+
| SDS Name:    "SD_250m"                                                      |
| Data Type:   INT16                                                          |
| Rank:        3                                                              |
| Dimensions:  ["40*nscans","Band_250m","4*SD_frames"]                        |
| Description: Raw, uncorrected, DNs for the SD sector for 250m bands.        |
|              These are copied directly from the middle L1A input granule.   |
+-----------------------------------------------------------------------------+
| SDS Name:    "SD_500m"                                                      |
| Data Type:   INT16                                                          |
| Rank:        3                                                              |
| Dimensions:  ["20*nscans","Band_500m","2*SD_frames"]                        |
| Description: Raw, uncorrected, DNs for the SD sector for 500m bands.        |
|              These are copied directly from the middle L1A input granule.   |
+-----------------------------------------------------------------------------+
| SDS Name:    "SD_1km_day"                                                   |
| Data Type:   INT16                                                          |
| Rank:        3                                                              |
| Dimensions:  ["10*nscans","Band_1km_day","SD_frames"]                       |
| Description: Raw, uncorrected, DNs for the SD sector for 1km_day bands.     |
|              These are copied directly from the middle L1A input granule.   |
+-----------------------------------------------------------------------------+
| SDS Name:    "SD_1km_night"                                                 |
| Data Type:   INT16                                                          |
| Rank:        3                                                              |
| Dimensions:  ["10*nscans","Band_1km_night","SD_frames"]                     |
| Description: Raw, uncorrected, DNs for the SD sector for 1km_night bands.   |
|              These are copied directly from the middle L1A input granule.   |
+-----------------------------------------------------------------------------+
| SDS Name:    "SRCA_250m"                                                    |
| Data Type:   INT16                                                          |
| Rank:        3                                                              |
| Dimensions:  ["40*nscans","Band_250m","4*SRCA_frames"]                      |
| Description: Raw, uncorrected, DNs for the SRCA sector for 250m bands.      |
|              These are copied directly from the middle L1A input granule.   |
+-----------------------------------------------------------------------------+
| SDS Name:    "SRCA_500m"                                                    |
| Data Type:   INT16                                                          |
| Rank:        3                                                              |
| Dimensions:  ["20*nscans","Band_500m","2*SRCA_frames"]                      |
| Description: Raw, uncorrected, DNs for the SRCA sector for 500m bands.      |
|              These are copied directly from the middle L1A input granule.   |
+-----------------------------------------------------------------------------+
| SDS Name:    "SRCA_1km_day"                                                 |
| Data Type:   INT16                                                          |
| Rank:        3                                                              |
| Dimensions:  ["10*nscans","Band_1km_day","SRCA_frames"]                     |
| Description: Raw, uncorrected, DNs for the SRCA sector for 1km_day bands.   |
|              These are copied directly from the middle L1A input granule.   |
+-----------------------------------------------------------------------------+
| SDS Name:    "SRCA_1km_night"                                               |
| Data Type:   INT16                                                          |
| Rank:        3                                                              |
| Dimensions:  ["10*nscans","Band_1km_night","SRCA_frames"]                   |
| Description: Raw, uncorrected, DNs for the SRCA sector for 1km_night bands. |
|              These are copied directly from the middle L1A input granule.   |
+-----------------------------------------------------------------------------+
| SDS Name:    "BB_250m"                                                      |
| Data Type:   INT16                                                          |
| Rank:        3                                                              |
| Dimensions:  ["40*nscans","Band_250m","4*BB_frames"]                        |
| Description: Raw, uncorrected, DNs for the BB sector for 250m bands.        |
|              These are copied directly from the middle L1A input granule.   |
+-----------------------------------------------------------------------------+
| SDS Name:    "BB_500m"                                                      |
| Data Type:   INT16                                                          |
| Rank:        3                                                              |
| Dimensions:  ["20*nscans","Band_500m","2*BB_frames"]                        |
| Description: Raw, uncorrected, DNs for the BB sector for 500m bands.        |
|              These are copied directly from the middle L1A input granule.   |
+-----------------------------------------------------------------------------+
| SDS Name:    "BB_1km_day"                                                   |
| Data Type:   INT16                                                          |
| Rank:        3                                                              |
| Dimensions:  ["10*nscans","Band_1km_day","BB_frames"]                       |
| Description: Raw, uncorrected, DNs for the BB sector for 1km_day bands.     |
|              These are copied directly from the middle L1A input granule.   |
+-----------------------------------------------------------------------------+
| SDS Name:    "BB_1km_night"                                                 |
| Data Type:   INT16                                                          |
| Rank:        3                                                              |
| Dimensions:  ["10*nscans","Band_1km_night","BB_frames"]                     |
| Description: Raw, uncorrected, DNs for the BB sector for 1km_night bands.   |
|              These are copied directly from the middle L1A input granule.   |
+-----------------------------------------------------------------------------+
| SDS Name:    "SV_250m"                                                      |
| Data Type:   INT16                                                          |
| Rank:        3                                                              |
| Dimensions:  ["40*nscans","Band_250m","4*SV_frames"]                        |
| Description: Raw, uncorrected, DNs for the SV sector for 250m bands.        |
|              These are copied directly from the middle L1A input granule.   |
+-----------------------------------------------------------------------------+
| SDS Name:    "SV_500m"                                                      |
| Data Type:   INT16                                                          |
| Rank:        3                                                              |
| Dimensions:  ["20*nscans","Band_500m","2*SV_frames"]                        |
| Description: Raw, uncorrected, DNs for the SV sector for 500m bands.        |
|              These are copied directly from the middle L1A input granule.   |
+-----------------------------------------------------------------------------+
| SDS Name:    "SV_1km_day"                                                   |
| Data Type:   INT16                                                          |
| Rank:        3                                                              |
| Dimensions:  ["10*nscans","Band_1km_day","SV_frames"]                       |
| Description: Raw, uncorrected, DNs for the SV sector for 1km_day bands.     |
|              These are copied directly from the middle L1A input granule.   |
+-----------------------------------------------------------------------------+
| SDS Name:    "SV_1km_night"                                                 |
| Data Type:   INT16                                                          |
| Rank:        3                                                              |
| Dimensions:  ["10*nscans","Band_1km_night","SV_frames"]                     |
| Description: Raw, uncorrected, DNs for the SV sector for 1km_night bands.   |
|              These are copied directly from the middle L1A input granule.   |
+-----------------------------------------------------------------------------+
| SDS Name:    "Moon Vector"                                                  |
| Data Type:   FLOAT32                                                        |
| Rank:        2                                                              |
| Dimensions:  ["nscans","vecdim"]                                            |
| Description: The (x,y,z) components of a unit vector that points toward     |
|              the moon within the MODIS instrument coordinate system.        |
|              These values are copied from the MOD03 geolocation file.       |
+-----------------------------------------------------------------------------+
| SDS Name:    "Moon in keep-out-box"                                         |
| Data Type:   INT8                                                           |
| Rank:        2                                                              |
| Dimensions:  ["nscans","num_bands"]                                         |
| Description: Flags that are calculated within the Level 1B code that        |
|              identify if the moon was within the SV keep-out box (KOB)      |
|              for a given band and scan (0: not inside, 1: inside).  The     |
|              limits of the KOB come from a Level 1B LUT and the moon        |
|              vector (described above) is used to make this determination.   |
+-----------------------------------------------------------------------------+
| SDS Name:    "SD Sun Azimuth"                                               |
| Data Type:   FLOAT32                                                        |
| Rank:        1                                                              |
| Dimensions:  ["nscans"]                                                     |
| Description: Sun vector azimuth angle in the SD frame of reference,         |
|              measured in a clock-wise rotation from the SD Y axis about     |
|              the SD Z axis.  These values are copied from the               |
|              MOD03 geolocation file (same units as that file).              |
+-----------------------------------------------------------------------------+
| SDS Name:    "SD Sun Zenith"                                                |
| Data Type:   FLOAT32                                                        |
| Rank:        1                                                              |
| Dimensions:  ["nscans"]                                                     |
| Description: Sun vector zenith angle in the SD frame of reference,          |
|              measured from the SD Z axis.  These values are copied from the |
|              MOD03 geolocation file (same units as that file).              |
+-----------------------------------------------------------------------------+
| SDS Name:    "DN_obc_avg_250m"                                              |
| Data Type:   FLOAT32                                                        |
| Rank:        3                                                              |
| Dimensions:  ["40*nscans","Band_250m","250m_subsamples"]                    |
| Description: Average DN to use as the zero point DN value in calibration    |
|              algorithms for 250m bands (bands 1,2).  These values are       |
|              calculated within Level 1B using the algorithm described in    |
|              Godden, et. al., "A common algorithm for handling the          |
|              instrument and electronic background for the reflective and    |
|              thermal bands", MCST document #M0656.                          |
+-----------------------------------------------------------------------------+
| SDS Name:    "DN_obc_var_250m"                                              |
| Data Type:   FLOAT32                                                        |
| Rank:        3                                                              |
| Dimensions:  ["40*nscans","Band_250m","250m_subsamples"]                    |
| Description: Variance of values used to compute DN_obc_avg_250m.  The       |
|              value computed does not include any outliers that may have     |
|              been rejected using the outlier rejection algorithem described |
|              in the reference noted above.                                  |
+-----------------------------------------------------------------------------+
| SDS Name:    "DN_obc_outlier_mask_250m"                                     |
| Data Type:   UINT32                                                         |
| Rank:        4                                                              |
| Dimensions:  ["40*nscans","Band_250m","250m_subsamples","2"]                |
| Description: Bit mask which identifies (with a 1) those pixels rejected by  |
|              the outlier algorithm when computing DN_obc_avg_250m.  The     |
|              most rapidly-varying dimension is used to provide a total of   |
|              64 bits, of which the first 50 are used.  The 50 values map    |
|              to the 50 frames available for providing data for each         |
|              average.  The 1st frame corresponds to the least significant   |
|              bit of the 1st word of the 2-word array.                       |
+-----------------------------------------------------------------------------+
| SDS Name:    "DN_obc_avg_500m"                                              |
| Data Type:   FLOAT32                                                        |
| Rank:        3                                                              |
| Dimensions:  ["20*nscans","Band_500m","500m_subsamples"]                    |
| Description: Average DN to use as the zero point DN value in calibration    |
|              algorithms for 500m bands.                                     |
+-----------------------------------------------------------------------------+
| SDS Name:    "DN_obc_var_500m"                                              |
| Data Type:   FLOAT32                                                        |
| Rank:        3                                                              |
| Dimensions:  ["20*nscans","Band_500m","500m_subsamples"]                    |
| Description: Variance of values used to compute DN_obc_avg_500m.            |
+-----------------------------------------------------------------------------+
| SDS Name:    "DN_obc_outlier_mask_500m"                                     |
| Data Type:   UINT32                                                         |
| Rank:        4                                                              |
| Dimensions:  ["20*nscans","Band_500m","500m_subsamples","2"]                |
| Description: Outlier mask for values used to compute DN_obc_avg_500m.       |
+-----------------------------------------------------------------------------+
| SDS Name:    "DN_obc_avg_1km_day"                                           |
| Data Type:   FLOAT32                                                        |
| Rank:        3                                                              |
| Dimensions:  ["10*nscans","Band_1km_day","1km_subsamples"]                  |
| Description: Average DN to use as the zero point DN value in calibration    |
|              algorithms for 1km_day bands.                                  |
+-----------------------------------------------------------------------------+
| SDS Name:    "DN_obc_var_1km_day"                                           |
| Data Type:   FLOAT32                                                        |
| Rank:        3                                                              |
| Dimensions:  ["10*nscans","Band_1km_day","1km_subsamples"]                  |
| Description: Variance of values used to compute DN_obc_avg_1km_day.         |
+-----------------------------------------------------------------------------+
| SDS Name:    "DN_obc_outlier_mask_1km_day"                                  |
| Data Type:   UINT32                                                         |
| Rank:        4                                                              |
| Dimensions:  ["10*nscans","Band_1km_day","1km_subsamples","2"]              |
| Description: Outlier mask for values used to compute DN_obc_avg_1km_day.    |
+-----------------------------------------------------------------------------+
| SDS Name:    "DN_obc_avg_1km_night"                                         |
| Data Type:   FLOAT32                                                        |
| Rank:        3                                                              |
| Dimensions:  ["10*nscans","Band_1km_night","1km_subsamples"]                |
| Description: Average DN to use as the zero point DN value in calibration    |
|              algorithms for 1km_night bands.                                |
+-----------------------------------------------------------------------------+
| SDS Name:    "DN_obc_var_1km_night"                                         |
| Data Type:   FLOAT32                                                        |
| Rank:        3                                                              |
| Dimensions:  ["10*nscans","Band_1km_night","1km_subsamples"]                |
| Description: Variance of values used to compute DN_obc_var_1km_night.       |
+-----------------------------------------------------------------------------+
| SDS Name:    "DN_obc_outlier_mask_1km_night"                                |
| Data Type:   UINT32                                                         |
| Rank:        4                                                              |
| Dimensions:  ["10*nscans","Band_1km_night","1km_subsamples","2"]            |
| Description: Outlier mask for values used to compute DN_obc_var_1km_night.  |
+-----------------------------------------------------------------------------+
| SDS Name:    "fpa_aem_config"                                               |
| Data Type:   INT8                                                           |
| Rank:        2                                                              |
| Dimensions:  ["nscans","10"]                                                |
| Description: See the MOD01 (Level 1A) file spec.  This SDS is copied        |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "science_state"                                                |
| Data Type:   INT8                                                           |
| Rank:        1                                                              |
| Dimensions:  ["nscans"]                                                     |
| Description: See the MOD01 (Level 1A) file spec.  This SDS is copied        |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "science_abnormal"                                             |
| Data Type:   INT8                                                           |
| Rank:        1                                                              |
| Dimensions:  ["nscans"]                                                     |
| Description: See the MOD01 (Level 1A) file spec.  This SDS is copied        |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "fpa_dcr_offset"                                               |
| Data Type:   INT8                                                           |
| Rank:        2                                                              |
| Dimensions:  ["nscans","550"]                                               |
| Description: See the MOD01 (Level 1A) file spec.  This SDS is copied        |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "raw_mir_enc"                                                  |
| Data Type:   INT16                                                          |
| Rank:        2                                                              |
| Dimensions:  ["nscans","78"]                                                |
| Description: See the MOD01 (Level 1A) file spec.  This SDS is copied        |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "raw_vs_def"                                                   |
| Data Type:   INT16                                                          |
| Rank:        2                                                              |
| Dimensions:  ["nscans","40"]                                                |
| Description: See the MOD01 (Level 1A) file spec.  This SDS is copied        |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "raw_vs_act"                                                   |
| Data Type:   INT16                                                          |
| Rank:        2                                                              |
| Dimensions:  ["nscans","24"]                                                |
| Description: See the MOD01 (Level 1A) file spec.  This SDS is copied        |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "raw_sci_eng"                                                  |
| Data Type:   INT8                                                           |
| Rank:        2                                                              |
| Dimensions:  ["nscans","212"]                                               |
| Description: See the MOD01 (Level 1A) file spec.  This SDS is copied        |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "raw_hk_telem"                                                 |
| Data Type:   INT8                                                           |
| Rank:        2                                                              |
| Dimensions:  ["nscans","128"]                                               |
| Description: See the MOD01 (Level 1A) file spec.  This SDS is copied        |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "raw_sc_ancil"                                                 |
| Data Type:   INT16                                                          |
| Rank:        2                                                              |
| Dimensions:  ["nscans","64"]                                                |
| Description: See the MOD01 (Level 1A) file spec.  This SDS is copied        |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "raw_param"                                                    |
| Data Type:   INT8                                                           |
| Rank:        2                                                              |
| Dimensions:  ["nscans","40"]                                                |
| Description: See the MOD01 (Level 1A) file spec.  This SDS is copied        |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "raw_pv_gains"                                                 |
| Data Type:   INT8                                                           |
| Rank:        2                                                              |
| Dimensions:  ["nscans","550"]                                               |
| Description: See the MOD01 (Level 1A) file spec.  This SDS is copied        |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "Scan number"                                                  |
| Data Type:   INT16                                                          |
| Rank:        1                                                              |
| Dimensions:  ["nscans"]                                                     |
| Description: See the MOD01 (Level 1A) file spec.  This SDS is copied        |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "Frame count array"                                            |
| Data Type:   INT16                                                          |
| Rank:        2                                                              |
| Dimensions:  ["nscans","6"]                                                 |
| Description: See the MOD01 (Level 1A) file spec.  This SDS is copied        |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "Scan Type"                                                    |
| Data Type:   CHAR8                                                          |
| Rank:        2                                                              |
| Dimensions:  ["nscans","10"]                                                |
| Description: See the MOD01 (Level 1A) file spec.  This SDS is copied        |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "SD start time"                                                |
| Data Type:   FLOAT64                                                        |
| Rank:        1                                                              |
| Dimensions:  ["nscans"]                                                     |
| Description: TAI time (sec.) starting the SD sector.  This SDS is copied    |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "SRCA start time"                                              |
| Data Type:   FLOAT64                                                        |
| Rank:        1                                                              |
| Dimensions:  ["nscans"]                                                     |
| Description: TAI time (sec.) starting the SRCA sector.  This SDS is copied  |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "BB start time"                                                |
| Data Type:   FLOAT64                                                        |
| Rank:        1                                                              |
| Dimensions:  ["nscans"]                                                     |
| Description: TAI time (sec.) starting the BB sector.  This SDS is copied    |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "SV start time"                                                |
| Data Type:   FLOAT64                                                        |
| Rank:        1                                                              |
| Dimensions:  ["nscans"]                                                     |
| Description: TAI time (sec.) starting the SV sector.  This SDS is copied    |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "EV start time"                                                |
| Data Type:   FLOAT64                                                        |
| Rank:        1                                                              |
| Dimensions:  ["nscans"]                                                     |
| Description: TAI time (sec.) starting the EV sector.  This SDS is copied    |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "SRCA calibration mode"                                        |
| Data Type:   INT16                                                          |
| Rank:        1                                                              |
| Dimensions:  ["nscans"]                                                     |
| Description: 0 = radiometric, 1 = spatial, 2 = spectral, 3 = off.           |
|              This SDS is copied directly from the middle L1A input granule. |
+-----------------------------------------------------------------------------+
| SDS Name:    "Packet scan count"                                            |
| Data Type:   INT16                                                          |
| Rank:        1                                                              |
| Dimensions:  ["nscans"]                                                     |
| Description: See the MOD01 (Level 1A) file spec.  This SDS is copied        |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "CCSDS Application Identifiers"                                |
| Data Type:   INT16                                                          |
| Rank:        2                                                              |
| Dimensions:  ["nscans","3"]                                                 |
| Description: See the MOD01 (Level 1A) file spec.  This SDS is copied        |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "Packet expedited data flag"                                   |
| Data Type:   INT16                                                          |
| Rank:        1                                                              |
| Dimensions:  ["nscans"]                                                     |
| Description: See the MOD01 (Level 1A) file spec.  This SDS is copied        |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "Mirror side"                                                  |
| Data Type:   INT16                                                          |
| Rank:        1                                                              |
| Dimensions:  ["nscans"]                                                     |
| Description: Mirror side index (0 or 1).  This SDS is copied                |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "Scan quality array"                                           |
| Data Type:   INT32                                                          |
| Rank:        2                                                              |
| Dimensions:  ["nscans","4"]                                                 |
| Description: See the MOD01 (Level 1A) file spec.  This SDS is copied        |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "SD sector Pixel quality"                                      |
| Data Type:   INT16                                                          |
| Rank:        3                                                              |
| Dimensions:  ["nscans","64","2"]                                            |
| Description: See the MOD01 (Level 1A) file spec.  This SDS is copied        |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "SRCA sector Pixel quality"                                    |
| Data Type:   INT16                                                          |
| Rank:        3                                                              |
| Dimensions:  ["nscans","64","2"]                                            |
| Description: See the MOD01 (Level 1A) file spec.  This SDS is copied        |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "BB sector Pixel quality"                                      |
| Data Type:   INT16                                                          |
| Rank:        3                                                              |
| Dimensions:  ["nscans","64","2"]                                            |
| Description: See the MOD01 (Level 1A) file spec.  This SDS is copied        |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "SV sector Pixel quality"                                      |
| Data Type:   INT16                                                          |
| Rank:        3                                                              |
| Dimensions:  ["nscans","64","2"]                                            |
| Description: See the MOD01 (Level 1A) file spec.  This SDS is copied        |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
| SDS Name:    "Earth sector Pixel quality"                                   |
| Data Type:   INT16                                                          |
| Rank:        3                                                              |
| Dimensions:  ["nscans","1400","2"]                                          |
| Description: See the MOD01 (Level 1A) file spec.  This SDS is copied        |
|              directly from the middle L1A input granule.                    |
+-----------------------------------------------------------------------------+
+-----------------------------------------------------------------------------+
| SDS Name:    "Bit QA Flags"                                                 |
| Data Type:   UINT32                                                         |
| Rank:        1                                                              |
| Dimensions:  ["nscans"]                                                     |
| Description: Same meanings as in the Level 1B EV products. (In those        |
|              products, the Bit QA Flags are one field of the "Level 1B      |
|              Swath Metadata" Vdata.  Here, the values are placed in an SDS.)|
|                                                                             |
|               Bit QA Flags are as follows (bit 0 is the least significant   |
|               bit in the word):                                             |
|               +========+==========================================+         |
|               | bit #  |   flag name/desciption (1=True, 0=False) |         |
|               +========+==========================================+         |
|               | bit 0  |   Moon within defined limits of SVP      |         |
|               +--------+------------------------------------------+         |
|               | bit 1  |   Spacecraft Maneuver                    |         |
|               +--------+------------------------------------------+         |
|               | bit 2  |   Sector Rotation                        |         |
|               +--------+------------------------------------------+         |
|               | bit 3  |   Negative Radiance Beyond Noise Level   |         |
|               +--------+------------------------------------------+         |
|               | bit 4  |   PC Ecal on                             |         |
|               +--------+------------------------------------------+         |
|               | bit 5  |   PV Ecal on                             |         |
|               +--------+------------------------------------------+         |
|               | bit 6  |   SD Door Open                           |         |
|               +--------+------------------------------------------+         |
|               | bit 7  |   SD Screen Down                         |         |
|               +--------+------------------------------------------+         |
|               | bit 8  |   NAD closed                             |         |
|               +--------+------------------------------------------+         |
|               | bit 9  |   SDSM On                                |         |
|               +--------+------------------------------------------+         |
|               | bit 10 |   Radcooler Heaters On                   |         |
|               +--------+------------------------------------------+         |
|               | bit 11 |   Day mode bands telemetered at night    |         |
|               +--------+------------------------------------------+         |
|               | bit 12 |   Linear Emissive Calibration            |         |
|               +--------+------------------------------------------+         |
|               | bit 13 |   DC Restore Change                      |         |
|               +--------+------------------------------------------+         |
|               | bit 14 |   (Unused)                               |         |
|               +--------+------------------------------------------+         |
|               | bit 15 |   BB Heater On                           |         |
|               +--------+------------------------------------------+         |
|               | bit 16 |   Missing Previous Granule               |         |
|               +--------+------------------------------------------+         |
|               | bit 17 |   Missing Subsequent Granule             |         |
|               +--------+------------------------------------------+         |
|               |        |   SRCA calibration mode                  |         |
|               | bits   |       +--------+--------+-------------+  |         |
|               |  18-19 |       | bit 18 | bit 19 |  Meaning    |  |         |
|               |        |       +--------+--------+-------------+  |         |
|               |        |       |    0   |    0   | Radiometric |  |         |
|               |        |       |    0   |    1   | Spatial     |  |         |
|               |        |       |    1   |    0   | Spectral    |  |         |
|               |        |       |    1   |    1   | undetermined|  |         |
|               |        |       +--------+--------+-------------+  |         |
|               +--------+------------------------------------------+         |
|               | bit 20 |   moon in keep out box, any RSB          |         |
|               +--------+------------------------------------------+         |
|               | bit 21 |   moon in keep out box, any TEB          |         |
|               +--------+------------------------------------------+         |
|               | bit 22 |   All SV data are bad for any RSB        |         |
|               +--------+------------------------------------------+         |
|               | bit 23 |   All BB data are bad for any RSB        |         |
|               +--------+------------------------------------------+         |
|               | bit 24 |   Dropped scan(s) between leading and    |         |
|               |        |      middle granules                     |         |
|               +--------+------------------------------------------+         |
|               | bit 25 |   Dropped scan(s) between middle and     |         |
|               |        |      trailing granules                   |         |
|               +--------+------------------------------------------+         |
|               | bit 26 |   Sci Abnormal                           |         |
|               +--------+------------------------------------------+         |
|               | bit 27 |   (Remaining bits reserved for           |         |
|               | ... 31 |    future use)                           |         |
|               +========+==========================================+         |
+=============================================================================+
                  Dimension names for the next few SDSs
      name                         value             meaning
 ----------------------------    ---------   --------------------------
 "number of scans"               (variable)  same as "nscans"
 "number of 250m bands"              2       same as "Band_250m"
 "detectors per 250m band"          40       number detectors/250m band
 "number of 500m bands"              5       same as "Band_500m"
 "detectors per 500m band"          20       number detectors/500m band
 "number of 1km reflective bands"   15       number of 1km RefSB bands
 "number of emissive bands"         16       number of 1km Emiss bands
 "detectors per 1km band"           10       number detectors/1km band
+=============================================================================+
| SDS Name:    "Noise in Thermal Detectors"                                   |
| Data Type:   UINT8                                                          |
| Rank:        2                                                              |
| Dimensions:  ["number of emissive bands","detectors per 1km band"]          |
| Description: A measure of noise in each thermal emissive band detector.     |
|              Each value is computed as the average NEdL of a detector,      |
|              divided by the pre-launch NEdL value (from a LUT).  The        |
|              computed ratio in the range of [0,9] is scaled to a UINT8      |
|              (range of [0,255]).                                            |
+-----------------------------------------------------------------------------+
| SDS Name:    "Change in relative responses of thermal detectors"            |
| Data Type:   UINT8                                                          |
| Rank:        2                                                              |
| Dimensions:  ["number of emissive bands","detectors per 1km band"]          |
| Description: A measure of the change in each thermal emissive band          |
|              detector.  Each value is computed as the average linear        |
|              calibration coefficient b1 for a detector, divided by the      |
|              pre-launch value (from a LUT).  The computed ratio in the      |
|              range of [0,9] is scaled to a UINT8 (range of [0,255])         |
+-----------------------------------------------------------------------------+
| SDS Name:    "DC Restore Change for Thermal Bands"                          |
| Data Type:   INT8                                                           |
| Rank:        3                                                              |
| Dimensions:  ["number of scans","number of emissive bands",                 |
|               "detectors per 1km band"]                                     |
| Description: Identifies if a change occurred from the previous scan         |
|              (1 = change, 0 = no change) for thermal emissive band          |
|              detectors.                                                     |
+-----------------------------------------------------------------------------+
| SDS Name:    "DC Restore Change for Reflective 250m Bands"                  |
| Data Type:   INT8                                                           |
| Rank:        3                                                              |
| Dimensions:  ["number of scans","number of 250m bands",                     |
|               "detectors per 250m band"]                                    |
| Description: Identifies if a change occurred from the previous scan         |
|              (1 = change, 0 = no change) for 250m reflective solar band     |
|              detectors.                                                     |
+-----------------------------------------------------------------------------+
| SDS Name:    "DC Restore Change for Reflective 500m Bands"                  |
| Data Type:   INT8                                                           |
| Rank:        3                                                              |
| Dimensions:  ["number of scans","number of 500m bands",                     |
|               "detectors per 500m band"]                                    |
| Description: Identifies if a change occurred from the previous scan         |
|              (1 = change, 0 = no change) for 500m reflective solar band     |
|              detectors.                                                     |
+-----------------------------------------------------------------------------+
| SDS Name:    "DC Restore Change for Reflective 1km Bands"                   |
| Data Type:   INT8                                                           |
| Rank:        3                                                              |
| Dimensions:  ["number of scans","number of 1km reflective bands",           |
|               "detectors per 1km band"]                                     |
| Description: Identifies if a change occurred from the previous scan         |
|              (1 = change, 0 = no change) for 1km reflective solar band      |
|              detectors.                                                     |
+=============================================================================+

III)  HDF Vdatas

This section contains HDF Vdatas.  All are copied directly from the middle
Level 1A input granule (none result from any processing within Level 1B).
The field names can be deciphered by using the following SBRS document:

    "MODIS command, telemetry, science and engineering description",
    Hughes, Santa Barbara Research Center, May 14, 1997 (CDRL 305).

+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle All Part 1"                      |
| Number of fields:   19                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "SS_CP1553_MAJCYC"               UINT16         1                      |
|      "SS_CP_CMD_ECHO"                 UINT16         1                      |
|      "SS_CP_LOG_EVENT"                UINT16         1                      |
|      "SS_CP_MODE"                     UINT16         1                      |
|      "SS_CP_MODECHG_GO"               UINT16         1                      |
|      "SS_CP_SPARE"                    UINT16         1                      |
|      "SS_CP_VALTLM_COL"               UINT16         1                      |
|      "SS_CP_VALTLM_FMT"               UINT16         1                      |
|      "SS_DR_NAD_STEP"                 UINT16         1                      |
|      "SS_DR_SPARE"                    UINT16         1                      |
|      "SS_DR_SVD_STEP"                 UINT16         1                      |
|      "CS_FR_BBRADTAB"                 UINT16         1                      |
|      "SS_FR_LOG_EVENT"                UINT16         1                      |
|      "IR_PS1_INPUT_CUR"               UINT16         1                      |
|      "IR_PS2_INPUT_CUR"               UINT16         1                      |
|      "TA_RC_LWIR_CFPA"                UINT16         1                      |
|      "TA_RC_SMIR_CFPA"                UINT16         1                      |
|      "SR_SA_APX_PERIOD"               UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle All Part 2"                      |
| Number of fields:   19                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "CR_SM_MIR_HOME_A"               UINT16         1                      |
|      "CR_SM_MIR_HOME_B"               UINT16         1                      |
|      "CS_SM_MIR_STEP"                 UINT16         1                      |
|      "CR_SR_GRAT_CH_A"                UINT16         1                      |
|      "CR_SR_GRAT_CH_B"                UINT16         1                      |
|      "CR_SR_GRAT_FH_A"                UINT16         1                      |
|      "CR_SR_GRAT_FH_B"                UINT16         1                      |
|      "CR_SR_IR_SRC_OFF"               UINT16         1                      |
|      "CR_SR_SISHTR_OFF"               UINT16         1                      |
|      "CS_FR_GAINTAB"                  UINT16         1                      |
|      "CR_SR_SLIT_HOMEA"               UINT16         1                      |
|      "CR_SR_SLIT_HOMEB"               UINT16         1                      |
|      "CR_SR_WHL_HOMEA"                UINT16         1                      |
|      "CR_SR_WHL_HOMEB"                UINT16         1                      |
|      "CS_SR_GRAT_STEP"                UINT16         1                      |
|      "CS_SR_SLIT_STEP"                UINT16         1                      |
|      "CS_SR_SRCWH_STEP"               UINT16         1                      |
|      "CS_SR_LAMPS"                    UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle All Part 3"                      |
| Number of fields:   11                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "SS_CP_LAST_EVENT"               UINT16         1                      |
|      "SS_FR_LAST_EVENT"               UINT16         1                      |
|      "SS_CP_TC1_DAYS"                 UINT16         1                      |
|      "SS_CP_TC2_MILLIS"               UINT16         1                      |
|      "SS_CP_TC3_MILLIS"               UINT16         1                      |
|      "SS_CP_TC4_MICROS"               UINT16         1                      |
|      "CS_FR_OFFSETTAB"                UINT16         1                      |
|      "SS_CP_MACRO_ID"                 UINT16         1                      |
|      "SS_CP_MACRO_ON"                 UINT16         1                      |
|      "SS_DR_SDD_STEP"                 UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 0 of 7"                          |
| Number of fields:   15                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "CR_BB_A_PWR_ON"                 UINT16         1                      |
|      "CR_BB_B_PWR_ON"                 UINT16         1                      |
|      "SS_BB_DCYCLE"                   UINT16         1                      |
|      "CS_BB_TEMP_SET"                 UINT16         1                      |
|      "IR_BB_HTRA_CURH"                UINT16         1                      |
|      "IR_BB_HTRB_CURH"                UINT16         1                      |
|      "TP_BB_TEMP01H"                  UINT16         1                      |
|      "TP_BB_TEMP02H"                  UINT16         1                      |
|      "TP_BB_TEMP03H"                  UINT16         1                      |
|      "TP_BB_TEMP04H"                  UINT16         1                      |
|      "TP_BB_TEMP05H"                  UINT16         1                      |
|      "TP_BB_TEMP06H"                  UINT16         1                      |
|      "TP_BB_TEMP07H"                  UINT16         1                      |
|      "TP_BB_TEMP08H"                  UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 1 of 7"                          |
| Number of fields:   24                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "TP_BB_TEMP09H"                  UINT16         1                      |
|      "TP_BB_TEMP10H"                  UINT16         1                      |
|      "TP_BB_TEMP11H"                  UINT16         1                      |
|      "TP_BB_TEMP12H"                  UINT16         1                      |
|      "CR_CE_A_ON"                     UINT16         1                      |
|      "CR_CE_B_ON"                     UINT16         1                      |
|      "CR_CPA_EEP_WRE_M"               UINT16         1                      |
|      "CR_CPB_EEP_WRE_M"               UINT16         1                      |
|      "CR_CP_A_ON_M"                   UINT16         1                      |
|      "CR_CP_B_ON_M"                   UINT16         1                      |
|      "CR_CP_SET_TMF_A"                UINT16         1                      |
|      "CS_CP_SPARE"                    UINT16         1                      |
|      "SS_CP_IMOK_ON"                  UINT16         1                      |
|      "SS_CP_SPARE_1"                  UINT16         1                      |
|      "SS_CP_SPARE_2"                  UINT16         1                      |
|      "SS_CP_RESET_SRC"                UINT16         1                      |
|      "SS_CP_SCAN_EST"                 UINT16         1                      |
|      "SS_CP_LOG_STATE"                UINT16         1                      |
|      "SS_FR_LOG_STATE"                UINT16         1                      |
|      "SS_FR_SCIABNORM"                UINT16         1                      |
|      "CS_FR_SCI_QLKSET"               UINT16         1                      |
|      "CS_FR_ENG_QLKSET"               UINT16         1                      |
|      "SS_CP_STATUS_04"                UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 2 of 7"                          |
| Number of fields:   9                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "SS_CP_STATUS_05"                UINT16         1                      |
|      "SS_CP_STATUS_06"                UINT16         1                      |
|      "SS_CP_STATUS_07"                UINT16         1                      |
|      "SS_CP_STATUS_08"                UINT16         1                      |
|      "SS_CP_STATUS_09"                UINT16         1                      |
|      "SS_CP_UART_RESET"               UINT16         1                      |
|      "SS_CP_UART_HUNT"                UINT16         1                      |
|      "SS_CP_UART_SYNC"                UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 3A of 7"                         |
| Number of fields:   23                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "SS_CP_UART_NORM"                UINT16         1                      |
|      "SS_CP_SPARE"                    UINT16         1                      |
|      "SS_CP_TMF_GONE"                 UINT16         1                      |
|      "CR_DR_DRV_ON"                   UINT16         1                      |
|      "CR_DR_FS_ENABL_M"               UINT16         1                      |
|      "CR_DR_FS_SW_CLSD"               UINT16         1                      |
|      "CR_DR_NAD_CLSD"                 UINT16         1                      |
|      "CR_DR_NAD_FS_ON"                UINT16         1                      |
|      "CR_DR_NAD_OPEN"                 UINT16         1                      |
|      "CR_DR_PRI_FS_SEL"               UINT16         1                      |
|      "CR_DR_SDD_CLSD"                 UINT16         1                      |
|      "CR_DR_SDD_DRV_A"                UINT16         1                      |
|      "CR_DR_SDD_OPEN"                 UINT16         1                      |
|      "CR_DR_SDFS_DRVON"               UINT16         1                      |
|      "CR_DR_SPARE"                    UINT16         1                      |
|      "CR_DR_SDS_OPEN"                 UINT16         1                      |
|      "CR_DR_SVD_CLSD"                 UINT16         1                      |
|      "CR_DR_SVD_FS_ON"                UINT16         1                      |
|      "CR_DR_SVD_OPEN"                 UINT16         1                      |
|      "CR_DR_UNLACH_AON"               UINT16         1                      |
|      "CR_DR_UNLACH_BON"               UINT16         1                      |
|      "CS_DR_SVD_AT_OG"                UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 3B of 7"                         |
| Number of fields:   22                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "TP_DR_NAD_FS"                   UINT16         1                      |
|      "TP_DR_SPARE"                    UINT16         1                      |
|      "TP_DR_SVD_FS"                   UINT16         1                      |
|      "CR_FI_A_ON"                     UINT16         1                      |
|      "CR_FI_PORT_A_ON"                UINT16         1                      |
|      "CR_FI_A_RESET"                  UINT16         1                      |
|      "CR_FI_B_ON"                     UINT16         1                      |
|      "CR_FI_PORT_B_ON"                UINT16         1                      |
|      "CR_FI_B_RESET"                  UINT16         1                      |
|      "CR_FO_BLK1_ON"                  UINT16         1                      |
|      "CR_FO_SPARE_1"                  UINT16         1                      |
|      "CR_FO_BLK2_ON"                  UINT16         1                      |
|      "CR_FO_SPARE_2"                  UINT16         1                      |
|      "CR_FO_BLK3_ON"                  UINT16         1                      |
|      "CR_FO_SPARE_3"                  UINT16         1                      |
|      "CR_FO_BLK4_ON"                  UINT16         1                      |
|      "CS_FR_SCI_NORMAL"               UINT16         1                      |
|      "SR_FO_BLK1_MODE"                UINT16         1                      |
|      "SR_FO_BLK2_MODE"                UINT16         1                      |
|      "SR_FO_BLK3_MODE"                UINT16         1                      |
|      "SR_FO_BLK4_MODE"                UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 3C of 7"                         |
| Number of fields:   17                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "CR_FRA_EEP_WRE"                 UINT16         1                      |
|      "CR_FRB_EEP_WRE"                 UINT16         1                      |
|      "CR_FR_A_ON"                     UINT16         1                      |
|      "CR_FR_A_RESET"                  UINT16         1                      |
|      "CR_FR_B_ON"                     UINT16         1                      |
|      "CR_FR_B_RESET"                  UINT16         1                      |
|      "CS_FR_DAY_RATE"                 UINT16         1                      |
|      "CS_FR_DELAY_BB"                 UINT16         1                      |
|      "CS_FR_DELAY_EA"                 UINT16         1                      |
|      "CS_FR_DELAY_SD"                 UINT16         1                      |
|      "CS_FR_DELAY_SP"                 UINT16         1                      |
|      "CS_FR_DELAY_SR"                 UINT16         1                      |
|      "CR_SR_LAMPS_LOW"                UINT16         1                      |
|      "SR_FR_A_MODELONG"               UINT16         1                      |
|      "SR_FR_B_MODELONG"               UINT16         1                      |
|      "SS_FR_SPARE"                    UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 4A of 7"                         |
| Number of fields:   23                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "SS_FR_SPARE_1"                  UINT16         1                      |
|      "SS_FR_PKT_TYPE"                 UINT16         1                      |
|      "SS_FR_RESET_SRC"                UINT16         1                      |
|      "SS_FR_SPARE_2"                  UINT16         1                      |
|      "CS_PC_3136_SRISE"               UINT16         1                      |
|      "CS_PC_3235_SFALL"               UINT16         1                      |
|      "CS_PC_3334_SRISE"               UINT16         1                      |
|      "CS_PC_3334_SFALL"               UINT16         1                      |
|      "CS_PC_3235_SRISE"               UINT16         1                      |
|      "CS_PC_3136_SFALL"               UINT16         1                      |
|      "CR_PCLWA_ECAL_ON"               UINT16         1                      |
|      "CR_PCLWB_ECAL_ON"               UINT16         1                      |
|      "CR_PCLW_A_ON"                   UINT16         1                      |
|      "CR_PCLW_B_ON"                   UINT16         1                      |
|      "CR_PS1SHDN_ENA_M"               UINT16         1                      |
|      "CR_PS2SHDN_ENA_M"               UINT16         1                      |
|      "CR_PVLWA_CSUB_ON"               UINT16         1                      |
|      "CR_PVLWA_ECAL_ON"               UINT16         1                      |
|      "CR_PVLWB_CSUB_ON"               UINT16         1                      |
|      "CR_PVLWB_ECAL_ON"               UINT16         1                      |
|      "CR_PVLW_A_ON"                   UINT16         1                      |
|      "CR_PVLW_B_ON"                   UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 4B of 7"                         |
| Number of fields:   26                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "CR_PVLW_S_DELAYH"               UINT16         1                      |
|      "CR_PVNIRA_ECALON"               UINT16         1                      |
|      "CR_PVNIRB_ECALON"               UINT16         1                      |
|      "CR_PVNIR_A_ON"                  UINT16         1                      |
|      "CR_PVNIR_B_ON"                  UINT16         1                      |
|      "CR_PVNIR_S_DELYH"               UINT16         1                      |
|      "CR_PVSMA_CSUB_ON"               UINT16         1                      |
|      "CR_PVSMA_ECAL_ON"               UINT16         1                      |
|      "CR_PVSMB_CSUB_ON"               UINT16         1                      |
|      "CR_PVSMB_ECAL_ON"               UINT16         1                      |
|      "CR_PVSM_A_ON"                   UINT16         1                      |
|      "CR_PVSM_B_ON"                   UINT16         1                      |
|      "CR_PVSM_S_DELAYH"               UINT16         1                      |
|      "CR_PVVISA_ECALON"               UINT16         1                      |
|      "CR_PVVISB_ECALON"               UINT16         1                      |
|      "CR_PVVIS_A_ON"                  UINT16         1                      |
|      "CR_PVVIS_B_ON"                  UINT16         1                      |
|      "CR_PVVIS_S_DELYH"               UINT16         1                      |
|      "CR_PV_A_MEM_RAM"                UINT16         1                      |
|      "CR_PV_B_MEM_RAM"                UINT16         1                      |
|      "CR_PV_ECAL_ENA_A"               UINT16         1                      |
|      "CR_PV_ECAL_ENA_B"               UINT16         1                      |
|      "VR_PVLW_VCALH"                  UINT16         1                      |
|      "CS_FR_PC_DCR_ON"                UINT16         1                      |
|      "CS_FR_PV_DCR_ON"                UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 5A of 7"                         |
| Number of fields:   16                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_PVNIR_VCALH"                 UINT16         1                      |
|      "VR_PVSM_VCALH"                  UINT16         1                      |
|      "VR_PVVIS_VCALH"                 UINT16         1                      |
|      "CR_RC_CFPA_T1SET"               UINT16         1                      |
|      "CR_RC_CFPA_T3SET"               UINT16         1                      |
|      "CR_RC_CSHTR_ON"                 UINT16         1                      |
|      "CR_RC_CSTLM_ON"                 UINT16         1                      |
|      "CR_RC_ISHTR_ON"                 UINT16         1                      |
|      "CR_RC_ISTLM_ON"                 UINT16         1                      |
|      "CR_RC_LWHTR_ON"                 UINT16         1                      |
|      "CR_RC_LWTLM_ON"                 UINT16         1                      |
|      "CR_RC_OSHTR_ON"                 UINT16         1                      |
|      "CR_RC_OSTLM_ON"                 UINT16         1                      |
|      "CR_RC_SMHTR_ON"                 UINT16         1                      |
|      "CR_RC_SMTLM_ON"                 UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 5B of 7"                         |
| Number of fields:   19                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "CR_SA_A_HI_GAIN"                UINT16         1                      |
|      "CR_SA_A_SCAN_ON"                UINT16         1                      |
|      "CR_SA_B_HI_GAIN"                UINT16         1                      |
|      "CR_SA_B_SCAN_ON"                UINT16         1                      |
|      "SR_SA_A_PH_LOCK"                UINT16         1                      |
|      "SR_SA_B_PH_LOCK"                UINT16         1                      |
|      "CR_SM_SDSM_A_ON"                UINT16         1                      |
|      "CR_SM_SDSM_B_ON"                UINT16         1                      |
|      "CR_SR_A_ON"                     UINT16         1                      |
|      "CR_SR_B_ON"                     UINT16         1                      |
|      "CR_SR_SISFB_RAD"                UINT16         1                      |
|      "CR_SR_L_SHDN_ENA"               UINT16         1                      |
|      "IR_SR_10WLA_CURH"               UINT16         1                      |
|      "IR_SR_10WLB_CURH"               UINT16         1                      |
|      "IR_SR_1WLA_CURH"                UINT16         1                      |
|      "IR_SR_1WLB_CURH"                UINT16         1                      |
|      "TA_SR_IR_SRC_A"                 UINT16         1                      |
|      "TA_SR_IR_SRC_B"                 UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 6 of 7"                          |
| Number of fields:   17                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "TP_SR_GRAT_ELEX"                UINT16         1                      |
|      "TP_SR_GRAT_MOTOR"               UINT16         1                      |
|      "TP_SR_LAMP_RING"                UINT16         1                      |
|      "TP_SR_MIR2_DET"                 UINT16         1                      |
|      "VR_SR_LAMPS_H"                  UINT16         1                      |
|      "VR_SR_SRC_A_RADH"               UINT16         1                      |
|      "VR_SR_SRC_B_RADH"               UINT16         1                      |
|      "CR_TG_A_ON"                     UINT16         1                      |
|      "CR_TG_A_RESET"                  UINT16         1                      |
|      "CR_TG_B_ON"                     UINT16         1                      |
|      "CR_TG_B_RESET"                  UINT16         1                      |
|      "CS_FR_ENC_DELTA"                UINT16         1                      |
|      "CS_SR_USE_L10WX1"               UINT16         1                      |
|      "CS_SR_USE_L10WX2"               UINT16         1                      |
|      "CS_SR_USE_L10WX3"               UINT16         1                      |
|      "CS_SR_USE_L1WX1"                UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 7 of 7"                          |
| Number of fields:   3                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "TA_SR_SRC_A_SIPD"               UINT16         1                      |
|      "TA_SR_SRC_B_SIPD"               UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 0 of 63"                         |
| Number of fields:   6                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "TP_AO_LWIR_LENS"                UINT16         1                      |
|      "TP_AO_LWIR_OBJ"                 UINT16         1                      |
|      "TP_AO_PX_NZ_CORN"               UINT16         1                      |
|      "TP_AO_SMIR_LENS"                UINT16         1                      |
|      "TP_AO_SMIR_OBJ"                 UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 1 of 63"                         |
| Number of fields:   8                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "TP_AO_VNDICH_HSG"               UINT16         1                      |
|      "TP_CE_CAL2"                     UINT16         1                      |
|      "TP_CP_A_1553"                   UINT16         1                      |
|      "TP_CP_B_1553"                   UINT16         1                      |
|      "VR_CP_N11V"                     UINT16         1                      |
|      "VR_CP_N5V"                      UINT16         1                      |
|      "VR_CP_P11V"                     UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 2 of 63"                         |
| Number of fields:   8                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_CP_P5V"                      UINT16         1                      |
|      "TP_DR_NAD"                      UINT16         1                      |
|      "TP_DR_SDD"                      UINT16         1                      |
|      "TP_DR_SVD"                      UINT16         1                      |
|      "TP_FR_A_ENGINE"                 UINT16         1                      |
|      "TP_FR_B_ENGINE"                 UINT16         1                      |
|      "TP_ME_CHAS_TOP"                 UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 3 of 63"                         |
| Number of fields:   8                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "TP_MF_CALBKHD_SR"               UINT16         1                      |
|      "TP_MF_CVR_OP_SR"                UINT16         1                      |
|      "TP_MF_NAD_APT_NX"               UINT16         1                      |
|      "TP_MF_NAD_APT_NY"               UINT16         1                      |
|      "TP_MF_NX_AOBKHD"                UINT16         1                      |
|      "TP_MF_PX_AOBKHD"                UINT16         1                      |
|      "TP_MF_SV_PORT"                  UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 4 of 63"                         |
| Number of fields:   8                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "TP_MF_TOP_BY_KM2"               UINT16         1                      |
|      "TP_MF_YZ_CALBKHD"               UINT16         1                      |
|      "TP_MF_Z_BKHD_BB"                UINT16         1                      |
|      "TA_PC_B31_MUX"                  UINT16         1                      |
|      "TA_PC_B32_MUX"                  UINT16         1                      |
|      "TA_PC_B33_MUX"                  UINT16         1                      |
|      "TA_PC_B34_MUX"                  UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 5 of 63"                         |
| Number of fields:   8                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "TA_PC_B35_MUX"                  UINT16         1                      |
|      "TA_PC_B36_MUX"                  UINT16         1                      |
|      "TP_PC_CLAM_MNT"                 UINT16         1                      |
|      "VR_PC_B31_GND"                  UINT16         1                      |
|      "VR_PC_B31_RN12V"                UINT16         1                      |
|      "VR_PC_B31_RN5V"                 UINT16         1                      |
|      "VR_PC_B31_RP12V"                UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 6 of 63"                         |
| Number of fields:   9                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_PC_B31_RP5V"                 UINT16         1                      |
|      "VR_PC_B32_GND"                  UINT16         1                      |
|      "VR_PC_B32_RN12V"                UINT16         1                      |
|      "VR_PC_B32_RN5V"                 UINT16         1                      |
|      "VR_PC_B32_RP12V"                UINT16         1                      |
|      "VR_PC_B32_RP5V"                 UINT16         1                      |
|      "VR_PC_B33_GND"                  UINT16         1                      |
|      "VR_PC_B33_RN12V"                UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 7 of 63"                         |
| Number of fields:   8                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_PC_B33_RN5V"                 UINT16         1                      |
|      "VR_PC_B33_RP12V"                UINT16         1                      |
|      "VR_PC_B33_RP5V"                 UINT16         1                      |
|      "VR_PC_B34_GND"                  UINT16         1                      |
|      "VR_PC_B34_RN12V"                UINT16         1                      |
|      "VR_PC_B34_RN5V"                 UINT16         1                      |
|      "VR_PC_B34_RP12V"                UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 8 of 63"                         |
| Number of fields:   8                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_PC_B34_RP5V"                 UINT16         1                      |
|      "VR_PC_B35_GND"                  UINT16         1                      |
|      "VR_PC_B35_RN12V"                UINT16         1                      |
|      "VR_PC_B35_RN5V"                 UINT16         1                      |
|      "VR_PC_B35_RP12V"                UINT16         1                      |
|      "VR_PC_B35_RP5V"                 UINT16         1                      |
|      "VR_PC_B36_GND"                  UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 9 of 63"                         |
| Number of fields:   8                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_PC_B36_RN12V"                UINT16         1                      |
|      "VR_PC_B36_RN5V"                 UINT16         1                      |
|      "VR_PC_B36_RP12V"                UINT16         1                      |
|      "VR_PC_B36_RP5V"                 UINT16         1                      |
|      "TP_PS1_CVTR_SW"                 UINT16         1                      |
|      "TP_PS1_DIODE_OUT"               UINT16         1                      |
|      "TP_PS1_DWNREG_SW"               UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 10 of 63"                        |
| Number of fields:   9                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "TP_PS1_PRELOAD"                 UINT16         1                      |
|      "TP_PS2_CVTR_SW"                 UINT16         1                      |
|      "TP_PS2_DIODE_OUT"               UINT16         1                      |
|      "TP_PS2_DWNREG_SW"               UINT16         1                      |
|      "TP_PS2_PRELOAD"                 UINT16         1                      |
|      "VR_PS1_N15V_A1ME"               UINT16         1                      |
|      "VR_PS1_N15V_A2AF"               UINT16         1                      |
|      "CS_CP_MODIS_MOD"                UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 11 of 63"                        |
| Number of fields:   9                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_PS1_N15V_A3AS"               UINT16         1                      |
|      "VR_PS1_N30V_A1ME"               UINT16         1                      |
|      "VR_PS1_N8V_A2"                  UINT16         1                      |
|      "VR_PS1_P15V_A1ME"               UINT16         1                      |
|      "VR_PS1_P15V_A2AF"               UINT16         1                      |
|      "VR_PS1_P15V_A3AS"               UINT16         1                      |
|      "VR_PS1_P30V_A1"                 UINT16         1                      |
|      "VR_PS1_P5_6V_D1"                UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 12 of 63"                        |
| Number of fields:   9                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_PS1_P88V_A1ME"               UINT16         1                      |
|      "VR_PS1_P8V_A2"                  UINT16         1                      |
|      "VR_PS1_P8V_N1ME"                UINT16         1                      |
|      "VR_PS2_N15V_A1ME"               UINT16         1                      |
|      "VR_PS2_N15V_A2AF"               UINT16         1                      |
|      "VR_PS2_N15V_A3AS"               UINT16         1                      |
|      "VR_PS2_N30V_A1ME"               UINT16         1                      |
|      "VR_PS2_N8V_A2"                  UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 13 of 63"                        |
| Number of fields:   9                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_PS2_P15V_A1ME"               UINT16         1                      |
|      "VR_PS2_P15V_A2AF"               UINT16         1                      |
|      "VR_PS2_P15V_A3AS"               UINT16         1                      |
|      "VR_PS2_P30V_A1"                 UINT16         1                      |
|      "VR_PS2_P5_6V_D1"                UINT16         1                      |
|      "VR_PS2_P88V_A1ME"               UINT16         1                      |
|      "VR_PS2_P8V_A2"                  UINT16         1                      |
|      "VR_PS2_P8V_N1ME"                UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 14 of 63"                        |
| Number of fields:   8                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "TA_PVLW_PWB4_10"                UINT16         1                      |
|      "TA_PVNIR_PWB2_8"                UINT16         1                      |
|      "TA_PVNIR_PWB3_9"                UINT16         1                      |
|      "TA_PVSM_PWB5_11"                UINT16         1                      |
|      "TA_PVSM_PWB6_12"                UINT16         1                      |
|      "TA_PVVIS_PWB1_7"                UINT16         1                      |
|      "VR_PVLW_P30V"                   UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 15 of 63"                        |
| Number of fields:   9                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_PVLW_RN11V"                  UINT16         1                      |
|      "VR_PVLW_RN5V"                   UINT16         1                      |
|      "VR_PVLW_RP11V"                  UINT16         1                      |
|      "VR_PVLW_RP5V"                   UINT16         1                      |
|      "VR_PVNIR_P30V"                  UINT16         1                      |
|      "VR_PVNIR_P5VD3_9"               UINT16         1                      |
|      "VR_PVNIR_RN11V28"               UINT16         1                      |
|      "VR_PVNIR_RN11V39"               UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 16 of 63"                        |
| Number of fields:   9                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_PVNIR_RN5V2_8"               UINT16         1                      |
|      "VR_PVNIR_RN5V3_9"               UINT16         1                      |
|      "VR_PVNIR_RP11V28"               UINT16         1                      |
|      "VR_PVNIR_RP11V39"               UINT16         1                      |
|      "VR_PVNIR_RP5V2_8"               UINT16         1                      |
|      "VR_PVNIR_RP5V3_9"               UINT16         1                      |
|      "VR_PVSM_P30V"                   UINT16         1                      |
|      "VR_PVSM_P5VD6_12"               UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 17 of 63"                        |
| Number of fields:   9                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_PVSM_RN11V511"               UINT16         1                      |
|      "VR_PVSM_RN11V612"               UINT16         1                      |
|      "VR_PVSM_RN5V5_11"               UINT16         1                      |
|      "VR_PVSM_RN5V6_12"               UINT16         1                      |
|      "VR_PVSM_RP11V511"               UINT16         1                      |
|      "VR_PVSM_RP11V612"               UINT16         1                      |
|      "VR_PVSM_RP5V5_11"               UINT16         1                      |
|      "VR_PVSM_RP5V6_12"               UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 18 of 63"                        |
| Number of fields:   8                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_PVVIS_P30V"                  UINT16         1                      |
|      "VR_PVVIS_RN11V"                 UINT16         1                      |
|      "VR_PVVIS_RN5V"                  UINT16         1                      |
|      "VR_PVVIS_RP11V"                 UINT16         1                      |
|      "VR_PVVIS_RP5V"                  UINT16         1                      |
|      "TA_RC_CS"                       UINT16         1                      |
|      "TA_RC_CS_OG"                    UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 19 of 63"                        |
| Number of fields:   6                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "TA_RC_IS"                       UINT16         1                      |
|      "TA_RC_IS_OG"                    UINT16         1                      |
|      "TA_RC_OS_OG"                    UINT16         1                      |
|      "TP_RC_MNT_RING"                 UINT16         1                      |
|      "TP_RC_SPARE"                    UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 20 of 63"                        |
| Number of fields:   7                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "TP_RC_SPARE"                    UINT16         1                      |
|      "SA_SPARE_1"                     UINT16         1                      |
|      "SA_SPARE_2"                     UINT16         1                      |
|      "TP_SA_A_MTR"                    UINT16         1                      |
|      "SS_SA_SPARE"                    UINT16         1                      |
|      "TP_SA_SPARE"                    UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 21 of 63"                        |
| Number of fields:   7                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "TP_SA_RCT1_MIR"                 UINT16         1                      |
|      "TP_SA_SPARE"                    UINT16         1                      |
|      "TP_SA_RCT2_MIR"                 UINT16         1                      |
|      "SA_SPARE"                       UINT16         1                      |
|      "VR_SA_A_MTR_TORQ"               UINT16         1                      |
|      "VR_SA_A_RN11V"                  UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 22 of 63"                        |
| Number of fields:   8                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_SA_A_RP11V"                  UINT16         1                      |
|      "SA_SPARE"                       UINT16         1                      |
|      "VR_SA_B_MTR_TORQ"               UINT16         1                      |
|      "VR_SA_B_RN11V"                  UINT16         1                      |
|      "VR_SA_B_RP11V"                  UINT16         1                      |
|      "TP_SD_SPARE"                    UINT16         1                      |
|      "TP_SM_DET_AMP3"                 UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 23 of 63"                        |
| Number of fields:   8                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "TP_SR_MONO_CHAS1"               UINT16         1                      |
|      "TP_SR_MONO_CHAS2"               UINT16         1                      |
|      "TP_SR_SNOUT"                    UINT16         1                      |
|      "TP_SS_SPARE"                    UINT16         1                      |
|      "VR_TC_CSCKT_PV"                 UINT16         1                      |
|      "VR_TC_ISCKT_PV"                 UINT16         1                      |
|      "VR_TC_LWCKT_NV"                 UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 24 of 63"                        |
| Number of fields:   9                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_TC_LWCKT_PV"                 UINT16         1                      |
|      "SS_TC_SPARE_1"                  UINT16         1                      |
|      "VR_TC_OSCKT_PV"                 UINT16         1                      |
|      "VR_TC_SMCKT_NV"                 UINT16         1                      |
|      "VR_TC_SMCKT_PV"                 UINT16         1                      |
|      "SS_TC_SPARE_2"                  UINT16         1                      |
|      "VR_TC_VISCKT_NV"                UINT16         1                      |
|      "VR_TC_VISCKT_PV"                UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 25 of 63"                        |
| Number of fields:   7                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "TP_TE_FOLD_MIR"                 UINT16         1                      |
|      "TP_TE_PRI_MIR"                  UINT16         1                      |
|      "TP_TE_SEC_MIR"                  UINT16         1                      |
|      "TP_TM_ANLG_CKT"                 UINT16         1                      |
|      "VR_TM_REF_ACT1_1"               UINT16         1                      |
|      "VR_TM_REF_ACT1_2"               UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 26 of 63"                        |
| Number of fields:   6                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_TM_REF_ACT1_3"               UINT16         1                      |
|      "VR_TM_REF_ACT2_1"               UINT16         1                      |
|      "VR_TM_REF_ACT2_2"               UINT16         1                      |
|      "VR_TM_REF_ACT3_1"               UINT16         1                      |
|      "VR_TM_REF_ACT4_1"               UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 27 of 63"                        |
| Number of fields:   6                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_TM_REF_ACT5_1"               UINT16         1                      |
|      "VR_TM_REF_ACT5_2"               UINT16         1                      |
|      "VR_TM_REF_ACT5_3"               UINT16         1                      |
|      "VR_TM_REF_ACT6_1"               UINT16         1                      |
|      "VR_TM_REF_ACT6_3"               UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 28 of 63"                        |
| Number of fields:   6                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_TM_REF_ACT7_1"               UINT16         1                      |
|      "VR_TM_REF_BB_1"                 UINT16         1                      |
|      "VR_TM_REF_BB_2"                 UINT16         1                      |
|      "VR_TM_REF_BB_3"                 UINT16         1                      |
|      "VR_TM_REF_PRT1"                 UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 29 of 63"                        |
| Number of fields:   6                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_TM_REF_PSV1"                 UINT16         1                      |
|      "VR_TM_REF_PSV2"                 UINT16         1                      |
|      "VR_TM_REF_PSV3"                 UINT16         1                      |
|      "VR_TM_REF_PSV4"                 UINT16         1                      |
|      "VR_TM_REF_PSV5"                 UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 30 of 63"                        |
| Number of fields:   6                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_TM_REF_ACTGND"               UINT16         1                      |
|      "TA_AO_VIS_FPA"                  UINT16         1                      |
|      "TA_AO_NIR_FPA"                  UINT16         1                      |
|      "VR_RC_LW_FPA_HTR"               UINT16         1                      |
|      "VR_RC_SM_FPA_HTR"               UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 31 of 63"                        |
| Number of fields:   5                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_TM_REF_PRT2"                 UINT16         1                      |
|      "VR_TM_REF_PSV6"                 UINT16         1                      |
|      "VR_TM_REF_PSV7"                 UINT16         1                      |
|      "VR_TM_REF_PSV8"                 UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Telemetry Major Cycle 32 of 63"                        |
| Number of fields:   5                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "IR_SA_A_ECDR_LED"               UINT16         1                      |
|      "IR_SA_B_ECDR_LED"               UINT16         1                      |
|      "VR_SA_A_ECDR_MON"               UINT16         1                      |
|      "VR_SA_B_ECDR_MON"               UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Current S/C Ancillary Data"                            |
| Number of fields:   32                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "PACKET_HEADER"                  UINT8          6                      |
|      "TIME_STAMP"                     UINT8          8                      |
|      "FLAG_BYTE"                      UINT8          1                      |
|      "TIME_CONVERSION"                INT32          1                      |
|      "S/C_POSITION_X"                 INT32          1                      |
|      "S/C_POSITION_Y"                 INT32          1                      |
|      "S/C_POSITION_Z"                 INT32          1                      |
|      "S/C_VELOCITY_X"                 INT32          1                      |
|      "S/C_VELOCITY_Y"                 INT32          1                      |
|      "S/C_VELOCITY_Z"                 INT32          1                      |
|      "RESERVED_ANGLE_ROLL"            INT8           1                      |
|      "ATTITUDE_ANGLE_ROLL"            INT16          1                      |
|      "RESERVED_ANGLE_PITCH"           INT8           1                      |
|      "ATTITUDE_ANGLE_PITCH"           INT16          1                      |
|      "RESERVED_ANGLE_YAW"             INT8           1                      |
|      "ATTITUDE_ANGLE_YAW"             INT16          1                      |
|      "RESERVED_RATE_ROLL"             INT8           1                      |
|      "ATTITUDE_RATE_ROLL"             INT16          1                      |
|      "RESERVED_RATE_PITCH"            INT8           1                      |
|      "ATTITUDE_RATE_PITCH"            INT16          1                      |
|      "RESERVED_RATE_YAW"              INT8           1                      |
|      "ATTITUDE_RATE_YAW"              INT16          1                      |
|      "MAGNETIC_COIL_CURRENT_X"        INT8           1                      |
|      "MAGNETIC_COIL_CURRENT_Y"        INT8           1                      |
|      "MAGNETIC_COIL_CURRENT_Z"        INT8           1                      |
|      "SOLAR_ARRAY_CURRENT"            UINT8          1                      |
|      "SOLAR_POSITION_X"               INT8           1                      |
|      "SOLAR_POSITION_Y"               INT8           1                      |
|      "SOLAR_POSITION_Z"               INT8           1                      |
|      "MOON_POSITION_X"                INT8           1                      |
|      "MOON_POSITION_Y"                INT8           1                      |
|      "MOON_POSITION_Z"                INT8           1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Prior S/C Ancillary Data"                              |
| Number of fields:   32                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "PACKET_HEADER"                  UINT8          6                      |
|      "TIME_STAMP"                     UINT8          8                      |
|      "FLAG_BYTE"                      UINT8          1                      |
|      "TIME_CONVERSION"                INT32          1                      |
|      "S/C_POSITION_X"                 INT32          1                      |
|      "S/C_POSITION_Y"                 INT32          1                      |
|      "S/C_POSITION_Z"                 INT32          1                      |
|      "S/C_VELOCITY_X"                 INT32          1                      |
|      "S/C_VELOCITY_Y"                 INT32          1                      |
|      "S/C_VELOCITY_Z"                 INT32          1                      |
|      "RESERVED_ANGLE_ROLL"            INT8           1                      |
|      "ATTITUDE_ANGLE_ROLL"            INT16          1                      |
|      "RESERVED_ANGLE_PITCH"           INT8           1                      |
|      "ATTITUDE_ANGLE_PITCH"           INT16          1                      |
|      "RESERVED_ANGLE_YAW"             INT8           1                      |
|      "ATTITUDE_ANGLE_YAW"             INT16          1                      |
|      "RESERVED_RATE_ROLL"             INT8           1                      |
|      "ATTITUDE_RATE_ROLL"             INT16          1                      |
|      "RESERVED_RATE_PITCH"            INT8           1                      |
|      "ATTITUDE_RATE_PITCH"            INT16          1                      |
|      "RESERVED_RATE_YAW"              INT8           1                      |
|      "ATTITUDE_RATE_YAW"              INT16          1                      |
|      "MAGNETIC_COIL_CURRENT_X"        INT8           1                      |
|      "MAGNETIC_COIL_CURRENT_Y"        INT8           1                      |
|      "MAGNETIC_COIL_CURRENT_Z"        INT8           1                      |
|      "SOLAR_ARRAY_CURRENT"            UINT8          1                      |
|      "SOLAR_POSITION_X"               INT8           1                      |
|      "SOLAR_POSITION_Y"               INT8           1                      |
|      "SOLAR_POSITION_Z"               INT8           1                      |
|      "MOON_POSITION_X"                INT8           1                      |
|      "MOON_POSITION_Y"                INT8           1                      |
|      "MOON_POSITION_Z"                INT8           1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Command Parameters"                                    |
| Number of fields:   70                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "SET_BB_HTR_TEMP"                UINT16         1                      |
|      "SET_CP_OPER_MODE"               UINT16         1                      |
|      "SET_FR_RATE"                    UINT16         1                      |
|      "SET_FR_SCI_APID"                UINT16         1                      |
|      "SET_SR_L10WX1"                  UINT16         1                      |
|      "SET_FR_PKT_TYPE"                UINT16         1                      |
|      "SET_PVVIS_VCAL"                 UINT16         1                      |
|      "SET_PVNIR_VCAL"                 UINT16         1                      |
|      "SET_PVSM_VCAL"                  UINT16         1                      |
|      "SET_PVLW_VCAL"                  UINT16         1                      |
|      "SET_PVVIS_ITWK_V"               UINT16         1                      |
|      "SET_PVNIR_ITWK_V"               UINT16         1                      |
|      "SET_PVSM_ITWK_V"                UINT16         1                      |
|      "SET_PVLW_ITWK_V"                UINT16         1                      |
|      "SET_PVSM_VDET_V"                UINT16         1                      |
|      "SET_PVLW_VDET_V"                UINT16         1                      |
|      "SET_FR_ENG_APID"                UINT16         1                      |
|      "SET_SR_L10WX2"                  UINT16         1                      |
|      "ENABLE_CP_IMOK"                 UINT16         1                      |
|      "SET_FR_BBRADTAB"                UINT16         1                      |
|      "SET_FR_OFFSETTAB"               UINT16         1                      |
|      "SET_FR_GAINTAB"                 UINT16         1                      |
|      "TEST_FR_BBRAD"                  UINT16         1                      |
|      "SET_FR_SCI_QLK"                 UINT16         1                      |
|      "SET_FR_SR_DELAY"                UINT16         1                      |
|      "SET_FR_BB_DELAY"                UINT16         1                      |
|      "SET_CP_TMF_BUS"                 UINT16         1                      |
|      "OPEN_DR_UL_LOCK"                UINT16         1                      |
|      "SET_FR_SD_DELAY"                UINT16         1                      |
|      "SET_FR_SP_DELAY"                UINT16         1                      |
|      "SET_PV_MEM"                     UINT16         1                      |
|      "SET_FR_ENG_QLK"                 UINT16         1                      |
|      "CP_SPARE"                       UINT16         1                      |
|      "SET_SR_3L10W"                   UINT16         1                      |
|      "SET_FR_ENC_DELTA"               UINT16         1                      |
|      "SET_PVSMIR_ECAL"                UINT16         1                      |
|      "SET_PVLW_ECAL"                  UINT16         1                      |
|      "TEST_FR_PCOFFSET"               UINT16         1                      |
|      "TEST_FR_PVOFFSET"               UINT16         1                      |
|      "TEST_FR_PVGAIN"                 UINT16         1                      |
|      "SET_PVVIS_NSTEP"                UINT16         1                      |
|      "SET_PVNIR_NSTEP"                UINT16         1                      |
|      "SET_PVSM_NSTEP"                 UINT16         1                      |
|      "SET_PVLW_NSTEP"                 UINT16         1                      |
|      "SET_FR_EA_DELAY"                UINT16         1                      |
|      "SET_PVSMIR_CSUB"                UINT16         1                      |
|      "SET_PVLW_CSUB"                  UINT16         1                      |
|      "SET_DR_SVD_UL"                  UINT16         1                      |
|      "SET_DR_NAD_UL"                  UINT16         1                      |
|      "SET_DR_SDD_UL"                  UINT16         1                      |
|      "SET_DR_SDD_FS"                  UINT16         1                      |
|      "SET_FR_PV_DCRCMP"               UINT16         1                      |
|      "SET_FR_PC_DCRCMP"               UINT16         1                      |
|      "FR_SPARE_1"                     UINT16         1                      |
|      "FR_SPARE_2"                     UINT16         1                      |
|      "SET_SR_SIPD_HTR"                UINT16         1                      |
|      "SET_SR_SIS_FB"                  UINT16         1                      |
|      "SET_SR_LOV_SHDN"                UINT16         1                      |
|      "SET_SR_LAMPLEVEL"               UINT16         1                      |
|      "SET_SR_LAMPS"                   UINT16         1                      |
|      "SET_CP_LOG_STATE"               UINT16         1                      |
|      "SET_FR_LOG_STATE"               UINT16         1                      |
|      "SET_SR_IR_SRC"                  UINT16         1                      |
|      "FR_SPARE_3"                     UINT16         1                      |
|      "SET_SR_L1WX1"                   UINT16         1                      |
|      "SET_PVVIS_ECAL"                 UINT16         1                      |
|      "SET_PVNIR_ECAL"                 UINT16         1                      |
|      "SET_FR_SCIABNORM"               UINT16         1                      |
|      "FILL_BITS"                      UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Engineering BB data"                                   |
| Number of fields:   16                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "IR_BB_HTRA_CURR"                UINT16         1                      |
|      "IR_BB_HTRB_CURR"                UINT16         1                      |
|      "TP_BB_TEMP01"                   UINT16         1                      |
|      "TP_BB_TEMP02"                   UINT16         1                      |
|      "TP_BB_TEMP03"                   UINT16         1                      |
|      "TP_BB_TEMP04"                   UINT16         1                      |
|      "TP_BB_TEMP05"                   UINT16         1                      |
|      "TP_BB_TEMP06"                   UINT16         1                      |
|      "TP_BB_TEMP07"                   UINT16         1                      |
|      "TP_BB_TEMP08"                   UINT16         1                      |
|      "TP_BB_TEMP09"                   UINT16         1                      |
|      "TP_BB_TEMP10"                   UINT16         1                      |
|      "TP_BB_TEMP11"                   UINT16         1                      |
|      "TP_BB_TEMP12"                   UINT16         1                      |
|      "TP_BB_TEMP_AVG"                 UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Engineering CP valid format"                           |
| Number of fields:   2                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "SS_CP_VALENG_FMT"               UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Engineering FAM AF01 mux"                              |
| Number of fields:   11                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_PC_B31C10_DCR"               UINT16         1                      |
|      "VR_PC_B31C09_DCR"               UINT16         1                      |
|      "VR_PC_B31C08_DCR"               UINT16         1                      |
|      "VR_PC_B31C07_DCR"               UINT16         1                      |
|      "VR_PC_B31C06_DCR"               UINT16         1                      |
|      "VR_PC_B31C05_DCR"               UINT16         1                      |
|      "VR_PC_B31C04_DCR"               UINT16         1                      |
|      "VR_PC_B31C03_DCR"               UINT16         1                      |
|      "VR_PC_B31C02_DCR"               UINT16         1                      |
|      "VR_PC_B31C01_DCR"               UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Engineering FAM AF02 mux"                              |
| Number of fields:   11                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_PC_B32C10_DCR"               UINT16         1                      |
|      "VR_PC_B32C09_DCR"               UINT16         1                      |
|      "VR_PC_B32C08_DCR"               UINT16         1                      |
|      "VR_PC_B32C07_DCR"               UINT16         1                      |
|      "VR_PC_B32C06_DCR"               UINT16         1                      |
|      "VR_PC_B32C05_DCR"               UINT16         1                      |
|      "VR_PC_B32C04_DCR"               UINT16         1                      |
|      "VR_PC_B32C03_DCR"               UINT16         1                      |
|      "VR_PC_B32C02_DCR"               UINT16         1                      |
|      "VR_PC_B32C01_DCR"               UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Engineering FAM AF03 mux"                              |
| Number of fields:   11                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_PC_B33C10_DCR"               UINT16         1                      |
|      "VR_PC_B33C09_DCR"               UINT16         1                      |
|      "VR_PC_B33C08_DCR"               UINT16         1                      |
|      "VR_PC_B33C07_DCR"               UINT16         1                      |
|      "VR_PC_B33C06_DCR"               UINT16         1                      |
|      "VR_PC_B33C05_DCR"               UINT16         1                      |
|      "VR_PC_B33C04_DCR"               UINT16         1                      |
|      "VR_PC_B33C03_DCR"               UINT16         1                      |
|      "VR_PC_B33C02_DCR"               UINT16         1                      |
|      "VR_PC_B33C01_DCR"               UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Engineering FAM AF04 mux"                              |
| Number of fields:   11                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_PC_B34C10_DCR"               UINT16         1                      |
|      "VR_PC_B34C09_DCR"               UINT16         1                      |
|      "VR_PC_B34C08_DCR"               UINT16         1                      |
|      "VR_PC_B34C07_DCR"               UINT16         1                      |
|      "VR_PC_B34C06_DCR"               UINT16         1                      |
|      "VR_PC_B34C05_DCR"               UINT16         1                      |
|      "VR_PC_B34C04_DCR"               UINT16         1                      |
|      "VR_PC_B34C03_DCR"               UINT16         1                      |
|      "VR_PC_B34C02_DCR"               UINT16         1                      |
|      "VR_PC_B34C01_DCR"               UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Engineering FAM AF05 mux"                              |
| Number of fields:   11                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_PC_B35C10_DCR"               UINT16         1                      |
|      "VR_PC_B35C09_DCR"               UINT16         1                      |
|      "VR_PC_B35C08_DCR"               UINT16         1                      |
|      "VR_PC_B35C07_DCR"               UINT16         1                      |
|      "VR_PC_B35C06_DCR"               UINT16         1                      |
|      "VR_PC_B35C05_DCR"               UINT16         1                      |
|      "VR_PC_B35C04_DCR"               UINT16         1                      |
|      "VR_PC_B35C03_DCR"               UINT16         1                      |
|      "VR_PC_B35C02_DCR"               UINT16         1                      |
|      "VR_PC_B35C01_DCR"               UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Engineering FAM AF06 mux"                              |
| Number of fields:   11                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_PC_B36C10_DCR"               UINT16         1                      |
|      "VR_PC_B36C09_DCR"               UINT16         1                      |
|      "VR_PC_B36C08_DCR"               UINT16         1                      |
|      "VR_PC_B36C07_DCR"               UINT16         1                      |
|      "VR_PC_B36C06_DCR"               UINT16         1                      |
|      "VR_PC_B36C05_DCR"               UINT16         1                      |
|      "VR_PC_B36C04_DCR"               UINT16         1                      |
|      "VR_PC_B36C03_DCR"               UINT16         1                      |
|      "VR_PC_B36C02_DCR"               UINT16         1                      |
|      "VR_PC_B36C01_DCR"               UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Engineering Reg Sample Delay"                          |
| Number of fields:   5                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "CR_PVLW_S_DELAY"                UINT16         1                      |
|      "CR_PVNIR_S_DELAY"               UINT16         1                      |
|      "CR_PVSM_S_DELAY"                UINT16         1                      |
|      "CR_PVVIS_S_DELAY"               UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Engineering LWIR data"                                 |
| Number of fields:   8                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_PVLW_ITWKA"                  UINT16         1                      |
|      "VR_PVLW_VCAL"                   UINT16         1                      |
|      "VR_PVLW_VDDA"                   UINT16         1                      |
|      "VR_PVLW_VDDD"                   UINT16         1                      |
|      "VR_PVLW_VDDOUT"                 UINT16         1                      |
|      "VR_PVLW_VDET"                   UINT16         1                      |
|      "VR_PVLW_VPWELL"                 UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Engineering NIR data"                                  |
| Number of fields:   10                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_PVNIR_ITWKA"                 UINT16         1                      |
|      "VR_PVNIR_VCAL"                  UINT16         1                      |
|      "VR_PVNIR_VD1"                   UINT16         1                      |
|      "VR_PVNIR_VDDA"                  UINT16         1                      |
|      "VR_PVNIR_VDDD"                  UINT16         1                      |
|      "VR_PVNIR_VDDOUT"                UINT16         1                      |
|      "VR_PVNIR_VDET"                  UINT16         1                      |
|      "VR_PVNIR_VGUARD"                UINT16         1                      |
|      "VR_PVNIR_VPWELL"                UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Engineering SMIR data"                                 |
| Number of fields:   8                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_PVSM_ITWKA"                  UINT16         1                      |
|      "VR_PVSM_VCAL"                   UINT16         1                      |
|      "VR_PVSM_VDDA"                   UINT16         1                      |
|      "VR_PVSM_VDDD"                   UINT16         1                      |
|      "VR_PVSM_VDDOUT"                 UINT16         1                      |
|      "VR_PVSM_VDET"                   UINT16         1                      |
|      "VR_PVSM_VPWELL"                 UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Engineering VIS data"                                  |
| Number of fields:   10                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_PVVIS_ITWKA"                 UINT16         1                      |
|      "VR_PVVIS_VCAL"                  UINT16         1                      |
|      "VR_PVVIS_VD1"                   UINT16         1                      |
|      "VR_PVVIS_VDDA"                  UINT16         1                      |
|      "VR_PVVIS_VDDD"                  UINT16         1                      |
|      "VR_PVVIS_VDDOUT"                UINT16         1                      |
|      "VR_PVVIS_VDET"                  UINT16         1                      |
|      "VR_PVVIS_VGUARD"                UINT16         1                      |
|      "VR_PVVIS_VPWELL"                UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Engineering View Sample"                               |
| Number of fields:   4                                                       |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "SS_SM_VIEW_SMPL1"               UINT16         1                      |
|      "SS_SM_VIEW_SMPL2"               UINT16         1                      |
|      "SS_SM_VIEW_SMPL3"               UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Engineering SDSM data"                                 |
| Number of fields:   28                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "VR_SM01_SMPL1"                  UINT16         1                      |
|      "VR_SM01_SMPL2"                  UINT16         1                      |
|      "VR_SM01_SMPL3"                  UINT16         1                      |
|      "VR_SM02_SMPL1"                  UINT16         1                      |
|      "VR_SM02_SMPL2"                  UINT16         1                      |
|      "VR_SM02_SMPL3"                  UINT16         1                      |
|      "VR_SM03_SMPL1"                  UINT16         1                      |
|      "VR_SM03_SMPL2"                  UINT16         1                      |
|      "VR_SM03_SMPL3"                  UINT16         1                      |
|      "VR_SM04_SMPL1"                  UINT16         1                      |
|      "VR_SM04_SMPL2"                  UINT16         1                      |
|      "VR_SM04_SMPL3"                  UINT16         1                      |
|      "VR_SM05_SMPL1"                  UINT16         1                      |
|      "VR_SM05_SMPL2"                  UINT16         1                      |
|      "VR_SM05_SMPL3"                  UINT16         1                      |
|      "VR_SM06_SMPL1"                  UINT16         1                      |
|      "VR_SM06_SMPL2"                  UINT16         1                      |
|      "VR_SM06_SMPL3"                  UINT16         1                      |
|      "VR_SM07_SMPL1"                  UINT16         1                      |
|      "VR_SM07_SMPL2"                  UINT16         1                      |
|      "VR_SM07_SMPL3"                  UINT16         1                      |
|      "VR_SM08_SMPL1"                  UINT16         1                      |
|      "VR_SM08_SMPL2"                  UINT16         1                      |
|      "VR_SM08_SMPL3"                  UINT16         1                      |
|      "VR_SM09_SMPL1"                  UINT16         1                      |
|      "VR_SM09_SMPL2"                  UINT16         1                      |
|      "VR_SM09_SMPL3"                  UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Engineering SRCA data"                                 |
| Number of fields:   38                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "SR_SPARE"                       UINT16         1                      |
|      "IR_SR_10WLA_CURR"               UINT16         1                      |
|      "IR_SR_10WLB_CURR"               UINT16         1                      |
|      "IR_SR_1WLA_CURR"                UINT16         1                      |
|      "IR_SR_1WLB_CURR"                UINT16         1                      |
|      "VR_SR_LAMPS"                    UINT16         1                      |
|      "VR_SR_SELF_CAL1"                UINT16         1                      |
|      "VR_SR_SELF_CAL2"                UINT16         1                      |
|      "VR_SR_SELF_CAL3"                UINT16         1                      |
|      "VR_SR_SPCT_NORM1"               UINT16         1                      |
|      "VR_SR_SPCT_NORM2"               UINT16         1                      |
|      "VR_SR_SPCT_NORM3"               UINT16         1                      |
|      "VR_SR_SRC_A_RAD"                UINT16         1                      |
|      "VR_SR_SRC_B_RAD"                UINT16         1                      |
|      "CR_SR_A_ONE"                    UINT16         1                      |
|      "CR_SR_B_ONE"                    UINT16         1                      |
|      "CR_SR_GRAT_CH_AE"               UINT16         1                      |
|      "CR_SR_GRAT_CH_BE"               UINT16         1                      |
|      "CR_SR_GRAT_FH_AE"               UINT16         1                      |
|      "CR_SR_GRAT_FH_BE"               UINT16         1                      |
|      "CR_SR_IR_SRCOFFE"               UINT16         1                      |
|      "CR_SR_LAMPS_LOWE"               UINT16         1                      |
|      "CR_SR_LSHDN_ENAE"               UINT16         1                      |
|      "CR_SR_SISFB_RADE"               UINT16         1                      |
|      "CR_SR_SISHTROFFE"               UINT16         1                      |
|      "CR_SR_SLIT_HOMAE"               UINT16         1                      |
|      "CR_SR_SLIT_HOMBE"               UINT16         1                      |
|      "CR_SR_WHL_HOMEAE"               UINT16         1                      |
|      "CR_SR_WHL_HOMEBE"               UINT16         1                      |
|      "CS_SR_GRAT_STEPE"               UINT16         1                      |
|      "CS_SR_LAMPSE"                   UINT16         1                      |
|      "CS_SR_SLIT_STEPE"               UINT16         1                      |
|      "CS_SR_SRCWH_STPE"               UINT16         1                      |
|      "CS_SR_USEL10WX1E"               UINT16         1                      |
|      "CS_SR_USEL10WX2E"               UINT16         1                      |
|      "CS_SR_USEL10WX3E"               UINT16         1                      |
|      "CS_SR_USEL1WX1E"                UINT16         1                      |
+-----------------------------------------------------------------------------+
| Vdata Name:         "Engineering Temperature data"                          |
| Number of fields:   20                                                      |
|                                                                             |
|           field name                   type        order                    |
|      --------------------             ------       -----                    |
|      "LAST_VALID_SCAN"                UINT16         1                      |
|      "TA_SR_IR_SRC_AE"                UINT16         1                      |
|      "TA_SR_IR_SRC_BE"                UINT16         1                      |
|      "TA_SR_SRC_A_SPDE"               UINT16         1                      |
|      "TA_SR_SRC_B_SPDE"               UINT16         1                      |
|      "TP_SR_GRAT_ELEXE"               UINT16         1                      |
|      "TP_SR_GRAT_MTRE"                UINT16         1                      |
|      "TP_SR_LAMP_RINGE"               UINT16         1                      |
|      "TP_SR_MIR2_DETE"                UINT16         1                      |
|      "TP_SR_MONO_CHS1E"               UINT16         1                      |
|      "TP_SR_MONO_CHS2E"               UINT16         1                      |
|      "TP_SR_SNOUTE"                   UINT16         1                      |
|      "TA_AO_NIR_FPAE"                 UINT16         1                      |
|      "TA_AO_VIS_FPAE"                 UINT16         1                      |
|      "TA_RC_LWIR_CFPAE"               UINT16         1                      |
|      "TA_RC_SMIR_CFPAE"               UINT16         1                      |
|      "TP_MF_CALBHD_SRE"               UINT16         1                      |
|      "TP_SA_A_MTRE"                   UINT16         1                      |
|      "TP_SA_RCT1_MIRE"                UINT16         1                      |
|      "TP_SA_RCT2_MIRE"                UINT16         1                      |
|      "VR_TM_REF_BB_1E"                UINT16         1                      |
|      "VR_TM_REF_BB_2E"                UINT16         1                      |
|      "VR_TM_REF_BB_3E"                UINT16         1                      |
+=============================================================================+


