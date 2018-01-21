+=============================================================================+
|                                                                             |
|              L1B EV 250m File Specification -- MODIS/AQUA                   |
|                           V6.1.9, 11/15/2010                                |
|       Effective for all higher numbered versions until superceded           |
|                                                                             |
+=============================================================================+

   This document specifies the format for the MODIS/AQUA Level 1B 250m
   earth-view (EV) product. The product is implemented in an HDF file.

   DETECTOR ORDER CONVENTION FOR DETECTOR DEPENDENT DATA:
   The "SBRS" detector order convention refers to the numbering of detectors
   on each focal plane assembly as defined by the manufacturer.  In this
   convention, the numbering of detectors is opposite in direction from the
   satellite track direction.  The "product" detector order convention is the
   reverse of the SBRS convention.  The product detector order convention was
   defined to allow consecutive scans to be laid down on an image with the
   last detector of one scan abutting the first detector of the next scan.
   ALL DATA ITEMS WITHIN THIS FILE FOLLOW THE "PRODUCT" DETECTOR ORDER
   CONVENTION.

  +=====================================================================+
  |                                                                     |
  |  Contents:                                                          |
  |                                                                     |
  |       I)  Global Metadata                                           |
  |                                                                     |
  |      II)  Instrument and Uncertainty SDSs                           |
  |                                                                     |
  |     III)  Band-Subsetting SDSs                                      |
  |                                                                     |
  |      IV)  Geolocation SDSs                                          |
  |                                                                     |
  +=====================================================================+
 
  ______
 /\     \ +---------------------------------------------------------------+
/  \     \| +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ |
\   \_____\ |  I)  Global Metadata                                      | |
 \  /     / +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ |
  \/_____/+---------------------------------------------------------------+

   Description:
     This section contains granule-level metadata (i.e. ECS metadata, product
     granule metadata, QA granule metadata, HDF-EOS swath metadata, and
     L1B swath metadata).

   Contents:
     1.1) ECS Standard Core Granule Metadata
     1.2) ECS Standard Archive Granule Metadata
     1.3) MODIS Level 1B Product Granule Metadata
     1.4) MODIS Level 1B QA Granule Metadata
          1.4.1) QA Granule Metadata stored as global attributes
          1.4.2) QA Granule Metadata stored as SDSs
     1.5) Level 1B HDF-EOS SWATH Metadata 
     1.6) Level 1B Swath Metadata

   For attributes dimensioned 38 (number of bands), the band order is:
       1, 2, 3, .. , 13lo, 13hi, 14lo, 14hi, 15, .. , 36.
   For the reflective bands, total of 22, the order is:
       1, 2, 3, .. , 13lo, 13hi, 14lo, 14hi, .. , 19, 26.
   The order of the 490 detectors follows the order of the 38 bands described
   earlier. Within each band, the order of detectors is product order.

+=============================================================================+
|                                                                             |
|  1.1) ECS Standard Core Granule Metadata                                    |
|       Stored (as One ECS PVL String) in global attribute: CoreMetadata.0    |
|                                                                             |
+=============================================================================+
|                                                                             |
| Name                                        Example                         |
+-----------------------------------------------------------------------------+
| LOCALGRANULEID                              "MYD02QKM.A2001103.1555.002.    |
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
| AUTOMATICQUALITYFLAG.1                      "Passed"                        |
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
| PARAMETERNAME.1                             "EV_250_RefSB"                  |
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
| SHORTNAME                                   "MYD02QKM"   (Aqua)             |
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
| PGEVERSION                                  "2.1.1"                         |
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
| ASSOCIATEDPLATFORMSHORTNAME.1               "Aqua"                          |
+-----------------------------------------------------------------------------+
| ASSOCIATEDINSTRUMENTSHORTNAME.1             "MODIS"                         |
+-----------------------------------------------------------------------------+
| ASSOCIATEDSENSORSHORTNAME.1                 "MODIS"                         |
+=============================================================================+


+=============================================================================+
|                                                                             |
|  1.2) ECS Standard Archive Granule Metadata                                 |
|       Stored (as HDF ECS PVL string) in global attribute: ArchiveMetadata.0 |
|                                                                             |
+=============================================================================+
|                                                                             |
| Name                                        Example                         |
+-----------------------------------------------------------------------------+
| ALGORITHMPACKAGEACCEPTANCEDATE              "1998-04-01"                    |
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
|  1.3) MODIS Level 1B Product Granule Metadata                               |
|       Stored as global attributes.                                          |
|                                                                             |
+=============================================================================+
|                                    |         |     |                        |
| Name                               | Type    |Count| Example                |
+------------------------------------+---------+-----+------------------------+
| Number of Scans                    | int32   |   1 |  203 (see note 1)      |
+------------------------------------+---------+-----+------------------------+
| Number of Day mode scans           | int32   |   1 |  203                   |
+------------------------------------+---------+-----+------------------------+
| Number of Night mode scans         | int32   |   1 |  0                     |
+------------------------------------+---------+-----+------------------------+
| Incomplete Scans                   | int32   |   1 |  14                    |
+------------------------------------+---------+-----+------------------------+
| Max Earth View Frames              | int32   |   1 |  1354                  |
+------------------------------------+---------+-----+------------------------+
| %Valid EV Observations             | float32 |  38 |  98.2, 87.1, ...       |
+------------------------------------+---------+-----+------------------------+
| %Saturated EV Observations         | float32 |  38 |  1.4, 0.2, ...         |
+------------------------------------+---------+-----+------------------------+
| Electronics Redundancy Vector      | uint32  |   2 |  (see table 1 below)   |
+------------------------------------+---------+-----+------------------------+
| Electronics Configuration Change   | uint32  |   2 |  (see table 1 below)   |
+------------------------------------+---------+-----+------------------------+
| Reflective LUT Serial Number and   | char8   |  21 | "R001 1998:01:28:12:00"|
| Date of Last Change                |         |     | (RNNN YYYY:MM:DD:HH:MM)|
+------------------------------------+---------+-----+------------------------+
| Emissive LUT Serial Number and     | char8   |  21 | "E001 1998:01:28:12:00"|
| Date of Last Change                |         |     | (ENNN YYYY:MM:DD:HH:MM)|
+------------------------------------+---------+-----+------------------------+
| QA LUT Serial Number and           | char8   |  21 | "Q001 1998:01:28:12:00"|
| Date of Last Change                |         |     | (QNNN YYYY:MM:DD:HH:MM)|
+------------------------------------+---------+-----+------------------------+
| Focal Plane Set Point State        | int8    |   1 | 0=Running uncontrolled |
|                                    |         |     | 1=Set Point is 83 deg  |
|                                    |         |     | 2=Set Point is 85 deg  |
|                                    |         |     | 3=Set Point is 88 deg  |
+=============================================================================+
Note 1.  Usage of "Number of Scans":  The value of the global attribute
         "Number of Scans" is used in the dimensioning of many of the data
         items within this product.  In later descriptions of the dimensions
         of data items, the shortened string "nscans" is used to represent
         the value of the attribute "Number of Scans".  This value will
         typically be 203 or 204 during operations.
         However, the file format supports any value greater than zero.

+============================================================================+
| Table 1.                                                                   |
| "Electronics Redundancy Vector" and "Electronics Configuration Change"     |
|    (bit 0 is the least significant bit in the uint32 word)                 |
| The Electronics Redundancy Vector gives the final valid value for          |
|    the corresponding telemetry point within the granule.                   |
|    Values are 0=OFF and 1=ON.                                              |
| The Electronics Configuration Change identifies if any change occurred     |
|    within the granule for the corresponding telemetry point.               |
|    Values are 0 = no change, 1 = change.                                   |
+----------------------------------------------------------------------------+
| Word | Bit |             System                        | Telemetry point   |
+======+=====+===========================================+===================+
|  0   |  0  |  BB - Blackbody A                         | CR_BB_A_PWR_ON    |
|  0   |  1  |  BB - Blackbody B                         | CR_BB_B_PWR_ON    |
|  0   |  2  |  CE - Calibrator Electronics A            | CR_CE_A_ON        |
|  0   |  3  |  CE - Calibrator Electronics B            | CR_CE_B_ON        |
|  0   |  4  |  CP - Command & Telemetry Processor A     | CR_CP_A_ON_M      |
|  0   |  5  |  CP - Command & Telemetry Processor B     | CR_CP_B_ON_M      |
|  0   |  6  |  FI - FDDI Formatter A                    | CR_FI_A_ON        |
|  0   |  7  |  FI - FDDI Formatter B                    | CR_FI_B_ON        |
|  0   |  8  |  FO - FIFO Memory, Block 1                | CR_FO_BLK1_ON     |
|  0   |  9  |  FO - FIFO Memory, Block 2                | CR_FO_BLK2_ON     |
|  0   | 10  |  FO - FIFO Memory, Block 3                | CR_FO_BLK3_ON     |
|  0   | 11  |  FO - FIFO Memory, Block 4                | CR_FO_BLK4_ON     |
|  0   | 12  |  FR - Formatter A                         | CR_FR_A_ON        |
|  0   | 13  |  FR - Formatter B                         | CR_FR_B_ON        |
|  0   | 14  |  FR - PC DCR                              | CS_FR_PC_DCR_ON   |
|  0   | 15  |  FR - PV DCR                              | CS_FR_PV_DCR_ON   |
|  0   | 16  |  PC - PC Bands ADCs A                     | CR_PCLW_A_ON      |
|  0   | 17  |  PC - PC Bands ADCs B                     | CR_PCLW_B_ON      |
|  0   | 18  |  PV - PV FPA, LWIR ADCs A                 | CR_PVLW_A_ON      |
|  0   | 19  |  PV - PV FPA, LWIR ADCs B                 | CR_PVLW_B_ON      |
|  0   | 20  |  PV - PV FPA, SMIR ADCs A                 | CR_PVSM_A_ON      |
|  0   | 21  |  PV - PV FPA, SMIR ADCs B                 | CR_PVSM_B_ON      |
|  0   | 22  |  PV - PV FPA, NIR ADCs A                  | CR_PVNIR_A_ON     |
|  0   | 23  |  PV - PV FPA, NIR ADCs B                  | CR_PVNIR_B_ON     |
|  0   | 24  |  PV - PV FPA, VIS ADC A                   | CR_PVVIS_A_ON     |
|  0   | 25  |  PV - PV FPA, VIS ADC B                   | CR_PVVIS_B_ON     |
|  0   | 26  |  PV - PV LWIR ECAL A                      | CR_PVLWA_ECAL_ON  |
|  0   | 27  |  PV - PV LWIR ECAL B                      | CR_PVLWB_ECAL_ON  |
|  0   | 28  |  PV - PV NIR ECAL A                       | CR_PVNIRA_ECALON  |
|  0   | 29  |  PV - PV NIR ECAL B                       | CR_PVNIRB_ECALON  |
|  0   | 30  |  PV - PV SMIR ECAL A                      | CR_PVSMA_ECAL_ON  |
|  0   | 31  |  PV - PV SMIR ECAL B                      | CR_PVSMB_ECAL_ON  |
+------+-----+-------------------------------------------+-------------------+
|  1   |  0  |  PV - PV VIS ECAL A                       | CR_PVVISA_ECALON  |
|  1   |  1  |  PV - PV VIS ECAL B                       | CR_PVVISB_ECALON  |
|  1   |  2  |  RC - SMIR Temp. Control Heater           | CR_RC_SMHTR_ON    |
|  1   |  3  |  RC - LWIR Temp. Control Heater           | CR_RC_LWHTR_ON    |
|  1   |  4  |  SA - Scan Mirror Assembly A              | CR_SA_A_SCAN_ON   |
|  1   |  5  |  SA - Scan Mirror Assembly B              | CR_SA_B_SCAN_ON   |
|  1   |  6  |  SM - SD Stability Monitor A              | CR_SM_SDSM_A_ON   |
|  1   |  7  |  SM - SD Stability Monitor B              | CR_SM_SDSM_B_ON   |
|  1   |  8  |  SR - Spectroradiometric Assy (SRCA) A    | CR_SR_A_ON        |
|  1   |  9  |  SR - Spectroradiometric Assy (SRCA) B    | CR_SR_B_ON        |
|  1   | 10  |  TG - Timing Generator A                  | CR_TG_A_ON        |
|  1   | 11  |  TG - Timing Generator B                  | CR_TG_B_ON        |
+------+-----+-------------------------------------------+-------------------+
|  1   |12-31|  (Unused)                                                     |
+============================================================================+
Note:  When the most significant bit of the second word is set to 1, this
denotes that at least one of the telemetry fields had no valid values found.

+=============================================================================+
|                                                                             |
|  1.4)  MODIS Level 1B QA Granule Metadata                                   |
|                                                                             |
+=============================================================================+

+=============================================================================+
|                                                                             |
|  1.4.1) MODIS Level 1B QA Granule Metadata                                  |
|         stored as global attributes.                                        |
|                                                                             |
+=============================================================================+
|                                         |       |     |                     |
| Name                                     | Type  |Count| Example            |
+------------------------------------------+-------+-----+--------------------+
| Doors and Screens Configuration          | int8  |   1 | 2                  |
+------------------------------------------+-------+-----+--------------------+
| Reflective Bands With Bad Data           | int8  |  22 | 1,0,0,1...         |
+------------------------------------------+-------+-----+--------------------+
| Emissive Bands With Bad Data             | int8  |  16 | 1,0,0,1...         |
+------------------------------------------+-------+-----+--------------------+
| Noise in Black Body Thermistors          | uint8 |  12 | 10                 |
+------------------------------------------+-------+-----+--------------------+
| Noise in Average BB Temperature          | uint8 |   1 | 10                 |
+------------------------------------------+-------+-----+--------------------+
| Noise in LWIR FPA Temperature            | uint8 |   1 | 10                 |
+------------------------------------------+-------+-----+--------------------+
| Noise in MWIR FPA Temperature            | uint8 |   1 | 10                 |
+------------------------------------------+-------+-----+--------------------+
| Noise in Scan Mirror Thermistor #1       | uint8 |   1 | 10                 |
+------------------------------------------+-------+-----+--------------------+
| Noise in Scan Mirror Thermistor #2       | uint8 |   1 | 10                 |
+------------------------------------------+-------+-----+--------------------+
| Noise in Scan Mirror Thermistor Average  | uint8 |   1 | 10                 |
+------------------------------------------+-------+-----+--------------------+
| Noise in Instrument Temperature          | uint8 |   1 | 10                 |
+------------------------------------------+-------+-----+--------------------+
| Noise in Cavity Temperature              | uint8 |   1 | 10                 |
+------------------------------------------+-------+-----+--------------------+
| Noise in Temperature of NIR FPA          | uint8 |   1 | 10                 |
+------------------------------------------+-------+-----+--------------------+
| Noise in Temperature of Vis FPA          | uint8 |   1 | 10                 |
+------------------------------------------+-------+-----+--------------------+
| Dead Detector List                       | int8  | 490 | 1=True or 0=False  |
+------------------------------------------+-------+-----+--------------------+
| Noisy Detector List                      | int8  | 490 | 1=True or 0=False  |
+------------------------------------------+-------+-----+--------------------+
| Dead Subframe List                       | int8  | 520 | 1=True or 0=False  |
+------------------------------------------+-------+-----+--------------------+
| Noisy Subframe List                      | int8  | 520 | 1=True or 0=False  |
+------------------------------------------+-------+-----+--------------------+
| Detector Quality Flag                    | uint8 | 490 | (See table below)  |
+------------------------------------------+-------+-----+--------------------+
| Detector Quality Flag2                   | uint8 | 180 | (See table below)  |
+------------------------------------------+-------+-----+--------------------+
| Earth-Sun Distance                       |float32|   1 | 1.0045             |
+------------------------------------------+-------+-----+--------------------+
| Solar Irradiance on RSB Detectors over pi|float32| 330 | 511.459,511.46,... |
+------------------------------------------+-------+-----+--------------------+
| % L1A EV All Scan Data are Missing       |float32|   1 | 1.301              |
+------------------------------------------+-------+-----+--------------------+
| % L1A EV RSB DN Not in Day Mode          |float32| 490 | 1.301, 0.012, ...  |
+------------------------------------------+-------+-----+--------------------+
| % L1A EV DN Missing Within Scan          |float32| 490 | 1.301, 0.012, ...  |
+------------------------------------------+-------+-----+--------------------+
| % Dead Detector EV Data                  |float32| 490 | 1.301, 0.012, ...  |
+------------------------------------------+-------+-----+--------------------+
| % Dead Subframe EV Data                  |float32| 180 | 1.301, 0.012, ...  |
+------------------------------------------+-------+-----+--------------------+
| % Sector Rotation EV Data                |float32| 490 | 1.301, 0.012, ...  |
+------------------------------------------+-------+-----+--------------------+
| % Saturated EV Data                      |float32| 490 | 1.301, 0.012, ...  |
+------------------------------------------+-------+-----+--------------------+
| % TEB EV Data With Moon in SVP           |float32| 490 | 1.301, 0.012, ...  |
+------------------------------------------+-------+-----+--------------------+
| % EV Data Where Cannot Compute BG DN     |float32| 490 | 1.301, 0.012, ...  |
+------------------------------------------+-------+-----+--------------------+
| % RSB EV Data With dn** Below Scale      |float32| 490 | 1.301, 0.012, ...  |
+------------------------------------------+-------+-----+--------------------+
| % EV Data Where Nadir Door Closed        |float32| 490 | 1.301, 0.012, ...  |
+------------------------------------------+-------+-----+--------------------+
| % EV Data Not Calibrated                 |float32| 490 | 1.301, 0.012, ...  |
+------------------------------------------+-------+-----+--------------------+
| Bit QA Flags Last Value                  | uint32|   1 | 0                  |
+------------------------------------------+-------+-----+--------------------+
| Bit QA Flags Change                      | uint32|   1 | 0                  |
+------------------------------------------+-------+-----+--------------------+
| Granule Average QA Values                |float32|  50 | 290.004,289.984,...|
+=============================================================================+

+================================================================+
|  Detector Quality Flag                                         |
|  bitmask, each bit has value: 1=True or 0=False                |
|                                                                |
|  Bit |                       Meaning                           |
|======+=========================================================+
|   0  |                    Noisy Detector                       |
| (LSB)|   (same value as in attribute "Noisy Detector List")    |
|------+---------------------------------------------------------+
|   1  |                     Dead Detector                       |
|      |   (same value as in attribute "Dead Detector List")     |
|------+---------------------------------------------------------+
|   2  |                  Out-of-Family Gain                     |
|------+---------------------------------------------------------+
|   3  |                     Dynamic Range                       |
|------+---------------------------------------------------------+
|   4  |        Detector saturates on calibration source         |
|------+---------------------------------------------------------+
|   5  |             High calibration fit residuals              |
|------+---------------------------------------------------------+
|   6  |                 Electrical Crosstalk                    |
|------+---------------------------------------------------------+
|   7  |                          TBD                            |
+======+=========================================================+

+================================================================+
|  Detector Quality Flag2                                        |
|  bitmask, each bit has value: 1=True or 0=False                |
|                                                                |
|  Bit |                       Meaning                           |
|======+=========================================================+
|   0  |                    Noisy Subframe 1                     |
| (LSB)|                                                         |
|------+---------------------------------------------------------+
|   1  |                    Noisy Subframe 2                     |
|------+---------------------------------------------------------+
|   2  |                    Noisy Subframe 3                     |
|------+---------------------------------------------------------+
|   3  |                    Noisy Subframe 4                     |
|------+---------------------------------------------------------+
|   4  |                    Dead Subframe 1                      |
|------+---------------------------------------------------------+
|   5  |                    Dead Subframe 2                      |
|------+---------------------------------------------------------+
|   6  |                    Dead Subframe 3                      |
|------+---------------------------------------------------------+
|   7  |                    Dead Subframe 4                      |
+======+=========================================================+

+=============================================================================+
|                                                                             |
|  1.4.2) MODIS Level 1B QA Granule Metadata                                  |
|         stored as SDSs.                                                     |
|                                                                             |
+=============================================================================+
|                                    |        |                  |            |
| Name                               |  Type  | Dimensions       | Example    |
+------------------------------------+--------+------------------+------------+
| Noise in Thermal Detectors         |  uint8 | [16][10]         |  8,6,4...  |
+------------------------------------+--------+------------------+------------+
| Change in relative responses of    |        |                  |            |
| thermal detectors                  |  uint8 | [16][10]         |  10,6,8... |
+------------------------------------+--------+------------------+------------+
| DC Restore Change for Thermal Bands|  int8  | [nscans][16][10] |  0,0,1 ... |
+------------------------------------+--------+------------------+------------+
| DC Restore Change for Reflective   |        |                  |            |
| 250m Bands                         |  int8  | [nscans][2][40]  |  0,0,1...  |
+------------------------------------+--------+------------------+------------+
| DC Restore Change for Reflective   |        |                  |            |
| 500m Bands                         |  int8  | [nscans][5][20]  |  0,0,1...  |
+------------------------------------+--------+------------------+------------+
| DC Restore Change for Reflective   |        |                  |            |
| 1km Bands                          |  int8  | [nscans][15][10] |  0,0,1...  |
+=============================================================================+


+=============================================================================+
|                                                                             |
|  1.5)  Level 1B HDF-EOS SWATH Metadata                                      |
|        Stored as global attribute: StructMetadata.0                         |
|                                                                             |
+=============================================================================+
|                                                                             |
|GROUP=SwathStructure                                                         |
| GROUP=SWATH_1                                                               |
|                                                                             |
|  SwathName="MODIS_SWATH_Type_L1B"                                           |
|                                                                             |
|   GROUP=Dimension                                                           |
|    +-------------+---------------------+-----------------+                  |
|    | Object      | DimensionName       | Size            |                  |
|    +-------------+---------------------+-----------------+                  |
|    | Dimension_1 | "Band_250M"         | 2               |                  |
|    | Dimension_2 | "40*nscans"         | 40*nscans       |                  |
|    | Dimension_3 | "4*Max_EV_frames"   | 5416            |                  |
|    | Dimension_4 | "10*nscans"         | 10*nscans       |                  |
|    | Dimension_5 | "Max_EV_frames"     | 1354            |                  |
|    +-------------+---------------------+-----------------+                  |
|                                                                             |
|   GROUP=DimensionMap                                                        |
|    +----------------+-----------------+-------------------+--------+------+ |
|    | Object         | GeoDimension    | DataDimension     | Offset | Inc. | |
|    +----------------+-----------------+-------------------+--------+------+ |
|    | DimensionMap_1 | "10*nscans"     | "40*nscans"       |   1    |  4   | |
|    | DimensionMap_2 | "Max_EV_frames" | "4*Max_EV_frames" |   0    |  4   | |
|    +----------------+-----------------+-------------------+--------+------+ |
|                                                                             |
|   GROUP=GeoField                                                            |
|    +------------+-------------+--------------+------------------------------+
|    | Object     | GeoFieldName| DataType     | DimList                      |
|    +------------+-------------+--------------+------------------------------+
|    | GeoField_1 | "Latitude"  | DFNT_FLOAT32 |("10*nscans","Max_EV_frames") |
|    | GeoField_2 | "Longitude" | DFNT_FLOAT32 |("10*nscans","Max_EV_frames") |
|    +------------+-------------+--------------+------------------------------+
|                                                                             |
|   GROUP=DataField                                                           |
|    +-------------+-----------------------+-------------+------------------+ |
|    | Object      | DataFieldName         | DataType    | DimList          | |
|    +-------------+-----------------------+-------------+------------------+ |
|    | DataField_1 | "EV_250_RefSB"        | DFNT_UINT16 |("Band_250M",     | |
|    |             |                       |             | "40*nscans",     | |
|    |             |                       |             |"4*Max_EV_frames")| |
|    +-------------+-----------------------+-------------+------------------+ |
|    | DataField_2 | "EV_250_RefSB_        | DFNT_UINT8  |("Band_250M",     | |
|    |             |  Uncert_Indexes"      |             | "40*nscans",     | |
|    |             |                       |             |"4*Max_EV_frames")| |
|    +-------------+-----------------------+-------------+------------------+ |
|    | DataField_3 | "Band_250M"           | DFNT_FLOAT32| ("Band_250M")    | |
|    +-------------+-----------------------+-------------+------------------+ |
|                                                                             |
+=============================================================================+


+=============================================================================+
|                                                                             |
|  1.6) Level 1B Swath Metadata                                               |
|                                                                             |
|  These data are stored as a single HDF Vdata table having the name          |
|  "Level 1B Swath Metadata".  The number of records is equal to nscans       |
|  (determined by the "Number of Scans" global attribute).  The table         |
|  below describes the fields of this Vdata object.                           |
|                                                                             |
+=============================================================================+
|                                |         |     |                            |
|  Field Name                    |  Type   |Order| Range                      |
+--------------------------------+---------+-----+----------------------------+
|  Scan Number                   |  int32  |  1  | 1 to 204                   |
+--------------------------------+---------+-----+----------------------------+
|  Complete Scan Flag            |  int32  |  1  | 1=Complete, 0=Incomplete   |
+--------------------------------+---------+-----+----------------------------+
|  Scan Type                     |  char8  |  4  | "D   "=day,                |
|                                |         |     | "N   "=night,              |
|                                |         |     | "M   "=mixed,              |
|                                |         |     | "O   "=other               |
+--------------------------------+---------+-----+----------------------------+
|  Mirror Side                   | int32   |  1  | 0  or  1                   |
+--------------------------------+---------+-----+----------------------------+
|  EV Sector Start Time          | float64 |  1  | TAI time (number of sec.   |
|                                |         |     |           since 1/1/93)    |
+--------------------------------+---------+-----+----------------------------+
|  EV_Frames                     | int32   |  1  | 1354    (constant)         |
+--------------------------------+---------+-----+----------------------------+
|  Nadir_Frame_Number            | int32   |  1  | 677     (constant)         |
+--------------------------------+---------+-----+----------------------------+
|  Latitude of Nadir Frame       | float32 |  1  | -90.0 to 90.0 degrees      |
+--------------------------------+---------+-----+----------------------------+
|  Longitude of Nadir Frame      | float32 |  1  | -180.0 to 180.0 degrees    |
+--------------------------------+---------+-----+----------------------------+
|  Solar Azimuth of Nadir Frame  | float32 |  1  | -180.0 to 180.0 degrees    |
+--------------------------------+---------+-----+----------------------------+
|  Solar Zenith of Nadir Frame   | float32 |  1  | 0.0 to 180.0 degrees       |
+--------------------------------+---------+-----+----------------------------+
|  No. OBC BB thermistor outliers| int32   |  1  | 0 to 12                    |
+--------------------------------+---------+-----+----------------------------+
|  Bit QA Flags                  | uint32  |  1  | (see table below)          |
+--------------------------------+---------+-----+----------------------------+
|  Sector Rotation Angle         | float32 |  1  | 0.0 to 360.0 degrees       |
+================================+=========+=====+============================+

+=============================================================================+
|                                                                             |
|  1.7) Level 1B Swath Metadata Stored as global attributes.                  |
|                                                                             |
+=============================================================================+
|                                                                             |
| Name                                                        |Type   |#|const|
+-------------------------------------------------------------+-------+-+-----+
| HDFEOS_FractionalOffset_40*nscans_MODIS_SWATH_Type_L1B      |float32|1| 0.5 |
+-------------------------------------------------------------+-------+-+-----+
| HDFEOS_FractionalOffset_4*Max_EV_frames_MODIS_SWATH_Type_L1B|float32|1| 0.0 |
+=============================================================+=======+=+=====+

  Notes
    1. Scan Type has the three "extra" characters to make it a 32 bit field
       (to make reading the Vdata structure easier).

    2. Bit QA Flags are as follows (bit 0 is the least significant
       bit in the word):
       +========+==========================================+
       | bit #  |   flag name/desciption (1=True, 0=False) |
       +========+==========================================+
       | bit 0  |   Moon within defined limits of SVP      |
       +--------+------------------------------------------+
       | bit 1  |   Spacecraft Maneuver                    |
       +--------+------------------------------------------+
       | bit 2  |   Sector Rotation                        |
       +--------+------------------------------------------+
       | bit 3  |   Negative Radiance Beyond Noise Level   |
       +--------+------------------------------------------+
       | bit 4  |   PC Ecal on                             |
       +--------+------------------------------------------+
       | bit 5  |   PV Ecal on                             |
       +--------+------------------------------------------+
       | bit 6  |   SD Door Open                           |
       +--------+------------------------------------------+
       | bit 7  |   SD Screen Down                         |
       +--------+------------------------------------------+
       | bit 8  |   NAD closed                             |
       +--------+------------------------------------------+
       | bit 9  |   SDSM On                                |
       +--------+------------------------------------------+
       | bit 10 |   Radcooler Heaters On                   |
       +--------+------------------------------------------+
       | bit 11 |   Day mode bands telemetered at night    |
       +--------+------------------------------------------+
       | bit 12 |   Linear Emissive Calibration            |
       +--------+------------------------------------------+
       | bit 13 |   DC Restore Change                      |
       +--------+------------------------------------------+
       | bit 14 |   (Unused)                               |
       +--------+------------------------------------------+
       | bit 15 |   BB Heater On                           |
       +--------+------------------------------------------+
       | bit 16 |   Missing Previous Granule               |
       +--------+------------------------------------------+
       | bit 17 |   Missing Subsequent Granule             |
       +--------+------------------------------------------+
       |        |   SRCA calibration mode                  |
       | bits   |       +--------+--------+-------------+  |
       |  18-19 |       | bit 18 | bit 19 |  Meaning    |  |
       |        |       +--------+--------+-------------+  |
       |        |       |    0   |    0   | Radiometric |  |
       |        |       |    0   |    1   | Spatial     |  |
       |        |       |    1   |    0   | Spectral    |  |
       |        |       |    1   |    1   | undetermined|  |
       |        |       +--------+--------+-------------+  |
       +--------+------------------------------------------+
       | bit 20 |   moon in keep out box, any RSB          |
       +--------+------------------------------------------+
       | bit 21 |   moon in keep out box, any TEB          |
       +--------+------------------------------------------+
       | bit 22 |   All SV data are bad for any RSB        |
       +--------+------------------------------------------+
       | bit 23 |   All BB data are bad for any RSB        |
       +--------+------------------------------------------+
       | bit 24 |   Dropped scan(s) between leading and    |
       |        |      middle granules                     |
       +--------+------------------------------------------+
       | bit 25 |   Dropped scan(s) between middle and     |
       |        |      trailing granules                   |
       +--------+------------------------------------------+
       | bit 26 |   Sci Abnormal                           |
       +--------+------------------------------------------+
       | bit 27 |   (Remaining bits reserved for           |
       | ... 31 |    future use)                           |
       +========+==========================================+

  ______
 /\     \ +---------------------------------------------------------------+
/  \     \| +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ |
\   \_____\ |   II)  Instrument and Uncertainty SDSs                    | |
 \  /     / +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ |
  \/_____/+---------------------------------------------------------------+

   Description:
     This section contains the scientific data sets (SDSs) in which the 
     "corrected count" data (scaled to 16 bit integers) is stored and the 
     corresponding uncertainty indices.

   SDS Indexing:
     Each SDS is an array through a space of five parameters: band, detector,
     scan, frame and sample.  However, the SDS is written as a 3-D array so
     that it can be added to the HDF-EOS "SWATH" data structure:

        SDS_name[band_index, scan_and_detect_index, frame_and_sample_index]

     Where "band_index" represents an index through the bands represented in
     the SDS, "scan_and_detect_index" represents a double index through all
     detectors and scans, and "frame_and_sample_index" represents a double
     index through all frames and samples in a frame.  In this discussion,
     an "index" assumes a convention of 0 to N-1.  The following formula can
     be used to obtain the scan_and_detect_index from the detector index and
     the scan index:

        scan_and_det_index  =  scan_index * num_detectors_at_resolution  
                               + detector_index

     where num_detectors_at_resolution = 40 for SDSs in the 250m product.
     Similarly, to obtain the frame_and_sample_index from the frame and
     sample indices:

        frame_and_sample_index = frame_index * num_samples_at_resolution
                                 + sample_index

     where num_samples_at_resolution = 4 for SDSs in the 250m product.

   Specific Data Values:
     Generally, any scaled integer (SI) value above 32767 represents unusable
     data.  The following values are reserved for the scaled integer and the
     uncertainty index (UI) to signal the reason why data are unusable:

              Reason for unusable data                      SI      UI
     ----------------------------------------------       -----    ----
     Fill Value (includes reflective band data
       at night and completely missing L1A scans)         65535    255
     L1A DN is missing within a scan                      65534     15
     Detector is saturated                                65533     15
     Cannot compute  zero point DN                        65532     15
     Detector is dead                                     65531*    15
     RSB dn** below the minimum of the scaling range      65530     15
     TEB radiance or RSB dn** exceeds the
       maximum of the scaling range                       65529     15
     Aggregation algorithm failure                        65528     15
     Rotation of Earth-View Sector from
       nominal science collection position                65527     15
     Calibration coefficient b1 could not be computed     65526     15
     Subframe is dead                                     65525     15
     Both sides of PCLW electronics on                    65524     15
     (reserved for future use)                         65501-65523  15
     NAD closed upper limit                               65500     15

     * The "Dead Detector List" attribute (Section 1.4.1) identifies those
     detectors that are declared non-functional (cannot be calibrated).
     In collection 5, a valid SI value will be interpolated from the nearest
     functional detectors within the same scan if it is possible, but
     the value will be set at 65531 in collection 6.

     For the case of nadir-aperture door (NAD) closed, the SI value will be
     set to a value > 32767 (denoting unusable data).  First, if appropriate,
     a specific unusable data value above will be assigned.  If one of the
     specific reasons does not apply, then the SI value will be calculated
     normally and then the most significant bit will be flipped, resulting
     in a value > 32767.  If the resulting scaled integer value exceeds the
     value of NAD closed upper limit, the value will be set to that value.

     Also, for any uncertainty value other than fill-value (255), the
     uncertainty is written into the 4 least significant bits of the 8-bit
     UI value (hence, 15 is the maximum 4-bit number).  

     To compute uncertainty in percent from the uncertainty index (UI), use
     the specified_uncertainty and scaling_factor for the given band in:

         % uncertainty = specified_uncertainty * exp(UI/scaling_factor)

     The band-dependent values of specified_uncertainty and scaling_factor
     are attached as attributes to each uncertainty index SDS.

   Contents:
     2.1)  EV_250_RefSB
     2.2)  EV_250_RefSB_Uncert_Indexes

+=============================================================================+
|  2.1)  EV_250_RefSB                                                         |
|                                                                             |
|  SDS Name                 Data Type   Dimensions                            |
+-----------------------------------------------------------------------------+
| "EV_250_RefSB"            uint16    (Band_250M, 40*nscans, 4*Max_EV_frames) |
|                                                                             |
|    Description:                                                             |
|      16 bit scaled integer array containing the 250m Earth View data.       |
|                                                                             |
|    SDS Attributes:                                                          |
|      +-----------------------------+---------+-----+------------------------+
|      | Attribute Name              | Type    |Count| Attribute Value        |
|      +-----------------------------+---------+-----+------------------------+
|      | units                       | string  |  1  | "none"                 |
|      +-----------------------------+---------+-----+------------------------+
|      | valid_range                 | uint16  |  2  | 0, 32767               |
|      +-----------------------------+---------+-----+------------------------+
|      | _FillValue                  | uint16  |  1  | 65535                  |
|      +-----------------------------+---------+-----+------------------------+
|      | long_name                   | string  |  1  |"Earth View 250M Reflect|
|      |                             |         |     |ive Solar Bands Scaled  |
|      |                             |         |     |Integers"               |
|      +-----------------------------+---------+-----+------------------------+
|      | band_names                  | string  |  1  | "1, 2"                 |
|      +-----------------------------+---------+-----+------------------------+
|      | radiance_scales             | float32 |  2  | x.f, x.f               |
|      +-----------------------------+---------+-----+------------------------+
|      | radiance_offsets            | float32 |  2  | x.f, x.f               |
|      +-----------------------------+---------+-----+------------------------+
|      | radiance_units              | string  |  1  |"Watts/m^2/micrometer/st|
|      |                             |         |     | eradian"               |
|      +-----------------------------+---------+-----+------------------------+
|      | reflectance_scales          | float32 |  2  | x.f, x.f               |
|      +-----------------------------+---------+-----+------------------------+
|      | reflectance_offsets         | float32 |  2  | x.f, x.f               |
|      +-----------------------------+---------+-----+------------------------+
|      | reflectance_units           | string  |  1  | "none"                 |
|      +-----------------------------+---------+-----+------------------------+
|      | corrected_counts_scales     | float32 |  2  | x.f, x.f               |
|      +-----------------------------+---------+-----+------------------------+
|      | corrected_counts_offsets    | float32 |  2  | x.f, x.f               |
|      +-----------------------------+---------+-----+------------------------+
|      | corrected_counts_units      | string  |  1  | "counts"               |
|      +-----------------------------+---------+-----+------------------------+
|                                                                             |
+=============================================================================+
  Note: "x.f" means the value of the attribute is variable (may change from  
        granule to granule). The values of all other attributes are constant.

+=============================================================================+
|  2.2)  EV_250_RefSB_Uncert_Indexes                                          |
|                                                                             |
|  SDS Name                                Data Type   Dimensions             |
+-----------------------------------------------------------------------------+
| "EV_250_RefSB_Uncert_Indexes"            uint8       (Band_250M, 40*nscans, |
|                                                       4*Max_EV_frames)      |
|                                                                             |
|    Description:                                                             |
|      8 bit unsigned integer array containing the uncertainty index values   |
|      corresponding to the 250m Earth View data stored in the "EV_250_RefSB" |
|      SDS.  The uncertainty index values are calculated from the percent     |
|      uncertainty in the reflectance product.                                |
|                                                                             |
|    SDS Attributes:                                                          |
|      +--------------------------+--------+-----+----------------------------+
|      | Attribute Name           | Type   |Count| Attribute Value            |
|      +--------------------------+--------+-----+----------------------------+
|      | long_name                | string |  1  |"Earth View 250M Reflective |
|      |                          |        |     | Solar Bands Uncertainty    |
|      |                          |        |     | Indexes"                   |
|      +--------------------------+--------+-----+----------------------------+
|      | units                    | string |  1  | "none"                     |
|      +--------------------------+--------+-----+----------------------------+
|      | valid_range              | uint8  |  2  | 0, 15                      |
|      +--------------------------+--------+-----+----------------------------+
|      | _FillValue               | uint8  |  1  | 255                        |
|      +--------------------------+--------+-----+----------------------------+
|      | specified_uncertainty    | float32|  2  | x.f, x.f                   |
|      +--------------------------+--------+-----+----------------------------+
|      | scaling_factor           | float32|  2  | x.f, x.f                   |
|      +--------------------------+--------+-----+----------------------------+
|      | uncertainty_units        | string |  1  | "percent"                  |
|      +--------------------------+--------+-----+----------------------------+
|                                                                             |
+=============================================================================+
  Note: "x.f" means the value of the attribute is variable (may change from
        granule to granule). The values of all other attributes are constant.

  ______
 /\     \ +---------------------------------------------------------------+
/  \     \| +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ |
\   \_____\ |   III)  Band-Subsetting SDSs                              | |
 \  /     / +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ |
  \/_____/+---------------------------------------------------------------+

   Description:
     This section contains the scientific data sets (SDSs) which hold
     the information used to accomplish band subsetting.  For the 250M
     product, there is 1 band-subsetting SDS.

+=============================================================================+
|  Band_250M                                                                  |
|                                                                             |
|  SDS Name                 Data Type        Dimensions                       |
+-----------------------------------------------------------------------------+
| "Band_250M"                float32         (Band_250M)                      |
|                                                                             |
|    Description:                                                             |
|      float32 array containing the band numbers used for band subsetting.    |
|      The values are:  1.0, 2.0                                              |
|                                                                             |
|    SDS Attributes:                                                          |
|      +-----------------+---------+-----+------------------------------------+
|      | Attribute Name  | Type    |Count| Attribute Value                    |
|      +-----------------+---------+-----+------------------------------------+
|      | long_name       | string  |  1  | "250M Band Numbers for Subsetting" |
|      +-----------------+---------+-----+------------------------------------+
|                                                                             |
+=============================================================================+

  ______
 /\     \ +---------------------------------------------------------------+
/  \     \| +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ |
\   \_____\ |  IV)  Geolocation SDSs                                    | |
 \  /     / +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ |
  \/_____/+---------------------------------------------------------------+

   Description:
     This section contains information about the geodetic position
     (latitude and longitude) used for visualization.

   Contents:
     4.1) Latitude
     4.2) Longitude

+=============================================================================+
|  4.1) Latitude                                                              |
|                                                                             |
|  SDS Name                    Data Type    Dimensions                        |
+-----------------------------------------------------------------------------+
| "Latitude"                   float32      (10*nscans, Max_EV_frames)        |
|                                                                             |
|    Description:                                                             |
|      32 bit signed floating point array containing the geodetic latitudes   |
|      for the center of the corresponding 1km Earth view frames. Given the   |
|      offset and increment defined in the DimensionMap of Section 1.5, the   |
|      latitude points map to Inc1*(s-1)+offset1+0.5 in track direction and   |
|      Inc2*(f-1)+offset2 in scan direction, where 0.5 is the fractional      |
|      offset defined in Section 1.7 and 's' and 'f' are scaneline and frame  |
|      number starting from 1. If the formula gives a value 5.5 in track      |
|      direction that means 1/2 between detector 5 and detector 6. (please    |
|      note that detector# starts from 0)                                     |
|                                                                             |
|    SDS Attributes:                                                          |
|      +--------------------------+--------+-----+----------------------------+
|      | Attribute Name           | Type   |Count| Attribute Value            |
|      +--------------------------+--------+-----+----------------------------+
|      | units                    | string |  1  | "degrees"                  |
|      +--------------------------+--------+-----+----------------------------+
|      | valid_range              | float32|  2  | -90.0, 90.0                |
|      +--------------------------+--------+-----+----------------------------+
|      | _FillValue               | float32|  1  | -999.9                     |
|      +--------------------------+--------+-----+----------------------------+
|                                                                             |
+=============================================================================+

+=============================================================================+
|  4.2) Longitude                                                             |
|                                                                             |
|  SDS Name                    Data Type    Dimensions                        |
+-----------------------------------------------------------------------------+
| "Longitude"                  float32      (10*nscans, Max_EV_frames)        |
|                                                                             |
|    Description:                                                             |
|      32 bit signed floating point array containing the geodetic longitudes  |
|      for the center of the corresponding 1km Earth view frames. Given the   |
|      offset and increment defined in the DimensionMap of Section 1.5, the   |
|      longitude points map to Inc1*(s-1)+offset1+0.5 in track direction and  |
|      Inc2*(f-1)+offset2 in scan direction, where 0.5 is the fractional      |
|      offset defined in Section 1.7 and 's' and 'f' are scaneline and frame  |
|      number starting from 1. If the formula gives a value 5.5 in track      |
|      direction that means 1/2 between detector 5 and detector 6. (please    |
|      note that detector# starts from 0)                                     |
|                                                                             |
|    SDS Attributes:                                                          |
|      +--------------------------+--------+-----+----------------------------+
|      | Attribute Name           | Type   |Count| Attribute Value            |
|      +--------------------------+--------+-----+----------------------------+
|      | units                    | string |  1  | "degrees"                  |
|      +--------------------------+--------+-----+----------------------------+
|      | valid_range              | float32|  2  | -180.0, 180.0              |
|      +--------------------------+--------+-----+----------------------------+
|      | _FillValue               | float32|  1  | -999.9                     |
|      +--------------------------+--------+-----+----------------------------+
|                                                                             |
+=============================================================================+


