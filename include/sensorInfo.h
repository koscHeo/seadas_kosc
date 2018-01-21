/*
 * sensorInfo.h
 *
 *  Created on: Oct 29, 2013
 *      Author: dshea
 */

#ifndef SENSOR_INFO_H_
#define SENSOR_INFO_H_

#include <sensorDefs.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

static char sensorName[][20] =
    {"SeaWiFS",
    "MOS",
    "OCTS",
    "AVHRR",
    "OSMI",
    "MODIST",
    "MODISA",
    "CZCS",
    "HMODIST",
    "HMODISA",
    "OCM1",
    "OCM2",
    "MERIS",
    "MERIS",
    "VIIRSN",
    "OCRVC",
    "HICO",
    "GOCI",
    "OLI",
    "Aquarius",
    "ORCA",
    "AVIRIS",
    "PRISM",
    "OLCI"
    };

static char instrumentName[][20] =
    {"SeaWiFS",
     "MOS",
     "OCTS",
     "AVHRR",
     "OSMI",
     "MODIS",
     "MODIS",
     "CZCS",
     "MODIS",
     "MODIS",
     "OCM",
     "OCM-2",
     "MERIS",
     "MERIS",
     "VIIRS",
     "OCRVC",
     "HICO",
     "GOCI",
     "OLI",
     "Aquarius",
     "ORCA",
     "AVIRIS",
     "PRISM",
     "OLCI"
    };

static char platformName[][20] =
    {"Orbview-2",
     "IRS-P3",
     "ADEOS",
     "AVHRR",
     "KOMPSAT",
     "Terra",
     "Aqua",
     "Nimbus-7",
     "Terra",
     "Aqua",
     "IRS-P4",
     "Oceansat-2",
     "Envisat",
     "Envisat",
     "Suomi-NPP",
     "OCRVC",
     "ISS",
     "COMS",
     "Landsat-8",
     "SAC-D",
     "PACE",
     "AVIRIS",
     "PRISM",
     "OLCI"
    };

static char sensorDir[][20] =
    {"seawifs",
    "mos",
    "octs",
    "avhrr",
    "osmi",
    "modist",
    "modisa",
    "czcs",
    "hmodist",
    "hmodisa",
    "ocm1",
    "ocm2",
    "meris",
    "meris",
    "viirsn",
    "ocrvc",
    "hico",
    "goci",
    "oli",
    "aquarius",
    "orca",
    "aviris",
    "prism",
    "olci"
    };

static char sensorSub[][20] =
    {"seawifs_gac",
    "seawifs_lac",
    "modisa",
    "hmodisa",
    "modist",
    "hmodist"};


int32_t rdsensorinfo(int32_t sensorID, int32_t evalmask, const char *pname, void **pval);
int sensorName2SensorId(const char* sensorName);
int instrumentPlatform2SensorID(const char* instrument, const char* platform);
int instrumentPlatform2subSensorID(const char* instrument, const char* platform);
const char* sensorId2SensorName(int sensorId);
const char* sensorId2PlatformName(int sensorId);
const char* instrumentPlatform2SensorName(const char* instrument, const char* platform);

#ifdef __cplusplus
}
#endif


#endif /* SENSOR_INFO_H_ */
