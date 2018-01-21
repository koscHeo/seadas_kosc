#include <sensorInfo.h>

#include <stddef.h>
#include <strings.h>

#include <genutils.h>

/**
 * lookup the ID for a sensor
 *
 * @param sensorName name of the sensor to lookup
 * @return sensor ID for the given sensor, -1 if not found.
 */
int sensorName2SensorId(const char* name) {
    int i;

    // convert the number if an int was sent in
    if(isValidInt(name)) {
        i = atoi(name);
        if(i>=0 && i<SENSOR_NUM) {
            return i;
        }
    }

    if(strcasecmp(name, "modisa") == 0)
      return HMODISA;
    if(strcasecmp(name, "modist") == 0)
      return HMODIST;

    for(i=0; i<SENSOR_NUM; i++) {
        if(strcasecmp(sensorName[i], name) == 0)
            return i;
    }

    return -1;
}

/**
 * get the sensor name of the sensor ID
 *
 * @param sensorId ID to lookup
 * @return the name of the sensor, NULL if not found
 */
const char* sensorId2SensorName(int sensorId) {
    if(sensorId < 0)
        return NULL;
    if(sensorId >= SENSOR_NUM)
        return NULL;
    return sensorName[sensorId];
}

/**
 * get the platform name of the sensor ID
 *
 * @param sensorId ID to lookup
 * @return the name of the sensor, NULL if not found
 */
const char* sensorId2PlatformName(int sensorId) {
    if(sensorId < 0)
        return NULL;
    if(sensorId >= SENSOR_NUM)
        return NULL;
    return platformName[sensorId];
}


/**
 * lookup the sensorID for a sensor given instrument and platform
 *
 * @param instrument instrument to use for sensorID look up
 * @param platform platform to use for sensorID lookup
 * @return the matching sensorID
 */
int instrumentPlatform2SensorID(const char* instrument, const char* platform) {

    if(strcasecmp(instrument, "MODIS") == 0) {
        if(strcasecmp(platform, "Terra") == 0) {
            return HMODIST;
        } else if(strcasecmp(platform, "Aqua") == 0) {
            return HMODISA;
        } else {
            return -1;
        }
    }

    int i;
    for(i=0; i<SENSOR_NUM; i++) {
        if(strcasecmp(instrumentName[i], instrument) == 0)
            return i;
    }

    return -1;
}

/**
 * lookup the sensorID for a sensor given instrument and platform
 *
 * @param instrument instrument to use for sensorID look up
 * @param platform platform to use for sensorID lookup
 * @return the matching sensorID
 */
int instrumentPlatform2subSensorID(const char* instrument, const char* platform) {

    if(strcasecmp(instrument, "MODIS") == 0) {
        if(strcasecmp(platform, "Aqua") == 0) {
            return MODIS_AQUA;
        } else if(strcasecmp(platform, "Terra") == 0) {
            return MODIS_TERRA;
        } else {
            return -1;
        }
    }
    return -1;
}


/**
 * get the sensor name given instrument and platform
 *
 * @param instrument instrument to use for sensorID look up
 * @param platform platform to use for sensorID lookup
 * @return the name of the sensor (internal memory), NULL if not found
 */
const char* instrumentPlatform2SensorName(const char* instrument, const char* platform) {
    int sensorId = instrumentPlatform2SensorID(instrument, platform);
    if(sensorId < 0)
        return NULL;
    return sensorName[sensorId];
}

/**
 * lookup the subsensorID for a sensor given instrument and platform
 *
 * @param sensorID
 * @param otherstring is USUALLY "spatialResolution"
 * @return the matching subsensorID
 */
//int get_subSensorID(int sensorID, const char* otherstring) {
//    int subsensorID;
//
//    switch (sensorID) {
//    case SEAWIFS:
//        if (strcmp(otherstring, "4.5 km") == 0)
//            subsensorID = SEAWIFS_GAC;
//        else
//            subsensorID = SEAWIFS_LAC;
//        break;
//    case HMODISA:
//        subsensorID = MODISA_FULL;
//        break;
//    case HMODIST:
//        subsensorID = MODIST_FULL;
//        break;
//    default:
//        subsensorID = -1;
//    }
//
//    return subsensorID;
//}

