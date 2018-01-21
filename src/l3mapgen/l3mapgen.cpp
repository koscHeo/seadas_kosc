/*
 * l3mapgen.cpp
 *
 *  Created on: Jan 15, 2014
 *      Author: dshea
 */
#include "l3mapgen.h"

#include "OutFile.h"
#include <L3FileSMI.h>

#include <iostream>

#include <genutils.h>
#include <string>

#include <proj_api.h>
#include <sensorInfo.h>
#include <productInfo.h>

using namespace std;
using namespace l3;

#define VERSION "1.0"

#define NUM_SEARCH_POINTS 51

enum MeasurementType { Avg, Stdev, Variance, Nobs, Nscenes, ObsTime, BinNum };

float percentDelta = 0.01;

// global params for program
int imageWidth = 0;
int imageHeight = 0;
bool doingQuality = false;
bool doingRGB = false;
string prodName;
vector<string> prodNameList;
vector<MeasurementType> prodMeasurementList;

clo_optionList_t* optionList = NULL;


void printStartInfo(OutFile* outFile) {
    if (want_verbose) {
        meta_l3bType* metaData = outFile->getMetadata();
        
        clo_printVersion();
        printf("ifile      : %s\n", clo_getString(optionList, "ifile"));
        printf("ofile      : %s\n", clo_getString(optionList, "ofile"));
        printf("oformat    : %s\n", clo_getString(optionList, "oformat"));
        if(clo_isSet(optionList, "ofile2")) {
            printf("ofile2     : %s\n", clo_getString(optionList, "ofile2"));
            printf("oformat2   : %s\n", clo_getString(optionList, "oformat2"));
        }
        printf("projection : %s\n", clo_getRawString(optionList, "projection"));
        printf("resolution : %.3fm\n", outFile->getResolution());
        if(doingRGB)
            printf("product_rgb: %s\n", prodName.c_str());
        else
            printf("product    : %s\n", prodName.c_str());
        printf("north      : %8.3f\n", metaData->north);
        printf("south      : %8.3f\n", metaData->south);
        printf("east       : %8.3f\n", metaData->east);
        printf("west       : %8.3f\n", metaData->west);
        float tmpf = clo_getFloat(optionList, "central_meridian");
        if(tmpf > -900.0) {
            printf("central_meridian : %8.3f\n", tmpf);
        }
        printf("image size : %d x %d\n", imageHeight, imageWidth);
        
        printf("\n");
    }
}

void printEndInfo(OutFile* outFile) {
    if (want_verbose) {
        meta_l3bType* metaData = outFile->getMetadata();
        
        printf("\n\n");
        printf("actual data min       : %f\n", outFile->getFileMinVal());
        printf("actual data max       : %f\n", outFile->getFileMaxVal());
        printf("num filled pixels     : %d\n", outFile->getNumFilledPixels());
        float tmp = outFile->getPercentFilledPixels();
        printf("percent filled pixels : %.2f%%\n", outFile->getPercentFilledPixels());
        
        printf("\n");
    }
}

InterpType interpStr2Type(const char* str) {
    string s = str;
    boost::trim(s);
    boost::to_lower(s);

    if(s.compare("bin") == 0)
       return Interp_Bin;
    if(s.compare("linear") == 0)
       return Interp_Linear;
    if(s.compare("area") == 0)
       return Interp_Area;
        
    return Interp_Nearest;
}

const char* interpType2Str(InterpType interp) {
    switch(interp) {
        case Interp_Nearest:
            return "nearest";
        case Interp_Bin:
            return "bin";
        case Interp_Linear:
            return "linear";
        case Interp_Area:
            return "area";
        default:
            return "unknown";
    }
}

/**
 * figure out if we want and can do quality processing
 * @param l3File input bin file
 * @param outFile output file
 */
bool setupQualityProcessing(L3File* l3File, OutFile* outFile, OutFile* outFile2) {
    doingQuality = true;
    clo_option_t* option = clo_findOption(optionList, "use_quality");
    if(clo_isOptionSet(option)) {
        if(clo_getOptionBool(option)) {
            if(l3File->hasQuality()) {
                doingQuality = true;
            } else {
                printf("-E- Quality processing was requested, but the input file does not have quality data.\n");
                exit(1);
            }
        } else {
            doingQuality = false;
        }
    } else {
        if(l3File->hasQuality())
            doingQuality = true;
        else
            doingQuality = false;
    }
    
    l3File->setQualityProcessing(doingQuality);
    outFile->setQualityProcessing(doingQuality);
    if(outFile2)
        outFile2->setQualityProcessing(doingQuality);
    return doingQuality;
}

void setupProduct(string &prodName, MeasurementType measure, OutFile* outFile, OutFile* outFile2) {
    // get the product info
    productInfo_t *p_info;
    p_info = allocateProductInfo();

    int sensorId = sensorName2SensorId(outFile->getMetadata()->sensor_name);
    if (sensorId == -1){
        printf("-E- Unknown sensor name %s\n", outFile->getMetadata()->sensor_name);
        exit(EXIT_FAILURE);
    }
    
    if (!findProductInfo(prodName.c_str(), sensorId, p_info)) {
      printf("-E- product %s not found in XML product table\n", prodName.c_str());
      exit(EXIT_FAILURE);
    }

    // now we have to fix the p_info structure
    // Avg, Stdev, Variance, Nobs, Nscenes, ObsTime, BinNum
    string tmpStr;
    switch(measure) {
        case Avg:
            break;
        case Stdev:
            tmpStr = p_info->description;
            free(p_info->description);
            tmpStr += " (Standard Deviation)";
            p_info->description = strdup(tmpStr.c_str());
            free(p_info->palette);
            p_info->palette = strdup(PRODUCT_DEFAULT_palette);
            free(p_info->dataType);
            p_info->dataType = strdup("float");
            tmpStr = p_info->suffix;
            free(p_info->suffix);
            tmpStr += "_stdev";
            p_info->suffix = strdup(tmpStr.c_str());
            free(p_info->displayScale);
            p_info->displayScale = strdup("linear");
            p_info->addOffset = 0.0;
            p_info->scaleFactor = 1.0;
            break;
        case Variance:
            tmpStr = p_info->description;
            free(p_info->description);
            tmpStr += " (Variance)";
            p_info->description = strdup(tmpStr.c_str());
            free(p_info->palette);
            p_info->palette = strdup(PRODUCT_DEFAULT_palette);
            free(p_info->dataType);
            p_info->dataType = strdup("float");
            tmpStr = p_info->suffix;
            free(p_info->suffix);
            tmpStr += "_var";
            p_info->suffix = strdup(tmpStr.c_str());
            free(p_info->displayScale);
            p_info->displayScale = strdup("linear");
            p_info->addOffset = 0.0;
            p_info->scaleFactor = 1.0;
            break;
        case Nobs:
            tmpStr = p_info->description;
            free(p_info->description);
            tmpStr += " (number of observations)";
            p_info->description = strdup(tmpStr.c_str());
            free(p_info->units);
            p_info->units = strdup("counts");
            free(p_info->palette);
            p_info->palette = strdup(PRODUCT_DEFAULT_palette);
            free(p_info->dataType);
            p_info->dataType = strdup("short");
            tmpStr = p_info->suffix;
            free(p_info->suffix);
            tmpStr += "_nobs";
            p_info->suffix = strdup(tmpStr.c_str());
            p_info->fillValue = PRODUCT_DEFAULT_fillValue;
            p_info->validMin = 0;
            p_info->validMax = 32767;
            free(p_info->displayScale);
            p_info->displayScale = strdup("linear");
            p_info->displayMin = PRODUCT_DEFAULT_displayMin;
            p_info->displayMax = PRODUCT_DEFAULT_displayMax;
            p_info->addOffset = 0.0;
            p_info->scaleFactor = 1.0;
            break;
        case Nscenes:
            tmpStr = p_info->description;
            free(p_info->description);
            tmpStr += " (number of scenes)";
            p_info->description = strdup(tmpStr.c_str());
            free(p_info->units);
            p_info->units = strdup("counts");
            free(p_info->palette);
            p_info->palette = strdup(PRODUCT_DEFAULT_palette);
            free(p_info->dataType);
            p_info->dataType = strdup("short");
            tmpStr = p_info->suffix;
            free(p_info->suffix);
            tmpStr += "_nscenes";
            p_info->suffix = strdup(tmpStr.c_str());
            p_info->fillValue = PRODUCT_DEFAULT_fillValue;
            p_info->validMin = 0;
            p_info->validMax = 32767;
            free(p_info->displayScale);
            p_info->displayScale = strdup("linear");
            p_info->displayMin = PRODUCT_DEFAULT_displayMin;
            p_info->displayMax = PRODUCT_DEFAULT_displayMax;
            p_info->addOffset = 0.0;
            p_info->scaleFactor = 1.0;
            break;
        case ObsTime:
            tmpStr = p_info->description;
            free(p_info->description);
            tmpStr += " (observation time, TAI93)";
            p_info->description = strdup(tmpStr.c_str());
            free(p_info->units);
            p_info->units = strdup("counts");
            free(p_info->palette);
            p_info->palette = strdup(PRODUCT_DEFAULT_palette);
            free(p_info->dataType);
            p_info->dataType = strdup("float");
            tmpStr = p_info->suffix;
            free(p_info->suffix);
            tmpStr += "_obs_time";
            p_info->suffix = strdup(tmpStr.c_str());
            p_info->fillValue = PRODUCT_DEFAULT_fillValue;
            p_info->validMin = PRODUCT_DEFAULT_validMin;
            p_info->validMax = PRODUCT_DEFAULT_validMin;
            free(p_info->displayScale);
            p_info->displayScale = strdup("linear");
            p_info->displayMin = PRODUCT_DEFAULT_displayMin;
            p_info->displayMax = PRODUCT_DEFAULT_displayMax;
            p_info->addOffset = 0.0;
            p_info->scaleFactor = 1.0;
            break;
        case BinNum:
            tmpStr = p_info->description;
            free(p_info->description);
            tmpStr += " (bin ID number)";
            p_info->description = strdup(tmpStr.c_str());
            free(p_info->units);
            p_info->units = strdup("dimensionless");
            free(p_info->palette);
            p_info->palette = strdup(PRODUCT_DEFAULT_palette);
            free(p_info->dataType);
            p_info->dataType = strdup("int");
            tmpStr = p_info->suffix;
            free(p_info->suffix);
            tmpStr += "_bin_num";
            p_info->suffix = strdup(tmpStr.c_str());
            p_info->fillValue = PRODUCT_DEFAULT_fillValue;
            p_info->validMin = PRODUCT_DEFAULT_validMin;
            p_info->validMax = PRODUCT_DEFAULT_validMin;
            free(p_info->displayScale);
            p_info->displayScale = strdup("linear");
            p_info->displayMin = PRODUCT_DEFAULT_displayMin;
            p_info->displayMax = PRODUCT_DEFAULT_displayMax;
            p_info->addOffset = 0.0;
            p_info->scaleFactor = 1.0;
            break;
    } // switch
    
    if(clo_getBool(optionList, "apply_pal")) {
        clo_option_t* option = clo_findOption(optionList, "palfile");
        if(clo_isOptionSet(option)) {
            outFile->setPalette(clo_getOptionString(option));
            if(outFile2)
                outFile2->setPalette(clo_getOptionString(option));
        } else {
            outFile->setPalette(p_info->palette);
            if(outFile2)
                outFile2->setPalette(p_info->palette);
        }
    }
    
    // set the default scale factors for RGB
    if(doingRGB) {
        p_info->displayMin = 0.01;
        p_info->displayMax = 0.9;
        if(p_info->displayScale)
            free(p_info->displayScale);
        p_info->displayScale = strdup("log");
    }

    // override default scale parameters if set on command line
    if(clo_isSet(optionList, "datamin")) {
        p_info->displayMin = clo_getFloat(optionList, "datamin");
    }
    if(clo_isSet(optionList, "datamax")) {
        p_info->displayMax = clo_getFloat(optionList, "datamax");
    }
    if(clo_isSet(optionList, "scale_type")) {
        if(p_info->displayScale)
            free(p_info->displayScale);
        p_info->displayScale = strdup(clo_getString(optionList, "scale_type"));
    }
    
    outFile->addProduct(p_info);
    if(outFile2)
        outFile2->addProduct(p_info);
    
    freeProductInfo(p_info);    
}

void writeRawFile(L3File* l3File, OutFile* outFile, OutFile* outFile2) {
    L3Bin* l3Bin;
    L3Row* l3Row;
    imageHeight = l3File->getNumRows();
    double resolution = EARTH_CIRCUMFERENCE / (imageHeight*2);
    imageWidth = imageHeight*2;
    int32_t start;
    int32_t numBins;
    int32_t baseBin;
    int32_t endBin;
    int32_t binNum;
    int32_t row, col;
    int i;
    float percent1 = 0.0;
    float percent2 = 0.0;
    uint32_t numFilledPixels = 0;

    meta_l3bType* metaData;

    metaData = outFile->getMetadata();
    sprintf(metaData->title, "%s%s", metaData->sensor_name, " Level-3 Raw Bin Mapped Image");

    outFile->setMapProjection("Integerized Sinusoidal");
    outFile->setResolution(resolution);
    outFile->setSize(imageWidth, imageHeight);
    if(outFile2) {
        strcpy(outFile2->getMetadata()->title, metaData->title);
        outFile2->setMapProjection("Integerized Sinusoidal");
        outFile2->setResolution(resolution);
        outFile2->setSize(imageWidth, imageHeight);
    }

    // setup all the product structures
    for (size_t i=0; i<prodNameList.size(); i++) {
        setupProduct(prodNameList[i], prodMeasurementList[i], outFile, outFile2);
    }

    printStartInfo(outFile);
    
    if(!outFile->open()) {
        printf("-E- Could not open ofile=\"%s\".\n", outFile->getFileName().c_str());
        exit(EXIT_FAILURE);
    }

    if(outFile2) {
        if(!outFile2->open()) {
            printf("-E- Could not open ofile2=\"%s\".\n", outFile2->getFileName().c_str());
            exit(EXIT_FAILURE);
        }
    }

    float centralMeridian = clo_getFloat(optionList, "central_meridian");
    if(centralMeridian < -900.0)
        centralMeridian = 0;
    i = 0;
    while(centralMeridian < 0.0) {
        centralMeridian += 360.0;
        i++;
        if(i>5) {
            printf("-E- central meridian is way off\n");
            exit(EXIT_FAILURE);
        }
    }
    i = 0;
    while(centralMeridian >= 360.0) {
        centralMeridian -= 360.0;
        i++;
        if(i>5) {
            printf("-E- central meridian is way off\n");
            exit(EXIT_FAILURE);
        }
    }

    for(row=imageHeight-1; row>=0; row--){
        if(want_verbose) {
            percent2 = (imageHeight - row) / (float)imageHeight;
            if(percent2-percent1 > percentDelta) {
                percent1 = percent2;
                printf("\r%2d%% complete", (int)(percent2*100));
                fflush(stdout);
            }
        }

        l3Row = l3File->getRow(row);
        numBins = l3File->getShape()->getNumCols(row);
        baseBin = l3File->getShape()->getBaseBin(row);
        endBin = baseBin + numBins;
        binNum = baseBin + numBins * centralMeridian / 360.0;
        start = (imageWidth - numBins) / 2;

        // clear out beginning empty pixels
        for(col=0; col<start; col++) {
            outFile->fillPixel(col);
            if(outFile2)
                outFile2->fillPixel(col);
        }

        // set pixel values
        for(int i=0; i<numBins; i++) {
            l3Bin = l3Row->getBin(binNum);
            if(l3Bin) {
                if(doingRGB) {
                    outFile->setPixelRGB(col, l3Bin->getMean(0), l3Bin->getMean(1), l3Bin->getMean(2));
                    if(outFile2)
                        outFile2->setPixelRGB(col, l3Bin->getMean(0), l3Bin->getMean(1), l3Bin->getMean(2));
                } else {
                    for(size_t prod=0; prod<prodNameList.size(); prod++) {
                        float val;
                        switch(prodMeasurementList[prod]) {
                            case Avg:
                                val = l3Bin->getMean(prod);
                                break;
                            case Stdev:
                                val = l3Bin->getStdev(prod);
                                break;
                            case Variance:
                                val = l3Bin->getVariance(prod);
                                break;
                            case Nobs:
                                val = l3Bin->getNobs();
                                break;
                            case Nscenes:
                                val = l3Bin->getNscenes();
                                break;
                            case ObsTime:
                                val = l3Bin->getObsTime();
                                break;
                            case BinNum:
                                val = l3Bin->getBinNum();
                                break;
                            default:
                                val = l3Bin->getMean(prod);
                        }
                        outFile->setPixel(col, val, prod);
                        if(outFile2)
                            outFile2->setPixel(col, val, prod);
                    }
                }
                numFilledPixels++;
            } else {
                outFile->missingPixel(col);
                if(outFile2)
                    outFile2->missingPixel(col);
            }
            col++;
            binNum++;
            if(binNum >= endBin)
                binNum = baseBin;
        }

        // clear out trailing empty pixels
        for(; col<imageWidth; col++) {
            outFile->fillPixel(col);
            if(outFile2)
                outFile2->fillPixel(col);
        }

        outFile->writeLine();
        if(outFile2)
            outFile2->writeLine();
    }

    outFile->setNumFilledPixels(numFilledPixels);
    outFile->close();
    if(outFile2) {
        outFile2->setNumFilledPixels(numFilledPixels);
        outFile2->close();
    }

}

void writeSmiFile(L3File* l3File, OutFile* outFile, OutFile* outFile2) 
{
    float percent1, percent2;
    int32_t numFilledPixels = 0;

    meta_l3bType* metaData;

    metaData = outFile->getMetadata();
    sprintf(metaData->title, "%s%s", metaData->sensor_name, " Level-3 Standard Mapped Image");
        
    float centralMeridian = clo_getFloat(optionList, "central_meridian");
    if(centralMeridian > -900.0) {
        metaData->west = centralMeridian - 180;
        metaData->east = centralMeridian + 180;
    } else {
        centralMeridian = 0;
    }
    
    double heightInDeg = metaData->north - metaData->south;
    double widthInDeg = metaData->east - metaData->west;
    double resolution = outFile->getResolution();
 
    // set up image parameters
    imageWidth = widthInDeg / 360.0 * EARTH_CIRCUMFERENCE / resolution;
    imageHeight = heightInDeg / 360.0 * EARTH_CIRCUMFERENCE / resolution;
    double deltaLon = widthInDeg / imageWidth;
    double deltaLat = heightInDeg / imageHeight;

    outFile->setSize(imageWidth, imageHeight);
    outFile->setMapProjection("Equidistant Cylindrical");
    if(outFile2) {
        strcpy(outFile2->getMetadata()->title, metaData->title);
        outFile2->getMetadata()->east = metaData->east;
        outFile2->getMetadata()->west = metaData->west;
        outFile2->setSize(imageWidth, imageHeight);
        outFile2->setMapProjection("Equidistant Cylindrical");        
    }
    
    // setup all the product structures
    for (size_t i=0; i<prodNameList.size(); i++) {
        setupProduct(prodNameList[i], prodMeasurementList[i], outFile, outFile2);
    }

    // set up quality processing
    setupQualityProcessing(l3File, outFile, outFile2);
    
    printStartInfo(outFile);

    if (!outFile->open()) {
        printf("-E- Could not open ofile=\"%s\".\n", outFile->getFileName().c_str());
        exit(1);
    }

    if(outFile2) {
        if (!outFile2->open()) {
            printf("-E- Could not open ofile2=\"%s\".\n", outFile2->getFileName().c_str());
            exit(1);
        }
    }
    
    InterpType interp = interpStr2Type(clo_getString(optionList, "interp"));
    
    // loop through output pixels
    double lat = metaData->north - (deltaLat / 2.0);
    double lon;
    L3Bin* l3Bin;
    
    percent1 = 0;
    for (int row = 0; row < imageHeight; row++) {
        if (want_verbose) {
            percent2 = (float) row / (float) imageHeight;
            if (percent2 - percent1 > percentDelta) {
                percent1 = percent2;
                printf("\r%2d%% complete", (int)(percent2*100));
                fflush(stdout);
            }
        }

        lon = metaData->west + (deltaLon / 2.0);
        for (int col = 0; col < imageWidth; col++) {
            switch(interp) {
                case Interp_Nearest:
                    l3Bin = l3File->getClosestBin(lat, lon);
                    break;
                case Interp_Bin: {
                    Point_t pMin(lon - (deltaLon / 2.0), lat - (deltaLat / 2.0));
                    Point_t pMax(lon + (deltaLon / 2.0), lat + (deltaLat / 2.0));
                    Box_t box(pMin, pMax);
                    l3Bin = l3File->getBinsInside(box);
                    break;
                }
                case Interp_Area: {
                    Point_t pMin(lon - (deltaLon / 2.0), lat - (deltaLat / 2.0));
                    Point_t pMax(lon + (deltaLon / 2.0), lat + (deltaLat / 2.0));
                    Box_t box(pMin, pMax);
                    l3Bin = l3File->getBinsInside(box, true);
                    break;
                }
                default:
                    printf("-E- interp = %s is not implemented.", interpType2Str(interp));
                    exit(EXIT_FAILURE);
            }
            if (l3Bin) {
                numFilledPixels++;
                if(doingRGB) {
                    outFile->setPixelRGB(col, l3Bin->getMean(0), l3Bin->getMean(1), l3Bin->getMean(2));
                    if(outFile2)
                        outFile2->setPixelRGB(col, l3Bin->getMean(0), l3Bin->getMean(1), l3Bin->getMean(2));
                } else {
                    for(size_t prod=0; prod<prodNameList.size(); prod++) {
                        float val;
                        switch(prodMeasurementList[prod]) {
                            case Avg:
                                val = l3Bin->getMean(prod);
                                break;
                            case Stdev:
                                val = l3Bin->getStdev(prod);
                                break;
                            case Variance:
                                val = l3Bin->getVariance(prod);
                                break;
                            case Nobs:
                                val = l3Bin->getNobs();
                                break;
                            case Nscenes:
                                val = l3Bin->getNscenes();
                                break;
                            case ObsTime:
                                val = l3Bin->getObsTime();
                                break;
                            case BinNum:
                                val = l3Bin->getBinNum();
                                break;
                            default:
                                val = l3Bin->getMean(prod);
                        }
                        outFile->setPixel(col, val, prod);
                        if(outFile2)
                            outFile2->setPixel(col, val, prod);
                    }
                }
                if(doingQuality) {
                    outFile->setQuality(col, l3Bin->getQuality());
                    if(outFile2)
                        outFile2->setQuality(col, l3Bin->getQuality());
                }
            } else {
                outFile->missingPixel(col);
                if(outFile2)
                    outFile2->missingPixel(col);
            }
            lon += deltaLon;
        } // for col
        outFile->writeLine();
        if(outFile2)
            outFile2->writeLine();
        lat -= deltaLat;
    } // for row
    outFile->setNumFilledPixels(numFilledPixels);
    outFile->close();
    if(outFile2) {
        outFile2->setNumFilledPixels(numFilledPixels);
        outFile2->close();
    }
}

void writeProj4File(L3File* l3File, char* projectionStr, 
        OutFile* outFile, OutFile* outFile2) {
    string projStr;
    float percent1, percent2;
    int32_t numFilledPixels = 0;
    meta_l3bType* metaData;

    // how about circumference of earth + 25%
    double limitMin = EARTH_CIRCUMFERENCE * -0.625;
    double limitMax = EARTH_CIRCUMFERENCE * 0.625;

    metaData = outFile->getMetadata();
    float centralMeridian = clo_getFloat(optionList, "central_meridian");
    double heightInDeg = metaData->north - metaData->south;
    double widthInDeg = metaData->east - metaData->west;
    double resolution = outFile->getResolution();

    if(strcasecmp(projectionStr, "mollweide") == 0) {
        // set the title
        sprintf(metaData->title, "%s%s", metaData->sensor_name, " Level-3 Mollweide Mapped Image");

        if(centralMeridian > -900.0) {
            metaData->west = centralMeridian - 180;
            metaData->east = centralMeridian + 180;
        } else {
            centralMeridian = (metaData->west + metaData->east) / 2.0;
        }

        char s[50];
        sprintf(s, "%f", centralMeridian);
        projStr = "+proj=moll +lon_0=";
        projStr += s;

        outFile->setMapProjection("Mollweide");
        if(outFile2) {
            strcpy(outFile2->getMetadata()->title, metaData->title);
            outFile2->getMetadata()->east = metaData->east;
            outFile2->getMetadata()->west = metaData->west;
            outFile2->setMapProjection("Mollweide");
        }
            
    } else if(strcasecmp(projectionStr, "lambert") == 0) {
        // set the title
        sprintf(metaData->title, "%s%s", metaData->sensor_name, " Level-3 Lambert Mapped Image");

        if(centralMeridian > -900.0) {
            metaData->west = centralMeridian - 180;
            metaData->east = centralMeridian + 180;
        } else {
            centralMeridian = (metaData->west + metaData->east) / 2.0;
        }

        char s[50];
        sprintf(s, "%f", centralMeridian);
        projStr = "+proj=lcc +lon_0=";
        projStr += s;

        outFile->setMapProjection("Lambert");
        if(outFile2) {
            strcpy(outFile2->getMetadata()->title, metaData->title);
            outFile2->getMetadata()->east = metaData->east;
            outFile2->getMetadata()->west = metaData->west;
            outFile2->setMapProjection("Lambert");
        }

    } else if(strcasecmp(projectionStr, "albersconic") == 0) {
        // set the title
        sprintf(metaData->title, "%s%s", metaData->sensor_name, " Level-3 Albers Equal-Area Conic Mapped Image");

        if(centralMeridian > -900.0) {
            metaData->west = centralMeridian - 180;
            metaData->east = centralMeridian + 180;
        } else {
            centralMeridian = (metaData->west + metaData->east) / 2.0;
        }

        char s[50];
        sprintf(s, "%f", centralMeridian);
        projStr = "+proj=aea +lon_0=";
        projStr += s;

        outFile->setMapProjection("Albersconic");
        if(outFile2) {
            strcpy(outFile2->getMetadata()->title, metaData->title);
            outFile2->getMetadata()->east = metaData->east;
            outFile2->getMetadata()->west = metaData->west;
            outFile2->setMapProjection("Albersconic");
        }
        
    } else if(strcasecmp(projectionStr, "mercator") == 0) {
        // set the title
        sprintf(metaData->title, "%s%s", metaData->sensor_name, " Level-3 Mercator Mapped Image");

        if(centralMeridian > -900.0) {
            metaData->west = centralMeridian - 180;
            metaData->east = centralMeridian + 180;
        } else {
            centralMeridian = (metaData->west + metaData->east) / 2.0;
        }

        char s[50];
        sprintf(s, "%f", centralMeridian);
        projStr = "+proj=merc +lon_0=";
        projStr += s;

        outFile->setMapProjection("Mercator");
        if(outFile2) {
            strcpy(outFile2->getMetadata()->title, metaData->title);
            outFile2->getMetadata()->east = metaData->east;
            outFile2->getMetadata()->west = metaData->west;
            outFile2->setMapProjection("Mercator");
        }
        
    } else if(strcasecmp(projectionStr, "ease2") == 0) {
        // set the title
        sprintf(metaData->title, "%s%s", metaData->sensor_name, " Level-3 Ease Grid 2 Mapped Image");

        if(centralMeridian > -900.0) {
            metaData->west = centralMeridian - 180;
            metaData->east = centralMeridian + 180;
        } else {
            centralMeridian = (metaData->west + metaData->east) / 2.0;
        }

        char s[50];
        sprintf(s, "%f", centralMeridian);
        projStr = "+proj=cea +lat_0=0 +lat_ts=30 +ellps=WGS84 +datum=WGS84 +units=m +lon_0=";
        projStr += s;

        outFile->setMapProjection("Ease2");
        if(outFile2) {
            strcpy(outFile2->getMetadata()->title, metaData->title);
            outFile2->getMetadata()->east = metaData->east;
            outFile2->getMetadata()->west = metaData->west;
            outFile2->setMapProjection("Ease2");
        }
        
    } else {
        // set the title
        sprintf(metaData->title, "%s%s", metaData->sensor_name, " Level-3 Proj4 Mapped Image");

        projStr = projectionStr;
        outFile->setMapProjection(projectionStr);
        if(outFile2) {
            strcpy(outFile2->getMetadata()->title, metaData->title);
            outFile2->setMapProjection(projectionStr);
        }
        
    }

    projPJ pj_new = pj_init_plus(projStr.c_str());
    projPJ pj_latlong = pj_latlong_from_proj(pj_new);

    // calculate the min and max
    double minX = limitMax;
    double minY = limitMax;
    double maxX = limitMin;
    double maxY = limitMin;
    double *tmpX, *tmpY;
    double lat, lon;
    double deltaLat = (metaData->north - metaData->south) / (NUM_SEARCH_POINTS-1);
    double deltaLon = (metaData->east - metaData->west) / (NUM_SEARCH_POINTS-1);

    tmpX = (double*) allocateMemory(NUM_SEARCH_POINTS*sizeof(double), "tmpX");
    tmpY = (double*) allocateMemory(NUM_SEARCH_POINTS*sizeof(double), "tmpY");

    lat = metaData->south;
    for(int j=0; j<NUM_SEARCH_POINTS; j++) {
        lon = metaData->west;
        for(int i=0; i<NUM_SEARCH_POINTS; i++) {
            tmpX[i] = lon * DEG_TO_RAD;       // convert to radians
            tmpY[i] = lat * DEG_TO_RAD;
            lon += deltaLon;
        }
        if(pj_transform(pj_latlong, pj_new, NUM_SEARCH_POINTS, 1, tmpX, tmpY, NULL )) {
            printf("Error - min/max proj4 transformation blew up\n");
            exit(1);
        }

        for(int i=0; i<NUM_SEARCH_POINTS; i++) {
            if(isnormal(tmpX[i]) && isnormal(tmpY[i])) {
                if(tmpX[i] < limitMax && tmpX[i] > limitMin) {
                    if(tmpX[i] < minX)
                        minX = tmpX[i];
                    if(tmpX[i] > maxX)
                        maxX = tmpX[i];
                }
                if(tmpY[i] < limitMax && tmpY[i] > limitMin) {
                    if(tmpY[i] < minY)
                        minY = tmpY[i];
                    if(tmpY[i] > maxY)
                        maxY = tmpY[i];
                }
            }
        }

        lat += deltaLat;
    }
    free(tmpX);
    free(tmpY);

    double startX = minX + resolution/2;
    double startY = maxY - resolution/2;
    imageWidth = (maxX-minX) / resolution;
    imageHeight = (maxY-minY) / resolution;
    double x, y;

    outFile->setSize(imageWidth, imageHeight);
    if(outFile2)
        outFile2->setSize(imageWidth, imageHeight);

    // setup all the product structures
    for (size_t i=0; i<prodNameList.size(); i++) {
        setupProduct(prodNameList[i], prodMeasurementList[i], outFile, outFile2);
    }
    
    // set up quality processing
    setupQualityProcessing(l3File, outFile, outFile2);

    printStartInfo(outFile);
    
    if(!outFile->open()) {
        printf("-E- Could not open ofile=\"%s\".\n", outFile->getFileName().c_str());
        exit(1);
    }

    if(outFile2) {
        if(!outFile2->open()) {
            printf("-E- Could not open ofile2=\"%s\".\n", outFile2->getFileName().c_str());
            exit(1);
        }
    }
    
    tmpX = (double*)allocateMemory(imageWidth*sizeof(double), "tmpX");
    tmpY = (double*)allocateMemory(imageWidth*sizeof(double), "tmpY");

    y = startY;
    percent1 = 0;
    for(int j=0; j<imageHeight; j++) {
        if(want_verbose) {
            percent2 = (float)j / (float)imageHeight;
            if(percent2-percent1 > percentDelta) {
                percent1 = percent2;
                printf("\r%2d%% complete", (int)(percent2*100));
                fflush(stdout);
            }
        }
        x = startX;
        for(int i=0; i<imageWidth; i++) {
            tmpX[i] = x;
            tmpY[i] = y;
            x += resolution;
        }
        if(pj_transform(pj_new, pj_latlong, imageWidth, 1, tmpX, tmpY, NULL )) {
            printf("Error - min/max proj4 transformation blew up\n");
            exit(1);
        }

        // make them all deg
        for(int i=0; i<imageWidth; i++) {
            tmpX[i] *= RAD_TO_DEG;
            tmpY[i] *= RAD_TO_DEG;

            if(tmpX[i] < metaData->west)
                tmpX[i] += 360.0;
        }

        int startPoint;
        int endPoint;
        double mid = tmpX[imageWidth/2];
        // find starting point in line
        for(startPoint=imageWidth/2; startPoint>=0; startPoint--) {
            lat = tmpX[startPoint];
            if(lat > mid)
                lat -= 360;
            if(lat < metaData->west)
                break;
        }
        startPoint++;

        // find ending point in line
        for(endPoint=imageWidth/2; endPoint<imageWidth; endPoint++) {
            lat = tmpX[endPoint];
            if(lat < mid)
                lat += 360;
            if(lat > metaData->east)
                break;
        }
        endPoint--;

        for(int i=0; i<imageWidth; i++) {
            if(i<startPoint || i>endPoint) {
                outFile->fillPixel(i);
                if(outFile2)
                    outFile2->fillPixel(i);
            } else if(!isnormal(tmpX[i]) || !isnormal(tmpY[i])) {
                outFile->fillPixel(i);
                if(outFile2)
                    outFile2->fillPixel(i);
            } else if(tmpY[i] > metaData->north || tmpY[i] < metaData->south ||
                    tmpX[i] > metaData->east || tmpX[i] < metaData->west) {
                outFile->fillPixel(i);
                if(outFile2)
                    outFile2->fillPixel(i);
            } else {
                L3Bin* l3Bin = l3File->getClosestBin(tmpY[i], tmpX[i]);
                if(l3Bin) {
                    numFilledPixels++;
                    if(doingRGB) {
                        outFile->setPixelRGB(i, l3Bin->getMean(0), l3Bin->getMean(1), l3Bin->getMean(2));
                        if(outFile2)
                            outFile2->setPixelRGB(i, l3Bin->getMean(0), l3Bin->getMean(1), l3Bin->getMean(2));
                    } else {
                        for(size_t prod=0; prod<prodNameList.size(); prod++) {
                            float val;
                            switch(prodMeasurementList[prod]) {
                                case Avg:
                                    val = l3Bin->getMean(prod);
                                    break;
                                case Stdev:
                                    val = l3Bin->getStdev(prod);
                                    break;
                                case Variance:
                                    val = l3Bin->getVariance(prod);
                                    break;
                                case Nobs:
                                    val = l3Bin->getNobs();
                                    break;
                                case Nscenes:
                                    val = l3Bin->getNscenes();
                                    break;
                                case ObsTime:
                                    val = l3Bin->getObsTime();
                                    break;
                                case BinNum:
                                    val = l3Bin->getBinNum();
                                    break;
                                default:
                                    val = l3Bin->getMean(prod);
                            }
                            outFile->setPixel(i, val, prod);
                            if(outFile2)
                                outFile2->setPixel(i, val, prod);
                        }
                    }
                } else
                    outFile->missingPixel(i);
                    if(outFile2)
                        outFile2->missingPixel(i);
            }
        } // for width
        outFile->writeLine();
        if(outFile2)
            outFile2->writeLine();
        y -= resolution;
    } // for height
    free(tmpX);
    free(tmpY);
    pj_free(pj_new);
    pj_free(pj_latlong);
    outFile->setNumFilledPixels(numFilledPixels);
    outFile->close();
    if(outFile2) {
        outFile2->setNumFilledPixels(numFilledPixels);
        outFile2->close();
    }
}

OutFile* makeOutputFile(const char* oformatStr, bool useColor, bool useRGB) {
    OutFile* outFile = NULL;

    const char* oformatStr2 = getFileFormatName(oformatStr);

    if(oformatStr2==NULL) {
        printf("-E- Unknown output file format \"%s\"\n", oformatStr);
        exit(EXIT_FAILURE);
    }

    string oformat = oformatStr2;
    if(oformat.compare("PPM") == 0) {
        if(useRGB) {
            outFile = new OutFile_ppm_rgb();
        } else {
            if(useColor) {
                outFile = new OutFile_ppm();
            } else {
                outFile = new OutFile_pgm();
            }
        }
    } else if(oformat.compare("PNG") == 0) {
        if(useRGB) {
            outFile = new OutFile_png_rgb();
        } else {
            outFile = new OutFile_png(useColor);
        }
    } else if(oformat.compare("TIFF") == 0) {
        if(useRGB) {
            outFile = new OutFile_tiff_rgb();
        } else {
            if(useColor) {
                outFile = new OutFile_tiff_color();
            } else {
                outFile = new OutFile_tiff_gray();
            }
        }
    } else if(oformat.compare("HDF4") == 0) {
    	outFile = new OutFile_hdf4();
    } else if(oformat.compare("netCDF4") == 0) {
        outFile = new OutFile_netcdf4();
    } else {
        printf("-E- Output file type %s not implemented\n", oformat.c_str());
        exit(EXIT_FAILURE);
    }
    return outFile;
}



//-------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    OutFile* outFile;
    OutFile* outFile2 = NULL;
    char* ifileName;
    char* ofileName;
    string oformat;
    int status;
    bool wroteFile = false;
    int i;
    char* tmpStr;
    
    char softwareVersion[200];
    sprintf(softwareVersion, "%d.%d.%d-r%d", L3MAPGEN_VERSION_MAJOR, L3MAPGEN_VERSION_MINOR,
            L3MAPGEN_VERSION_PATCH_LEVEL, SVN_REVISION);

    optionList = clo_createList();

    l3mapgen_init_options(optionList, softwareVersion);
    if(argc == 1) {
        clo_printUsage(optionList);
        exit(1);
    }
    l3mapgen_read_options(optionList, argc, argv);

    if(clo_getBool(optionList, "quiet")) {
        want_verbose = 0;
    }

    ifileName = clo_getString(optionList, "ifile");
    ofileName = clo_getString(optionList, "ofile");

    // try SMI input file
    int oldVerbose = want_verbose;
    want_verbose = 0;
    L3File* l3File = new L3FileSMI();
    if(!l3File->open(ifileName)) {
        
        // try real L3 bin format
        delete l3File;
        l3File = new L3File();
        if(!l3File->open(ifileName)) {
            printf("-E- Could not open ifile=\"%s\".\n", ifileName);
            exit(EXIT_FAILURE);
        }
    }
    want_verbose = oldVerbose;

    if(clo_getBool(optionList, "use_rgb")) {
        doingRGB = true;
        prodName = clo_getRawString(optionList, "product_rgb");
    } else {
        doingRGB = false;
        if(clo_isSet(optionList, "product")) {
            prodName = clo_getRawString(optionList, "product");
        } else {
            prodName.clear();
            for(int i=0; i<l3File->getNumProducts(); i++) {
                if(i != 0)
                    prodName += ",";
                prodName += l3File->getProductName(i);
            }
        }
    }
    boost::split(prodNameList, prodName, boost::is_any_of(","));
    string cleanProdName;
    vector<string> parts;
    for(size_t i=0; i<prodNameList.size(); i++) {
        if(i!=0)
            cleanProdName += ",";
        boost::split(parts, prodNameList[i], boost::is_any_of(":"));
        if(parts.size() == 1) {
            cleanProdName += parts[0];
            prodMeasurementList.push_back(Avg);
        } else if(parts.size() == 2) {
            prodNameList[i] = parts[0]; // get rid of the modifier
            cleanProdName += parts[0];
            if(parts[1].compare("avg") == 0)
                prodMeasurementList.push_back(Avg);
            else if(parts[1].compare("stdev") == 0)
                prodMeasurementList.push_back(Stdev);
            else if(parts[1].compare("var") == 0)
                prodMeasurementList.push_back(Variance);
            else if(parts[1].compare("nobs") == 0)
                prodMeasurementList.push_back(Nobs);
            else if(parts[1].compare("nscenes") == 0)
                prodMeasurementList.push_back(Nscenes);
            else if(parts[1].compare("obs_time") == 0)
                prodMeasurementList.push_back(ObsTime);
            else if(parts[1].compare("bin_num") == 0)
                prodMeasurementList.push_back(BinNum);
            else {
                printf("-E- measurement type \"%s\" not understood for product \"%s\".\n", parts[1].c_str(), parts[0].c_str());
                exit(EXIT_FAILURE);
            }
        } else {
            printf("-E- product name not understood \"%s\".\n", prodNameList[i].c_str());
            exit(EXIT_FAILURE);
        }
    }
    
    if(!l3File->setActiveProductList(cleanProdName.c_str())) {
        printf("-E- Could not find product=\"%s\".\n", cleanProdName.c_str());
        exit(EXIT_FAILURE);
    }

    l3File->setFudgeFactor(clo_getFloat(optionList, "fudge"));

    // copy the L3 meta data since we will modify it a bit for the output file.
    meta_l3bType metaData = *l3File->getMetaData();

    outFile = makeOutputFile(clo_getString(optionList, "oformat"), clo_getBool(optionList, "apply_pal"), doingRGB);
    if(clo_isSet(optionList, "ofile2")) {
        outFile2 = makeOutputFile(clo_getString(optionList, "oformat2"), clo_getBool(optionList, "apply_pal"), doingRGB);
    }
        
    // this will be the # of meters across 1 pixel in the center of the scene
    outFile->setResolution(clo_getString(optionList, "resolution"));
    if(outFile2)
        outFile2->setResolution(outFile->getResolution());
    
    // projection
    clo_option_t* projectionOption = clo_findOption(optionList, "projection");
    char* projectionStr = clo_getOptionRawString(projectionOption);

    // check the metadata of the bin file
    if(metaData.north == metaData.south) {
        printf("-E- north and south metadata are equal.\n");
        exit(110);
    }
    if(metaData.east == metaData.west) {
        printf("-E- east and west metadata are equal.\n");
        exit(110);
    }
        
    // default to whole globe for SMI files
    if((strcmp(projectionStr, "smi") == 0) || 
        (strcmp(projectionStr, "raw") == 0)) {
         metaData.north = 90.0;
         metaData.south = -90.0;
         metaData.east = 180.0;
         metaData.west = -180.0;
    }
    // read in north, south, east, west from command line
    float tmpf = clo_getFloat(optionList, "north");
    if(tmpf > -900.0) {
        metaData.north = tmpf;
    }
    tmpf = clo_getFloat(optionList, "south");
    if(tmpf > -900.0) {
        metaData.south = tmpf;
    }
    tmpf = clo_getFloat(optionList, "east");
    if(tmpf > -900.0) {
        metaData.east = tmpf;
    }
    tmpf = clo_getFloat(optionList, "west");
    if(tmpf > -900.0) {
        metaData.west = tmpf;
    }

    if(metaData.north <= metaData.south) {
        printf("-E- north must be greater than south.\n");
        exit(EXIT_FAILURE);
    }
    if(metaData.east <= metaData.west) {
        metaData.east += 360.0;
    }
    if(metaData.east <= metaData.west) {
        printf("-E- east must be greater than west.\n");
        exit(EXIT_FAILURE);
    }
    double heightInDeg = metaData.north - metaData.south;
    if(heightInDeg > 180.0) {
        printf("-E- height in degrees must be less than or equal to 180.\n");
        exit(EXIT_FAILURE);
    }
    double widthInDeg = metaData.east - metaData.west;
    if(widthInDeg > 360.0) {
        printf("-E- width in degrees must be less than or equal to 360.\n");
        exit(EXIT_FAILURE);
    }

    // set other fields in the metadata
    strcpy(metaData.soft_name, "l3mapgen");
    strcpy(metaData.soft_ver, softwareVersion);
    if ((tmpStr = strrchr(ifileName, '/')) != NULL)
        tmpStr++;
    else
        tmpStr = ifileName;
    strcpy(metaData.infiles, tmpStr);
    metaData.proc_con[0] = 0;
    for(i=0; i<argc; i++) {
        strcat(metaData.proc_con, argv[i]);
        strcat(metaData.proc_con, " ");
    }
    strcpy(metaData.pversion, clo_getString(optionList,"pversion"));

    // set input parameters
    metaData.input_parms[0] = 0;
    int numOptions = clo_getNumOptions(optionList);
    for(int i=0; i<numOptions; i++) {
        clo_option_t* option = clo_getOption(optionList, i);
        if(option) {
            if(option->key[0] != '-') {
                char* val = option->valStr;
                if(val == NULL)
                    val = option->defaultVal;
                if(val) {
                    strcat(metaData.input_parms, option->key);
                    strcat(metaData.input_parms, "=");
                    strcat(metaData.input_parms, val);
                    strcat(metaData.input_parms, "|");
                }
            }
        }
    }

    // set processing time
    get_time(metaData.ptime);

    outFile->setMetaData(&metaData);
    outFile->setFileName(ofileName);
    outFile->setDeflate(clo_getInt(optionList, "deflate"));
    if(outFile2) {
        outFile2->setMetaData(&metaData);
        outFile2->setFileName(clo_getString(optionList, "ofile2"));
        outFile2->setDeflate(clo_getInt(optionList, "deflate"));
    }
    
    if(strcasecmp(projectionStr, "raw") == 0) {
        writeRawFile(l3File, outFile, outFile2);
    } else if(strcasecmp(projectionStr, "smi") == 0) {
        writeSmiFile(l3File, outFile, outFile2);
    } else if(strcasecmp(projectionStr, "platecarree") == 0) {
        writeSmiFile(l3File, outFile, outFile2);
    } else {
        writeProj4File(l3File, projectionStr, outFile, outFile2);
    }

    float threshold = clo_getFloat(optionList, "threshold");
    if(threshold > 0.0) {
        if(outFile->getPercentFilledPixels() < threshold) {
            printf("\nPercent filled pixels (%.1f) is below the threshold (%.1f)\n",
                    outFile->getPercentFilledPixels(), threshold);
            printf("Deleting output file.\n");
            
            string cmd = "rm -f ";
            cmd += outFile->getFileName();
            system(cmd.c_str());

            if(outFile2) {
                cmd = "rm -f ";
                cmd += outFile2->getFileName();
                system(cmd.c_str());
            }
            
            exit(EXIT_FAILURE);
        }
    }            
    
    printEndInfo(outFile);
    
    free(outFile);
    if(outFile2)
        free(outFile2);
    free(l3File);

    return EXIT_SUCCESS;
}
