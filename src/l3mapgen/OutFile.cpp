
#include <xtiffio.h>

#include "OutFile.h"

#include <stdio.h>
#include <math.h>
#include <genutils.h>
#include <string>
#include <float.h>

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <geo_tiffp.h>


using namespace std;

//-------------------------------------------------------------------------------------------
// OutFile::ProductStuff
//-------------------------------------------------------------------------------------------
OutFile::ProductStuff::ProductStuff(int32_t width, const productInfo_t* productInfo) {
    this->width = width;
    this->productInfo = allocateProductInfo();
    copyProductInfo(this->productInfo, productInfo);
    dataStorage = UByteDS;
    scaleType = Linear;
    scale = 1.0;
    offset = 0.0;
    minOutputVal = 0;
    maxOutputVal = 255;
    minVal = minOutputVal;
    maxVal = maxOutputVal;
    missingValue = 255;
    lineData = (double*) allocateMemory(width * sizeof(double), "OutFile::ProductStuff::lineData");
}

OutFile::ProductStuff::ProductStuff(const OutFile::ProductStuff& pStuff) {
    width = pStuff.width;
    productInfo = allocateProductInfo();
    copyProductInfo(productInfo, pStuff.productInfo);
    dataStorage = pStuff.dataStorage;
    scaleType = pStuff.scaleType;
    scale = pStuff.scale;
    offset = pStuff.offset;
    minOutputVal = pStuff.minOutputVal;
    maxOutputVal = pStuff.maxOutputVal;
    minVal = pStuff.minVal;
    maxVal = pStuff.maxVal;
    missingValue = pStuff.missingValue;
    lineData = (double*) allocateMemory(width * sizeof(double), "OutFile::ProductStuff::lineData");
}

OutFile::ProductStuff::~ProductStuff() {
    freeProductInfo(productInfo);
    free(lineData);
}

/**
 * set the scale factors.  note that default minOutputVal=0, maxOutputVal=255
 * @param min min geophysical value
 * @param max max geophysical value
 * @param log do you want log10 scaling
 */
void OutFile::ProductStuff::setScale(double min, double max, ScaleType scaleType) {
    this->scaleType = scaleType;
    minVal = min;
    maxVal = max;
    switch(scaleType) {
    case Linear:
        offset = min;
        scale = (max - offset) / (maxOutputVal - minOutputVal);
        break;
    case Log:
        offset = log10(min);
        scale = (log10(max) - offset) / (maxOutputVal - minOutputVal);
        break;
    case ArcTan:
        offset = min;
        scale = max;
        minVal = calcPhysicalVal(minOutputVal);
        maxVal = calcPhysicalVal(maxOutputVal);
        break;
    default:
        printf("-E- OutFile::setScale - invalid scaleType = %d\n", (int)scaleType);
        exit(EXIT_FAILURE);
    }
}

void OutFile::ProductStuff::setScale(double min, double max, ScaleType scaleType,
        double minOutput, double maxOutput) {
    minOutputVal = minOutput;
    maxOutputVal = maxOutput;
    setScale(min, max, scaleType);
}

/**
 * set the scale factors.  note that default minOutputVal=0, maxOutputVal=255
 * @param scale slope
 * @param offset intercept
 * @param scaleType type of scaling to calculate
 */
void OutFile::ProductStuff::setScaleOffset(double scale, double offset, ScaleType scaleType) {
    this->scaleType = scaleType;
    this->scale = scale;
    this->offset = offset;
	
    // have to set these so calcPhysicalVal does not limit the physical val
    minVal = 0-FLT_MAX;
    maxVal = FLT_MAX;
    minVal = calcPhysicalVal(minOutputVal);
    maxVal = calcPhysicalVal(maxOutputVal);
}

void OutFile::ProductStuff::setScaleOffset(double scale, double offset, ScaleType scaleType,
        double minOutput, double maxOutput) {
    minOutputVal = minOutput;
    maxOutputVal = maxOutput;
    setScaleOffset(scale, offset, scaleType);
}

double OutFile::ProductStuff::calcOutputVal(double val) const {
    double outVal;

    if(val == badPixelValue)
        return missingValue;
        
    // don't scale if out output type is floating point
    if(dataStorage == FloatDS || dataStorage == DoubleDS)
        return val;
    
    switch(scaleType) {
    case Linear:
        outVal = (val - offset) / scale;
        break;
    case Log:
        if(val < 0)
            return minOutputVal;
        else
            outVal = (log10(val) - offset) / scale;
        break;
    case ArcTan:
        outVal = scale * ( (atan(0.5 * val - offset) / atan(offset)) + 1);
        break;
    default:
        printf("-E- OutFile::ProductStuff::calcOutputVal - invalid scaleType = %d\n", (int)scaleType);
        exit(EXIT_FAILURE);
    }

    if(outVal < minOutputVal)
        return minOutputVal;
    if(outVal > maxOutputVal)
        return maxOutputVal;
    return outVal;
}

double OutFile::ProductStuff::calcPhysicalVal(double val) const {
    double physicalVal;

    if(val == missingValue)
        return badPixelValue;
        
    // don't scale if out output type is floating point
    if(dataStorage == FloatDS || dataStorage == DoubleDS)
        return val;
    
    switch(scaleType) {
    case Linear:
        physicalVal = val * scale + offset;
        break;
    case Log:
        physicalVal = pow(10, val * scale + offset);
        break;
    case ArcTan:
        physicalVal = (tan((val / scale - 1) * atan(offset)) + offset) / 0.5;
        break;
    default:
        printf("-E- OutFile::calcPhysicalVal - invalid scaleType = %d\n", (int)scaleType);
        exit(EXIT_FAILURE);
    }

    if(physicalVal < minVal)
        return minVal;
    if(physicalVal > maxVal)
        return maxVal;
    return physicalVal;
}

//-------------------------------------------------------------------------------------------
// OutFile
//-------------------------------------------------------------------------------------------

OutFile::OutFile() {
    width = 0;
    height = 0;
    qualityData = NULL;
    currentLine = 0;
    colorType = Grayscale;
    fileMinVal = DBL_MAX;
    fileMaxVal = 0-DBL_MAX;
    resolution = 0;
    deflate = 0;

    red = (uint8_t*) allocateMemory(256, "red");
    green = (uint8_t*) allocateMemory(256, "green");
    blue = (uint8_t*) allocateMemory(256, "blue");

    // init to grayscale
    for(int i=0; i<256; i++) {
        red[i] = i;
        green[i] = i;
        blue[i] = i;
    }

    metaData = (meta_l3bType*) allocateMemory(sizeof(meta_l3bType), "OutFile::metaData");
    metaData->north = 90.0;
    metaData->south = -90.0;
    metaData->east = 180.0;
    metaData->west = -180.0;

    mapProjection = "Undefined";
}

OutFile::~OutFile() {
    free(red);
    free(green);
    free(blue);
    for(size_t i=0; i<productStuff.size(); i++) {
        delete productStuff[i];
    }
    productStuff.clear();
    free(metaData);
    if(qualityData)
        free(qualityData);
}

void OutFile::setSize(int32_t width, int32_t height) {
    this->width = width;
    this->height = height;
    if(qualityData) {
        free(qualityData);
        qualityData = (uint8_t*) allocateMemory(width, "OutFile::qualityData");
    }
    currentLine = 0;

    for(size_t i=0; i<productStuff.size(); i++) {
        delete productStuff[i];
    }
    productStuff.clear();
}

int32_t OutFile::getWidth() const {
    return width;
}

int32_t OutFile::getHeight() const {
    return height;
}

void OutFile::setFileName(string fileName) {
    this->fileName = fileName;
}

void OutFile::setPixel(int32_t x, double val, int32_t prod) {
    if(x<0 || x>=width) {
        fprintf(stderr, "-E- OutFile::setPixel - x=%d is not within range, width=%d.\n",
                x, width);
        exit(EXIT_FAILURE);
    }
    productStuff[prod]->lineData[x] = val;
    if(val > fileMaxVal)
        fileMaxVal = val;
    if(val < fileMinVal)
        fileMinVal = val;
}

void OutFile::setPixelRGB(int32_t x, float red, float green, float blue) {
    fprintf(stderr, "-E- OutFile::setPixelRGB - RGB not implemented with this file type.\n");
    exit(EXIT_FAILURE);
}

void OutFile::setQuality(int32_t x, uint8_t val) {
    if(x<0 || x>=width) {
        fprintf(stderr, "-E- OutFile::setQuality - x=%d is not within range, width=%d.\n",
                x, width);
        exit(EXIT_FAILURE);
    }
    if(!qualityData) {
        fprintf(stderr, "-E- OutFile::setQuality - qualityData id NULL.\n");
        exit(EXIT_FAILURE);
    }
    qualityData[x] = val;
}

void OutFile::fillPixel(int32_t x) {
    if(x<0 || x>=width) {
        fprintf(stderr, "-E- OutFile::fillPixel - x=%d is not within range, width=%d.\n",
                x, width);
        exit(EXIT_FAILURE);
    }
    for(size_t prod=0; prod<productStuff.size(); prod++)
        productStuff[prod]->lineData[x] = badPixelValue;
    if(qualityData)
        qualityData[x] = qualityUnused;
}

void OutFile::missingPixel(int32_t x) {
    if(x<0 || x>=width) {
        fprintf(stderr, "-E- OutFile::missingPixel - x=%d is not within range, width=%d.\n",
                x, width);
        exit(EXIT_FAILURE);
    }
    for(size_t prod=0; prod<productStuff.size(); prod++)
        productStuff[prod]->lineData[x] = badPixelValue;
    if(qualityData)
        qualityData[x] = qualityUnused;
}

bool OutFile::setPalette(const char* paletteName) {
    char* dataRoot;
    string paletteFileName;
    short r[256], g[256], b[256];

    if ((dataRoot = getenv("OCDATAROOT")) == NULL) {
      printf("OCDATAROOT environment variable is not defined.\n");
      return(EXIT_FAILURE);
    }
    paletteFileName = dataRoot;
    paletteFileName += "/common/palette/";
    paletteFileName += paletteName;
    paletteFileName += ".pal";

    if (getlut_file((char*)paletteFileName.c_str(), r, g, b)) {
        fprintf(stderr, "Error reading palette file %s\n", paletteFileName.c_str());
        return false;
    }
//    if (input.mask) {
//        r[252] = 128;
//        g[252] = 128;
//        b[252] = 128;
//        r[253] = 160;
//        g[253] = 82;
//        b[253] = 45;
//        r[254] = 255;
//        g[254] = 255;
//        b[254] = 255;
//        r[255] = 0;
//        g[255] = 0;
//        b[255] = 0;
//    }
    for (int i = 0; i < 256; i++) {
        red[i] = r[i];
        green[i] = g[i];
        blue[i] = b[i];
    }

    return true;
}

void OutFile::setMetaData(meta_l3bType* metaData) {
    *this->metaData = *metaData;
}

/**
 * Add a product for display type output files
 * @param productInfo info structure to copy
 * @return the index for the new product
 */
int32_t OutFile::addProduct(productInfo_t* productInfo) {
    ProductStuff* stuff = new ProductStuff(width, productInfo);
    
    // setup display scaling
    if(!strcmp(productInfo->displayScale, "linear"))
        stuff->setScale(productInfo->displayMin, productInfo->displayMax, Linear);
    else if(!strcmp(productInfo->displayScale, "log"))
        stuff->setScale(productInfo->displayMin, productInfo->displayMax, Log);
    else if(!strcmp(productInfo->displayScale, "arctan"))
        stuff->setScale(productInfo->displayMin, productInfo->displayMax, ArcTan);
    else {
        printf("-E- OutFile::addProduct - invalid displayScale = %s\n", productInfo->displayScale);
        exit(EXIT_FAILURE);
    }
    
    productStuff.push_back(stuff);
    return productStuff.size() - 1;
}

int32_t OutFile::addProductNonDisplay(productInfo_t* productInfo) {
    ProductStuff* stuff = new ProductStuff(width, productInfo);

    if(!strcmp(productInfo->dataType, "byte")) {
        stuff->dataStorage = ByteDS;
        stuff->minOutputVal = SCHAR_MIN;
        stuff->maxOutputVal = SCHAR_MAX;
    } else if(!strcmp(productInfo->dataType, "ubyte")) {
        stuff->dataStorage = UByteDS;
        stuff->minOutputVal = 0;
        stuff->maxOutputVal = UCHAR_MAX;
    } else if(!strcmp(productInfo->dataType, "short")) {
        stuff->dataStorage = ShortDS;
        stuff->minOutputVal = SHRT_MIN;
        stuff->maxOutputVal = SHRT_MAX;
    } else if(!strcmp(productInfo->dataType, "ushort")) {
        stuff->dataStorage = UShortDS;
        stuff->minOutputVal = 0;
        stuff->maxOutputVal = USHRT_MAX;
    } else if(!strcmp(productInfo->dataType, "int")) {
        stuff->dataStorage = IntDS;
        stuff->minOutputVal = INT_MIN;
        stuff->maxOutputVal = INT_MAX;
    } else if(!strcmp(productInfo->dataType, "uint")) {
        stuff->dataStorage = UIntDS;
        stuff->minOutputVal = 0;
        stuff->maxOutputVal = UINT_MAX;
    } else if(!strcmp(productInfo->dataType, "float")) {
        stuff->dataStorage = FloatDS;
        stuff->minOutputVal = 0;
        stuff->maxOutputVal = 0;
    } else if(!strcmp(productInfo->dataType, "double")) {
        stuff->dataStorage = DoubleDS;
        stuff->minOutputVal = 0;
        stuff->maxOutputVal = 0;
    } else {
        printf("-E- OutFile::addProductNonDisplay - invalid data type = %s\n", productInfo->dataType);
        exit(EXIT_FAILURE);
    }

    // setup scaling
    stuff->setScaleOffset(productInfo->scaleFactor, productInfo->addOffset, Linear);
    stuff->missingValue = productInfo->fillValue;

    productStuff.push_back(stuff);
    return productStuff.size() - 1;
}

void OutFile::setMapProjection(string projection) {
	mapProjection = projection;
}

void OutFile::setNumFilledPixels(int32_t num) {
    if(metaData) {
        metaData->data_bins = num;
    }
}

int32_t OutFile::getNumFilledPixels() {
    if(metaData)
        return metaData->data_bins;
    else
        return -1;
}

float OutFile::getPercentFilledPixels() {
    if(metaData) {
        float numPix = width*height;
        return metaData->data_bins / numPix * 100.0;
    } else
        return -1;
}

void OutFile::resetFileMinMax() {
    fileMinVal = DBL_MAX;
    fileMaxVal = 0-DBL_MAX;
}

void OutFile::setResolution(string resolutionStr) {

    // this will be the # of meters across 1 pixel in the center of the scene
    boost::trim(resolutionStr);
    boost::to_lower(resolutionStr);
    if(resolutionStr.compare("90km") == 0)
        resolution = EARTH_CIRCUMFERENCE / 432.0;
    else if(resolutionStr.compare("36km") == 0)
        resolution = EARTH_CIRCUMFERENCE / 1080.0;
    else if(resolutionStr.compare("18km") == 0)
        resolution = EARTH_CIRCUMFERENCE / 2160.0;
    else if(resolutionStr.compare("9km") == 0)
        resolution = EARTH_CIRCUMFERENCE / 4320.0;
    else if(resolutionStr.compare("4km") == 0)
        resolution = EARTH_CIRCUMFERENCE / 8640.0;
    else if(resolutionStr.compare("2km") == 0)
        resolution = EARTH_CIRCUMFERENCE / 17280.0;
    else if(resolutionStr.compare("1km") == 0)
        resolution = EARTH_CIRCUMFERENCE / 34560.0;
    else if(resolutionStr.compare("hkm") == 0)
        resolution = EARTH_CIRCUMFERENCE / 69120.0;
    else if(resolutionStr.compare("qkm") == 0)
        resolution = EARTH_CIRCUMFERENCE / 138240.0;
    else if(resolutionStr.compare("smi") == 0)
        resolution = EARTH_CIRCUMFERENCE / 4096.0;
    else if(resolutionStr.compare("smi4") == 0)
        resolution = EARTH_CIRCUMFERENCE / 8192.0;
    else if(resolutionStr.compare("land") == 0)
        resolution = EARTH_CIRCUMFERENCE / 8640.0;
    else if(resolutionStr.compare("thirddeg") == 0)
        resolution = EARTH_CIRCUMFERENCE / 1080.0;
    else if(boost::ends_with(resolutionStr, "km")) {
        string val = resolutionStr.substr(0, resolutionStr.length() - 2);
        resolution = atof(val.c_str()) * 1000.0;
    } else if(boost::ends_with(resolutionStr, "deg")) {
        string val = resolutionStr.substr(0, resolutionStr.length() - 3);
        resolution = atof(val.c_str()) / 360.0 * EARTH_CIRCUMFERENCE;
    } else {
        resolution = atof(resolutionStr.c_str());
    }
}
    
void OutFile::setQualityProcessing(bool val) {
    if(val) {
        // if width is not set yet allocate some dummy memory to flag we want
        // to do SST quality processing.
        if(width<=0) {
            if(!qualityData)
                qualityData = (uint8_t*) allocateMemory(2, "OutFile::qualityData");
        } else {
            if(qualityData)
                free(qualityData);
            qualityData = (uint8_t*) allocateMemory(width, "OutFile::qualityData");
        }
    } else {
        if(qualityData) {
            free(qualityData);
            qualityData = NULL;
        }
    }
}

bool OutFile::getQualityProcessing() {
    if(qualityData)
        return true;
    else
        return false;
}

//-------------------------------------------------------------------------------------------
// OutFile_pgm
//-------------------------------------------------------------------------------------------

OutFile_pgm::OutFile_pgm() :
        OutFile() {
    outfp = NULL;
    fileData = NULL;
}

OutFile_pgm::~OutFile_pgm() {
    if(outfp)
        fclose(outfp);
    if(fileData)
        free(fileData);
}

void OutFile_pgm::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if(fileData)
        free(fileData);
    fileData = (uint8_t*) allocateMemory(width, "OutFile_pgm::fileData");
}

bool OutFile_pgm::open() {
    outfp = fopen(fileName.c_str(), "w");
    if (!outfp)
        return false;

    /* Write pgm header */
    fprintf(outfp, "P5\n");
    fprintf(outfp, "%d %d\n", width, height);
    fprintf(outfp, "255\n");

    return true;
}

void OutFile_pgm::writeLine() {
    for(int i=0; i<width; i++)
        fileData[i] = (uint8_t) productStuff[0]->calcOutputVal(productStuff[0]->lineData[i]);
    fwrite(fileData, 1, width, outfp);
    currentLine++;
}

bool OutFile_pgm::close() {
    fclose(outfp);
    outfp = NULL;
    return true;
}

//-------------------------------------------------------------------------------------------
// OutFile_ppm
//-------------------------------------------------------------------------------------------

void OutFile_ppm::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if(fileData)
        free(fileData);
    fileData = (uint8_t*) allocateMemory(width*3, "OutFile_ppm::fileData");
}

bool OutFile_ppm::open() {
    outfp = fopen(fileName.c_str(), "w");
    if (!outfp)
        return false;

    /*
     * Write ppm file header
     */
     fprintf(outfp, "P6\n");
     fprintf(outfp, "%d %d\n", width, height);
     fprintf(outfp, "255\n");

    return true;
}

void OutFile_ppm::writeLine() {
    int j = 0;
    for(int i=0; i<width; i++) {
        uint8_t val = (uint8_t) round(productStuff[0]->calcOutputVal(productStuff[0]->lineData[i]));
        fileData[j++] = red[val];
        fileData[j++] = green[val];
        fileData[j++] = blue[val];
    }
    fwrite(fileData, 1, width*3, outfp);
    currentLine++;
}

//-------------------------------------------------------------------------------------------
// OutFile_ppm_rgb
//-------------------------------------------------------------------------------------------

OutFile_ppm_rgb::OutFile_ppm_rgb() : OutFile() {
    outfp = NULL;
    fileData = NULL;
}

OutFile_ppm_rgb::~OutFile_ppm_rgb() {
    if(outfp)
        fclose(outfp);
    if(fileData)
        free(fileData);
}

void OutFile_ppm_rgb::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if(fileData)
        free(fileData);
    fileData = (uint8_t*) allocateMemory(width*3, "OutFile_ppm_rgb::fileData");
}

bool OutFile_ppm_rgb::open() {
    outfp = fopen(fileName.c_str(), "w");
    if (!outfp)
        return false;
    
    /*
     * Write ppm file header
     */
    fprintf(outfp, "P6\n");
    fprintf(outfp, "%d %d\n", width, height);
    fprintf(outfp, "255\n");
    return true;
}

bool OutFile_ppm_rgb::close() {
    fclose(outfp);
    outfp = NULL;
    return true;
}

void OutFile_ppm_rgb::setPixel(int32_t x, double val, int32_t prod) {
    fprintf(stderr, "-E- OutFile_ppm_rgb::setPixel - only RGB is implemented with this file type.\n");
    exit(EXIT_FAILURE);
}

void OutFile_ppm_rgb::setPixelRGB(int32_t x, float red, float green, float blue) {
    if(x<0 || x>=width) {
        fprintf(stderr, "-E- OutFile_ppm_rgb::setPixel - x=%d is not within range, width=%d.\n",
                x, width);
        exit(EXIT_FAILURE);
    }
        
    uint8_t *ptr = fileData + x*3;

    *ptr = round(productStuff[0]->calcOutputVal(red));
    ptr++;
    *ptr = round(productStuff[0]->calcOutputVal(green));
    ptr++;
    *ptr = round(productStuff[0]->calcOutputVal(blue));
    
    // do this to keep the file min/max reasonable
    if(red > fileMaxVal)
        fileMaxVal = red;
    if(red < fileMinVal)
        fileMinVal = red;
}

void OutFile_ppm_rgb::fillPixel(int32_t x) {
    if(x<0 || x>=width) {
        fprintf(stderr, "-E- OutFile_ppm_rgb::fillPixel - x=%d is not within range, width=%d.\n",
                x, width);
        exit(EXIT_FAILURE);
    }
    uint8_t *ptr = fileData + x*3;
    *ptr = 255;
    ptr++;
    *ptr = 255;
    ptr++;
    *ptr = 255;
}

void OutFile_ppm_rgb::missingPixel(int32_t x) {
    if(x<0 || x>=width) {
        fprintf(stderr, "-E- OutFile_ppm_rgb::missingPixel - x=%d is not within range, width=%d.\n",
                x, width);
        exit(EXIT_FAILURE);
    }
    uint8_t *ptr = fileData + x*3;
    *ptr = 255;
    ptr++;
    *ptr = 255;
    ptr++;
    *ptr = 255;
}

void OutFile_ppm_rgb::writeLine() {
    fwrite(fileData, 1, width*3, outfp);
    currentLine++;
}

//-------------------------------------------------------------------------------------------
// OutFile_png
//-------------------------------------------------------------------------------------------

OutFile_png::OutFile_png(bool color) : OutFile() {
    isColor = color;
    outfp = NULL;
    fileData = NULL;
    info_ptr = NULL;
    png_ptr = NULL;
}

OutFile_png::~OutFile_png() {
    if(outfp)
        fclose(outfp);
    if(fileData)
        free(fileData);
}

void OutFile_png::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if(fileData)
        free(fileData);
    fileData = (uint8_t*) allocateMemory(width, "OutFile_png::fileData");
}

bool OutFile_png::open() {
    outfp = fopen(fileName.c_str(), "w");
    if (!outfp)
        return false;

    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
            NULL, NULL, NULL);
    if (!png_ptr) {
        fprintf(stderr, "-E- Unable to create PNG write structure.\n");
        exit(EXIT_FAILURE);
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
        png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
        fprintf(stderr, "-E- Unable to create PNG info structure.\n");
        exit(EXIT_FAILURE);
    }
    if (setjmp(png_jmpbuf(png_ptr))) {
        png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
        fprintf(stderr, "-E- Unable to call PNG setjmp().\n");
        exit(EXIT_FAILURE);
    }
    png_init_io(png_ptr, outfp);

    if (isColor) {

        // color palette
        png_set_IHDR(png_ptr, info_ptr, width, height,
                8, PNG_COLOR_TYPE_PALETTE, PNG_INTERLACE_NONE,
                PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

        uint8_t pal[256*3];
        for (int i = 0; i < 256; i++) {
            pal[i * 3] = red[i];
            pal[i * 3 + 1] = green[i];
            pal[i * 3 + 2] = blue[i];
        }
        png_set_PLTE(png_ptr, info_ptr, (png_const_colorp) pal, 256);
    } else {

        // Grayscale
        png_set_IHDR(png_ptr, info_ptr, width, height,
                8, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE,
                PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    }
    png_write_info(png_ptr, info_ptr);

    return true;
}

void OutFile_png::writeLine() {
    for(int i=0; i<width; i++)
        fileData[i] = (uint8_t) round(productStuff[0]->calcOutputVal(productStuff[0]->lineData[i]));
    png_write_row(png_ptr, (png_bytep)fileData);
    currentLine++;
}

bool OutFile_png::close() {
    png_write_end(png_ptr, info_ptr);
    png_destroy_write_struct(&png_ptr, &info_ptr);

    fclose(outfp);
    outfp = NULL;
    return true;
}

//-------------------------------------------------------------------------------------------
// OutFile_png_rgb
//-------------------------------------------------------------------------------------------

OutFile_png_rgb::OutFile_png_rgb() : OutFile() {
    outfp = NULL;
    fileData = NULL;
    info_ptr = NULL;
    png_ptr = NULL;
}

OutFile_png_rgb::~OutFile_png_rgb() {
    if(outfp)
        fclose(outfp);
    if(fileData)
        free(fileData);
}

void OutFile_png_rgb::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if(fileData)
        free(fileData);
    fileData = (uint8_t*) allocateMemory(width*3, "OutFile_png_rgb::fileData");
}

bool OutFile_png_rgb::open() {
    outfp = fopen(fileName.c_str(), "w");
    if (!outfp)
        return false;

    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
            NULL, NULL, NULL);
    if (!png_ptr) {
        fprintf(stderr, "-E- Unable to create PNG write structure.\n");
        exit(EXIT_FAILURE);
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
        png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
        fprintf(stderr, "-E- Unable to create PNG info structure.\n");
        exit(EXIT_FAILURE);
    }
    if (setjmp(png_jmpbuf(png_ptr))) {
        png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
        fprintf(stderr, "-E- Unable to call PNG setjmp().\n");
        exit(EXIT_FAILURE);
    }
    png_init_io(png_ptr, outfp);
    png_set_IHDR(png_ptr, info_ptr, width, height,
            8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
            PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    png_write_info(png_ptr, info_ptr);

    return true;
}

bool OutFile_png_rgb::close() {
    if (setjmp(png_jmpbuf(png_ptr))) {
        fprintf(stderr, "-E- OutFile_png_rgb::close - Unable to call PNG setjmp().\n");
        exit(EXIT_FAILURE);
    }
    png_write_end(png_ptr, info_ptr);
    png_destroy_write_struct(&png_ptr, &info_ptr);

    fclose(outfp);
    outfp = NULL;
    return true;
}

void OutFile_png_rgb::setPixel(int32_t x, double val, int32_t prod) {
    fprintf(stderr, "-E- OutFile_png_rgb::setPixel - only RGB is implemented with this file type.\n");
    exit(EXIT_FAILURE);
}

void OutFile_png_rgb::setPixelRGB(int32_t x, float red, float green, float blue) {
    if(x<0 || x>=width) {
        fprintf(stderr, "-E- OutFile_png_rgb::setPixel - x=%d is not within range, width=%d.\n",
                x, width);
        exit(EXIT_FAILURE);
    }
        
    uint8_t *ptr = fileData + x*3;

    *ptr = round(productStuff[0]->calcOutputVal(red));
    ptr++;
    *ptr = round(productStuff[0]->calcOutputVal(green));
    ptr++;
    *ptr = round(productStuff[0]->calcOutputVal(blue));
    
    // do this to keep the file min/max reasonable
    if(red > fileMaxVal)
        fileMaxVal = red;
    if(red < fileMinVal)
        fileMinVal = red;
}

void OutFile_png_rgb::fillPixel(int32_t x) {
    if(x<0 || x>=width) {
        fprintf(stderr, "-E- OutFile_png_rgb::fillPixel - x=%d is not within range, width=%d.\n",
                x, width);
        exit(EXIT_FAILURE);
    }
    uint8_t *ptr = fileData + x*3;
    *ptr = 255;
    ptr++;
    *ptr = 255;
    ptr++;
    *ptr = 255;
}

void OutFile_png_rgb::missingPixel(int32_t x) {
    if(x<0 || x>=width) {
        fprintf(stderr, "-E- OutFile_png_rgb::missingPixel - x=%d is not within range, width=%d.\n",
                x, width);
        exit(EXIT_FAILURE);
    }
    uint8_t *ptr = fileData + x*3;
    *ptr = 255;
    ptr++;
    *ptr = 255;
    ptr++;
    *ptr = 255;
}

void OutFile_png_rgb::writeLine() {

    if (setjmp(png_jmpbuf(png_ptr))) {
        fprintf(stderr, "-E- OutFile_png_rgb::writeLine - Unable to call PNG setjmp().\n");
        exit(EXIT_FAILURE);
    }

    png_write_row(png_ptr, (png_bytep)fileData);
    currentLine++;
}

//-------------------------------------------------------------------------------------------
// OutFile_tiff_color
//-------------------------------------------------------------------------------------------

OutFile_tiff_color::OutFile_tiff_color() : OutFile() {
    fileData = NULL;
    tiff = NULL;
    gtif = NULL;
}

OutFile_tiff_color::~OutFile_tiff_color() {
    if(fileData)
        free(fileData);
}

void OutFile_tiff_color::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if(fileData)
        free(fileData);
    fileData = (uint8_t*) allocateMemory(width, "OutFile_tiff_color::fileData");
}

bool OutFile_tiff_color::open() {
    currentLine = 0;
    tiff = XTIFFOpen(fileName.c_str(), "w");
    if (tiff == NULL) {
        fprintf(stderr, "-E- Could not open outgoing TIFF image\n");
        exit(EXIT_FAILURE);
    }
    gtif = GTIFNew(tiff);
    if (gtif == NULL) {
        fprintf(stderr, "-E- Could not create geoTIFF structure\n");
        exit(EXIT_FAILURE);
    }

     // calc geo TIFF  tags
     double tiepoints[6] = {0, 0, 0, 0, 0, 0};
     double pixscale[3] = {0, 0, 0};

     double north, south, east, west;
     if(metaData) {
         north = metaData->north;
         south = metaData->south;
         east = metaData->east;
         west = metaData->west;
     } else {
         north = 90.0;
         south = -90.0;
         east = 180.0;
         west = -180.0;
     }

     // pixel width
     pixscale[0] = (east - west) / width;

     // pixel height
     pixscale[1] = (north - south) / height;

     // set top left corner pixel lat, lon
     tiepoints[3] = west;
     tiepoints[4] = north;

     TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, width);
     TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, height);
     TIFFSetField(tiff, GTIFF_PIXELSCALE, 3, pixscale);
     TIFFSetField(tiff, GTIFF_TIEPOINTS, 6, tiepoints);

     // Colormap
     uint16_t* rr = (uint16_t*) malloc(256 * sizeof (uint16_t));
     uint16_t* gg = (uint16_t*) malloc(256 * sizeof (uint16_t));
     uint16_t* bb = (uint16_t*) malloc(256 * sizeof (uint16_t));
     if (rr == NULL || gg == NULL || bb == NULL) {
         fprintf(stderr, "-E- Could not allocate memory for TIFF color map\n");
         exit(EXIT_FAILURE);
     }

     // need a colormap of shorts not bytes
     for (int i = 0; i < 256; i++) {
         rr[i] = red[i] << 8;
         gg[i] = green[i] << 8;
         bb[i] = blue[i] << 8;
     }

     TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
     TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_PALETTE);
     TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 8);
     TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 1);
     TIFFSetField(tiff, TIFFTAG_COLORMAP, rr, gg, bb);

     free(rr);
     free(gg);
     free(bb);

     // write geo TIFF keys
     GTIFKeySet(gtif, GTModelTypeGeoKey, TYPE_SHORT, 1, ModelGeographic);
     GTIFKeySet(gtif, GTRasterTypeGeoKey, TYPE_SHORT, 1, RasterPixelIsArea);
     GTIFKeySet(gtif, GeographicTypeGeoKey, TYPE_SHORT, 1, GCS_WGS_84);

     GTIFWriteKeys(gtif);

     return true;
}

void OutFile_tiff_color::writeLine() {
    for(int i=0; i<width; i++)
        fileData[i] = (uint8_t) round(productStuff[0]->calcOutputVal(productStuff[0]->lineData[i]));

    if (TIFFWriteScanline(tiff, (void*)fileData, currentLine) == 0) {
        fprintf(stderr, "Could not write TIFF image\n");
        exit(EXIT_FAILURE);
    }
    currentLine++;
}

bool OutFile_tiff_color::close() {
    GTIFFree(gtif);
    XTIFFClose(tiff);
    return true;
}

//-------------------------------------------------------------------------------------------
// OutFile_tiff_gray
//-------------------------------------------------------------------------------------------

OutFile_tiff_gray::OutFile_tiff_gray() : OutFile() {
    fileData = NULL;
    tiff = NULL;
    gtif = NULL;
}

OutFile_tiff_gray::~OutFile_tiff_gray() {
    if(fileData)
        free(fileData);
}

void OutFile_tiff_gray::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if(fileData)
        free(fileData);
    fileData = (float*) allocateMemory(width*sizeof(float), "OutFile_tiff_gray::fileData");
}

bool OutFile_tiff_gray::open() {
    currentLine = 0;
    tiff = XTIFFOpen(fileName.c_str(), "w");
     if (tiff == NULL) {
         fprintf(stderr, "-E- Could not open outgoing TIFF image\n");
         exit(EXIT_FAILURE);
     }
     gtif = GTIFNew(tiff);
     if (gtif == NULL) {
         fprintf(stderr, "-E- Could not create geoTIFF structure\n");
         exit(EXIT_FAILURE);
     }

     // calc geo TIFF  tags
     double tiepoints[6] = {0, 0, 0, 0, 0, 0};
     double pixscale[3] = {0, 0, 0};

     double north, south, east, west;
     if(metaData) {
         north = metaData->north;
         south = metaData->south;
         east = metaData->east;
         west = metaData->west;
     } else {
         north = 90.0;
         south = -90.0;
         east = 180.0;
         west = -180.0;
     }

     // pixel width
     pixscale[0] = (east - west) / width;

     // pixel height
     pixscale[1] = (north - south) / height;

     // set top left corner pixel lat, lon
     tiepoints[3] = west;
     tiepoints[4] = north;

     TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, width);
     TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, height);
     TIFFSetField(tiff, GTIFF_PIXELSCALE, 3, pixscale);
     TIFFSetField(tiff, GTIFF_TIEPOINTS, 6, tiepoints);

     // Grayscale
     TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
     TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
     TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 32);
     TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 1);
     TIFFSetField(tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);

     // write geo TIFF keys
     GTIFKeySet(gtif, GTModelTypeGeoKey, TYPE_SHORT, 1, ModelGeographic);
     GTIFKeySet(gtif, GTRasterTypeGeoKey, TYPE_SHORT, 1, RasterPixelIsArea);
     GTIFKeySet(gtif, GeographicTypeGeoKey, TYPE_SHORT, 1, GCS_WGS_84);

     GTIFWriteKeys(gtif);

     return true;
}

void OutFile_tiff_gray::writeLine() {
    for(int i=0; i<width; i++)
        fileData[i] = productStuff[0]->lineData[i];

    if (TIFFWriteScanline(tiff, (void*)fileData, currentLine) == 0) {
        fprintf(stderr, "Could not write TIFF image\n");
        exit(EXIT_FAILURE);
    }
    currentLine++;
}

bool OutFile_tiff_gray::close() {
    GTIFFree(gtif);
    XTIFFClose(tiff);
    return true;
}

//-------------------------------------------------------------------------------------------
// OutFile_tiff_rgb
//-------------------------------------------------------------------------------------------

OutFile_tiff_rgb::OutFile_tiff_rgb() : OutFile() {
    fileData = NULL;
    tiff = NULL;
    gtif = NULL;
}

OutFile_tiff_rgb::~OutFile_tiff_rgb() {
    if(fileData)
        free(fileData);
}

void OutFile_tiff_rgb::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if(fileData)
        free(fileData);
    fileData = (uint8_t*) allocateMemory(width*3, "OutFile_tiff_rgb::fileData");
}

bool OutFile_tiff_rgb::open() {
    currentLine = 0;
    tiff = XTIFFOpen(fileName.c_str(), "w");
    if (tiff == NULL) {
        fprintf(stderr, "-E- Could not open outgoing TIFF image\n");
        exit(EXIT_FAILURE);
    }
    gtif = GTIFNew(tiff);
    if (gtif == NULL) {
        fprintf(stderr, "-E- Could not create geoTIFF structure\n");
        exit(EXIT_FAILURE);
    }

    // calc geo TIFF  tags
    double tiepoints[6] = {0, 0, 0, 0, 0, 0};
    double pixscale[3] = {0, 0, 0};

    double north, south, east, west;
    if (metaData) {
        north = metaData->north;
        south = metaData->south;
        east = metaData->east;
        west = metaData->west;
    } else {
        north = 90.0;
        south = -90.0;
        east = 180.0;
        west = -180.0;
    }

    // pixel width
    pixscale[0] = (east - west) / width;

    // pixel height
    pixscale[1] = (north - south) / height;

    // set top left corner pixel lat, lon
    tiepoints[3] = west;
    tiepoints[4] = north;

    TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField(tiff, GTIFF_PIXELSCALE, 3, pixscale);
    TIFFSetField(tiff, GTIFF_TIEPOINTS, 6, tiepoints);

    // RGB
    TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 3);
    TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(tiff, TIFFTAG_COMPRESSION, COMPRESSION_DEFLATE);

    // write geo TIFF keys
    GTIFKeySet(gtif, GTModelTypeGeoKey, TYPE_SHORT, 1, ModelGeographic);
    GTIFKeySet(gtif, GTRasterTypeGeoKey, TYPE_SHORT, 1, RasterPixelIsArea);
    GTIFKeySet(gtif, GeographicTypeGeoKey, TYPE_SHORT, 1, GCS_WGS_84);

    GTIFWriteKeys(gtif);

    return true;
}

bool OutFile_tiff_rgb::close() {
    GTIFFree(gtif);
    XTIFFClose(tiff);
    return true;
}

void OutFile_tiff_rgb::setPixel(int32_t x, double val, int32_t prod) {
    fprintf(stderr, "-E- OutFile_tiff_rgb::setPixel - only RGB is implemented with this file type.\n");
    exit(EXIT_FAILURE);
}

void OutFile_tiff_rgb::setPixelRGB(int32_t x, float red, float green, float blue) {
    if(x<0 || x>=width) {
        fprintf(stderr, "-E- OutFile_tiff_rgb::setPixel - x=%d is not within range, width=%d.\n",
                x, width);
        exit(EXIT_FAILURE);
    }
        
    uint8_t *ptr = fileData + x*3;

    *ptr = round(productStuff[0]->calcOutputVal(red));
    ptr++;
    *ptr = round(productStuff[0]->calcOutputVal(green));
    ptr++;
    *ptr = round(productStuff[0]->calcOutputVal(blue));
    
    // do this to keep the file min/max reasonable
    if(red > fileMaxVal)
        fileMaxVal = red;
    if(red < fileMinVal)
        fileMinVal = red;
}

void OutFile_tiff_rgb::fillPixel(int32_t x) {
    if(x<0 || x>=width) {
        fprintf(stderr, "-E- OutFile_tiff_rgb::fillPixel - x=%d is not within range, width=%d.\n",
                x, width);
        exit(EXIT_FAILURE);
    }
    uint8_t *ptr = fileData + x*3;
    *ptr = 255;
    ptr++;
    *ptr = 255;
    ptr++;
    *ptr = 255;
}

void OutFile_tiff_rgb::missingPixel(int32_t x) {
    if(x<0 || x>=width) {
        fprintf(stderr, "-E- OutFile_tiff_rgb::missingPixel - x=%d is not within range, width=%d.\n",
                x, width);
        exit(EXIT_FAILURE);
    }
    uint8_t *ptr = fileData + x*3;
    *ptr = 255;
    ptr++;
    *ptr = 255;
    ptr++;
    *ptr = 255;
}

void OutFile_tiff_rgb::writeLine() {
   if (TIFFWriteScanline(tiff, (void*)fileData, currentLine) == 0) {
        fprintf(stderr, "Could not write TIFF image line\n");
        exit(EXIT_FAILURE);
    }
    currentLine++;
}

//-------------------------------------------------------------------------------------------
// OutFile_hdf4
//-------------------------------------------------------------------------------------------

OutFile_hdf4::OutFile_hdf4() : OutFile() {
    fileData = NULL;
    sdfid = -1;
    sdsid = -1;
    quality_sdsid = -1;
    hdfDataType = DFNT_FLOAT32;
}

OutFile_hdf4::~OutFile_hdf4() {
    if(fileData)
        free(fileData);
}

void OutFile_hdf4::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if(fileData)
        free(fileData);
    fileData = NULL;
}

bool OutFile_hdf4::open() {
    const char* tmpStr;
    float tmpFloat;
    int32 tmpInt32;

    currentLine = 0;

    sdfid = SDstart(fileName.c_str(), DFACC_CREATE);

    if (sdfid < 0) {
        printf("-E- Could not create HDF4 file %s\n", fileName.c_str());
        exit(EXIT_FAILURE);
    }

    get_time(metaData->ptime);

    string prodName;
    size_t pos = fileName.find_last_of('/');
    if(pos == string::npos)
        prodName = fileName;
    else
        prodName = fileName.substr(pos+1);
    DPTB(SDsetattr(sdfid, "Product Name", DFNT_CHAR, prodName.size() + 1, (VOIDP) prodName.c_str()));
    DPTB(SDsetattr(sdfid, "Sensor Name", DFNT_CHAR, strlen(metaData->sensor_name) + 1, (VOIDP) metaData->sensor_name));
    DPTB(SDsetattr(sdfid, "Sensor", DFNT_CHAR, strlen(metaData->sensor) + 1, (VOIDP) metaData->sensor));
    DPTB(SDsetattr(sdfid, "Title", DFNT_CHAR, strlen(metaData->title) + 1, (VOIDP) metaData->title));
    DPTB(SDsetattr(sdfid, "Data Center", DFNT_CHAR, strlen(metaData->data_center) + 1, (VOIDP) metaData->data_center));
    DPTB(SDsetattr(sdfid, "Station Name", DFNT_CHAR, strlen(metaData->station) + 1, (VOIDP) metaData->station));
    DPTB(SDsetattr(sdfid, "Station Latitude", DFNT_FLOAT32, 1, (VOIDP) &metaData->station_lat));
    DPTB(SDsetattr(sdfid, "Station Longitude", DFNT_FLOAT32, 1, (VOIDP) &metaData->station_lon));
    DPTB(SDsetattr(sdfid, "Mission", DFNT_CHAR, strlen(metaData->mission) + 1, (VOIDP) metaData->mission));
    DPTB(SDsetattr(sdfid, "Mission Characteristics", DFNT_CHAR, strlen(metaData->mission_char) + 1, (VOIDP) metaData->mission_char));
    DPTB(SDsetattr(sdfid, "Sensor Characteristics", DFNT_CHAR, strlen(metaData->sensor_char) + 1, (VOIDP) metaData->sensor_char));
    DPTB(SDsetattr(sdfid, "Product Type", DFNT_CHAR, strlen(metaData->prod_type) + 1, (VOIDP) metaData->prod_type));
    DPTB(SDsetattr(sdfid, "Processing Version", DFNT_CHAR, strlen(metaData->pversion) + 1, (VOIDP) metaData->pversion));
    DPTB(SDsetattr(sdfid, "Software Name", DFNT_CHAR, strlen(metaData->soft_name) + 1, (VOIDP) metaData->soft_name));
    DPTB(SDsetattr(sdfid, "Software Version", DFNT_CHAR, strlen(metaData->soft_ver) + 1, (VOIDP) metaData->soft_ver));
    DPTB(SDsetattr(sdfid, "Processing Time", DFNT_CHAR, strlen(metaData->ptime) + 1, (VOIDP) metaData->ptime));
    DPTB(SDsetattr(sdfid, "Input Files", DFNT_CHAR, strlen(metaData->infiles) + 1, (VOIDP) metaData->infiles));
    DPTB(SDsetattr(sdfid, "Processing Control", DFNT_CHAR, strlen(metaData->proc_con) + 1, (VOIDP) metaData->proc_con));
    DPTB(SDsetattr(sdfid, "Input Parameters", DFNT_CHAR, strlen(metaData->input_parms) + 1, (VOIDP) metaData->input_parms));
    DPTB(SDsetattr(sdfid, "L2 Flag Names", DFNT_CHAR, strlen(metaData->flag_names) + 1, (VOIDP) metaData->flag_names));

    short syear, sday, eyear, eday;
    double ssec, esec;
    int32 smsec, emsec;
    unix2yds(metaData->startTime, &syear, &sday, &ssec);
    smsec = (int32) (ssec * 1000.0);
    unix2yds(metaData->endTime, &eyear, &eday, &esec);
    emsec = (int32) (esec * 1000.0);
    DPTB(SDsetattr(sdfid, "Period Start Year", DFNT_INT16, 1, (VOIDP) &syear));
    DPTB(SDsetattr(sdfid, "Period Start Day", DFNT_INT16, 1, (VOIDP) &sday));
    DPTB(SDsetattr(sdfid, "Period End Year", DFNT_INT16, 1, (VOIDP) &eyear));
    DPTB(SDsetattr(sdfid, "Period End Day", DFNT_INT16, 1, (VOIDP) &eday));
    tmpStr = ydhmsf(metaData->startTime,'G');
    DPTB(SDsetattr(sdfid, "Start Time", DFNT_CHAR, strlen(tmpStr) + 1, (VOIDP) tmpStr));
    tmpStr = ydhmsf(metaData->endTime,'G');
    DPTB(SDsetattr(sdfid, "End Time", DFNT_CHAR, strlen(tmpStr) + 1, (VOIDP) tmpStr));
    DPTB(SDsetattr(sdfid, "Start Year", DFNT_INT16, 1, (VOIDP) &syear));
    DPTB(SDsetattr(sdfid, "Start Day", DFNT_INT16, 1, (VOIDP) &sday));
    DPTB(SDsetattr(sdfid, "Start Millisec", DFNT_INT32, 1, (VOIDP) &smsec));
    DPTB(SDsetattr(sdfid, "End Year", DFNT_INT16, 1, (VOIDP) &eyear));
    DPTB(SDsetattr(sdfid, "End Day", DFNT_INT16, 1, (VOIDP) &eday));
    DPTB(SDsetattr(sdfid, "End Millisec", DFNT_INT32, 1, (VOIDP) &emsec));

    DPTB(SDsetattr(sdfid, "Start Orbit", DFNT_INT32, 1, (VOIDP) &metaData->start_orb));
    DPTB(SDsetattr(sdfid, "End Orbit", DFNT_INT32, 1, (VOIDP) &metaData->end_orb));
    DPTB(SDsetattr(sdfid, "Orbit", DFNT_INT32, 1, (VOIDP) &metaData->orbit));
    DPTB(SDsetattr(sdfid, "Map Projection", DFNT_CHAR, mapProjection.size()+1, mapProjection.c_str()));
    DPTB(SDsetattr(sdfid, "Latitude Units", DFNT_CHAR, 14, (VOIDP) "degrees North"));
    DPTB(SDsetattr(sdfid, "Longitude Units", DFNT_CHAR, 13, (VOIDP) "degrees East"));
    DPTB(SDsetattr(sdfid, "Northernmost Latitude", DFNT_FLOAT32, 1, (VOIDP) &metaData->north));
    DPTB(SDsetattr(sdfid, "Southernmost Latitude", DFNT_FLOAT32, 1, (VOIDP) &metaData->south));
    DPTB(SDsetattr(sdfid, "Westernmost Longitude", DFNT_FLOAT32, 1, (VOIDP) &metaData->west));
    DPTB(SDsetattr(sdfid, "Easternmost Longitude", DFNT_FLOAT32, 1, (VOIDP) &metaData->east));

    float latStep = (metaData->north - metaData->south) / (float)getHeight();
    DPTB(SDsetattr(sdfid, "Latitude Step", DFNT_FLOAT32, 1, (VOIDP) &latStep));
    float lonStep;
    if(metaData->east < metaData->west)
        lonStep = (360 + metaData->east - metaData->west) / (float)getWidth();
    else
        lonStep = (metaData->east - metaData->west) / (float)getWidth();
    DPTB(SDsetattr(sdfid, "Longitude Step", DFNT_FLOAT32, 1, (VOIDP) &lonStep));

    tmpFloat = metaData->south + latStep/2.0;
    DPTB(SDsetattr(sdfid, "SW Point Latitude", DFNT_FLOAT32, 1, (VOIDP) &tmpFloat));
    tmpFloat = metaData->west + lonStep/2.0;
    DPTB(SDsetattr(sdfid, "SW Point Longitude", DFNT_FLOAT32, 1, (VOIDP) &tmpFloat));

    DPTB(SDsetattr(sdfid, "Data Bins", DFNT_INT32, 1, (VOIDP) &metaData->data_bins));

    tmpInt32 = getHeight();
    DPTB(SDsetattr(sdfid, "Number of Lines", DFNT_INT32, 1, (VOIDP) &tmpInt32));
    tmpInt32 = getWidth();
    DPTB(SDsetattr(sdfid, "Number of Columns", DFNT_INT32, 1, (VOIDP) &tmpInt32));
    DPTB(SDsetattr(sdfid, "Parameter", DFNT_CHAR, strlen(productStuff[0]->productInfo->description) + 1, (VOIDP) productStuff[0]->productInfo->description));
    DPTB(SDsetattr(sdfid, "Measure", DFNT_CHAR, 5, (VOIDP) "Mean"));
    DPTB(SDsetattr(sdfid, "Units", DFNT_CHAR, strlen(productStuff[0]->productInfo->units) + 1, (VOIDP) productStuff[0]->productInfo->units));

    // we only use linear scaling for data storage
    tmpStr = "linear";
    DPTB(SDsetattr(sdfid, "Scaling", DFNT_CHAR, strlen(tmpStr) + 1, (VOIDP) tmpStr));
    tmpStr = "(Slope*l3m_data) + Intercept = Parameter value";
    DPTB(SDsetattr(sdfid, "Scaling Equation", DFNT_CHAR, strlen(tmpStr) + 1, (VOIDP) tmpStr));

    float tmpScale;
    float tmpOffset;
    const char* imageScalingApplied;
    
    switch(productStuff[0]->dataStorage) {
    case FloatDS:
    case DoubleDS:
        tmpScale = 1.0;
        tmpOffset = 0.0;
        imageScalingApplied = "No";
        break;
    default:
        tmpScale = productStuff[0]->scale;
        tmpOffset = productStuff[0]->offset;
        if(tmpScale==1.0 && tmpOffset==0.0)
            imageScalingApplied = "No";
        else
            imageScalingApplied = "Yes";
    }

    DPTB(SDsetattr(sdfid, "Slope", DFNT_FLOAT32, 1, (VOIDP) &tmpScale));
    DPTB(SDsetattr(sdfid, "Intercept", DFNT_FLOAT32, 1, (VOIDP) &tmpOffset));

    tmpFloat = productStuff[0]->productInfo->validMin;
    DPTB(SDsetattr(sdfid, "Data Minimum", DFNT_FLOAT32, 1, (VOIDP) &tmpFloat));
    tmpFloat = productStuff[0]->productInfo->validMax;
    DPTB(SDsetattr(sdfid, "Data Maximum", DFNT_FLOAT32, 1, (VOIDP) &tmpFloat));
    tmpFloat = productStuff[0]->productInfo->displayMin;
    DPTB(SDsetattr(sdfid, "Suggested Image Scaling Minimum", DFNT_FLOAT32, 1, (VOIDP) &tmpFloat));
    tmpFloat = productStuff[0]->productInfo->displayMax;
    DPTB(SDsetattr(sdfid, "Suggested Image Scaling Maximum", DFNT_FLOAT32, 1, (VOIDP) &tmpFloat));

    tmpStr = strdup(productStuff[0]->productInfo->displayScale);
    upcase((char*)tmpStr);
    DPTB(SDsetattr(sdfid, "Suggested Image Scaling Type", DFNT_CHAR, strlen(tmpStr) + 1, (VOIDP) tmpStr));
    free((void*)tmpStr);

    DPTB(SDsetattr(sdfid, "Suggested Image Scaling Applied", DFNT_CHAR, strlen(imageScalingApplied) + 1, (VOIDP) imageScalingApplied));
    DPTB(SDsetattr(sdfid, "_lastModified", DFNT_CHAR, strlen(metaData->ptime) + 1, (VOIDP) metaData->ptime));

    // delete file data
    if(fileData)
        free(fileData);

    int32 dims[2];
    dims[0] = height;
    dims[1] = width;

    if(!strcmp(productStuff[0]->productInfo->dataType, "byte")) {
    	hdfDataType = DFNT_INT8;
    	sdsid = SDcreate(sdfid, "l3m_data", hdfDataType, 2, dims);
    	int8 tmp = productStuff[0]->missingValue;
    	DPTB(SDsetattr(sdsid, "Fill", hdfDataType, 1, (VOIDP) &tmp));
        fileData = allocateMemory(width*sizeof(int8_t), "OutFile_hdf4::open fileData");
    } else if(!strcmp(productStuff[0]->productInfo->dataType, "ubyte")) {
        hdfDataType = DFNT_UINT8;
        sdsid = SDcreate(sdfid, "l3m_data", hdfDataType, 2, dims);
        uint8 tmp = productStuff[0]->missingValue;
        DPTB(SDsetattr(sdsid, "Fill", hdfDataType, 1, (VOIDP) &tmp));
        fileData = allocateMemory(width*sizeof(uint8_t), "OutFile_hdf4::open fileData");
    } else if(!strcmp(productStuff[0]->productInfo->dataType, "short")) {
    	hdfDataType = DFNT_INT16;
    	sdsid = SDcreate(sdfid, "l3m_data", hdfDataType, 2, dims);
    	int16 tmp = productStuff[0]->missingValue;
    	DPTB(SDsetattr(sdsid, "Fill", hdfDataType, 1, (VOIDP) &tmp));
        fileData = allocateMemory(width*sizeof(int16_t), "OutFile_hdf4::open fileData");
    } else if(!strcmp(productStuff[0]->productInfo->dataType, "ushort")) {
    	hdfDataType = DFNT_UINT16;
    	sdsid = SDcreate(sdfid, "l3m_data", hdfDataType, 2, dims);
    	uint16 tmp = productStuff[0]->missingValue;
    	DPTB(SDsetattr(sdsid, "Fill", hdfDataType, 1, (VOIDP) &tmp));
        fileData = allocateMemory(width*sizeof(uint16_t), "OutFile_hdf4::open fileData");
   } else if(!strcmp(productStuff[0]->productInfo->dataType, "int")) {
    	hdfDataType = DFNT_INT32;
    	sdsid = SDcreate(sdfid, "l3m_data", hdfDataType, 2, dims);
    	int32 tmp = productStuff[0]->missingValue;
    	DPTB(SDsetattr(sdsid, "Fill", hdfDataType, 1, (VOIDP) &tmp));
        fileData = allocateMemory(width*sizeof(int32_t), "OutFile_hdf4::open fileData");
    } else if(!strcmp(productStuff[0]->productInfo->dataType, "uint")) {
    	hdfDataType = DFNT_UINT32;
    	sdsid = SDcreate(sdfid, "l3m_data", hdfDataType, 2, dims);
    	uint32 tmp = productStuff[0]->missingValue;
    	DPTB(SDsetattr(sdsid, "Fill", hdfDataType, 1, (VOIDP) &tmp));
        fileData = allocateMemory(width*sizeof(uint32_t), "OutFile_hdf4::open fileData");
   } else if(!strcmp(productStuff[0]->productInfo->dataType, "float")) {
        hdfDataType = DFNT_FLOAT32;
        sdsid = SDcreate(sdfid, "l3m_data", hdfDataType, 2, dims);
        float32 tmp = productStuff[0]->missingValue;
        DPTB(SDsetattr(sdsid, "Fill", hdfDataType, 1, (VOIDP) &tmp));
        fileData = allocateMemory(width*sizeof(float), "OutFile_hdf4::open fileData");
    } else if(!strcmp(productStuff[0]->productInfo->dataType, "double")) {
        hdfDataType = DFNT_FLOAT64;
        sdsid = SDcreate(sdfid, "l3m_data", hdfDataType, 2, dims);
        float64 tmp = productStuff[0]->missingValue;
        DPTB(SDsetattr(sdsid, "Fill", hdfDataType, 1, (VOIDP) &tmp));
        fileData = allocateMemory(width*sizeof(double), "OutFile_hdf4::open fileData");
   } else {
        printf("-E- Data type %s, not supported\n", productStuff[0]->productInfo->dataType);
        exit(EXIT_FAILURE);
    }

    tmpStr = "linear";
    DPTB(SDsetattr(sdsid, "Scaling", DFNT_CHAR, strlen(tmpStr) + 1, (VOIDP) tmpStr));
    tmpStr = "(Slope*l3m_data) + Intercept = Parameter value";
    DPTB(SDsetattr(sdsid, "Scaling Equation", DFNT_CHAR, strlen(tmpStr) + 1, (VOIDP) tmpStr));
    DPTB(SDsetattr(sdsid, "Slope", DFNT_FLOAT32, 1, (VOIDP) &tmpScale));
    DPTB(SDsetattr(sdsid, "Intercept", DFNT_FLOAT32, 1, (VOIDP) &tmpOffset));

    // create the SST quality data set
    if(qualityData) {
        quality_sdsid = SDcreate(sdfid, "l3m_qual", DFNT_UINT8, 2, dims);
        int32 validRange[2];
        validRange[0] = 0;
        validRange[1] = 2;
        DPTB(SDsetattr(quality_sdsid, "valid_range", DFNT_INT32, 2, (VOIDP) validRange));
    }
    
    return true;
}

void OutFile_hdf4::writeLine() {

    switch(productStuff[0]->dataStorage) {
    case ByteDS:
        for(int i=0; i<width; i++)
            ((int8_t*)fileData)[i] = round(productStuff[0]->calcOutputVal(productStuff[0]->lineData[i]));
        break;
    case UByteDS:
        for(int i=0; i<width; i++)
            ((uint8_t*)fileData)[i] = round(productStuff[0]->calcOutputVal(productStuff[0]->lineData[i]));
        break;
    case ShortDS:
        for(int i=0; i<width; i++)
            ((int16_t*)fileData)[i] = round(productStuff[0]->calcOutputVal(productStuff[0]->lineData[i]));
        break;
    case UShortDS:
        for(int i=0; i<width; i++)
            ((uint16_t*)fileData)[i] = round(productStuff[0]->calcOutputVal(productStuff[0]->lineData[i]));
        break;
    case IntDS:
        for(int i=0; i<width; i++)
            ((int32_t*)fileData)[i] = round(productStuff[0]->calcOutputVal(productStuff[0]->lineData[i]));
        break;
    case UIntDS:
        for(int i=0; i<width; i++)
            ((uint32_t*)fileData)[i] = round(productStuff[0]->calcOutputVal(productStuff[0]->lineData[i]));
        break;
    case FloatDS:
        for(int i=0; i<width; i++)
            ((float*)fileData)[i] = productStuff[0]->lineData[i];
        break;
    case DoubleDS:
        for(int i=0; i<width; i++)
            ((double*)fileData)[i] = productStuff[0]->lineData[i];
        break;
    default:
        printf("-E- OutFile_hdf4::writeLine - unrecognized data type = %d\n", productStuff[0]->dataStorage);
        exit(EXIT_FAILURE);
    }

    int32 start[2];
    int32 count[2];

    start[0] = currentLine;
    start[1] = 0;
    count[0] = 1;
    count[1] = width;

    if((SDwritedata(sdsid, start, NULL, count, (VOIDP)fileData)) < 0) {
        printf("\n-E- OutFile_hdf4::writeLine(): SDwritedata unsuccessful\n");
        exit(EXIT_FAILURE);
    }

    if(qualityData) {
        if((SDwritedata(quality_sdsid, start, NULL, count, (VOIDP)qualityData)) < 0) {
            printf("\n-E- OutFile_hdf4::writeLine(): SDwritedata unsuccessful\n");
            exit(EXIT_FAILURE);
        }
    }
    
    currentLine++;
}

bool OutFile_hdf4::close() {

    if(fileData)
        free(fileData);

    if(metaData) {
        DPTB(SDsetattr(sdfid, "Data Bins", DFNT_INT32, 1, (VOIDP) &metaData->data_bins));
    }
    float tmpFloat = getFileMinVal();
    DPTB(SDsetattr(sdfid, "Data Minimum", DFNT_FLOAT32, 1, (VOIDP) &tmpFloat));
    tmpFloat = getFileMaxVal();
    DPTB(SDsetattr(sdfid, "Data Maximum", DFNT_FLOAT32, 1, (VOIDP) &tmpFloat));

    if(qualityData) {
        SDendaccess(quality_sdsid);
	quality_sdsid = -1;
    }
    SDendaccess(sdsid);
    sdsid = -1;
    SDend(sdfid);
    sdfid = -1;

    /*-----------------------------------------------------------------------*/
    /*  write map_palette */

    // make sure a palette has been loaded
    if(red && green && blue) {
        uint8_t data[768];
        int j=0;
        for (int i = 0; i < 256; i++) {
            data[j++] = red[i];
            data[j++] = green[i];
            data[j++] = blue[i];
        }


        if ((DFPaddpal(fileName.c_str(), (VOIDP)data)) < 0) {
            printf("-E- OutFile_hdf4::close - Error writing map_palette.\n");
            return false;
        }

        uint16 pal_ref;
        if ((pal_ref = DFPlastref()) > 0) {
            if ((DFANputlabel(fileName.c_str(), DFTAG_IP8, pal_ref, (char*)"palette")) < 0) {
                printf("-E- OutFile_hdf4::close - Error writing palette label\n");
                return false;
            }
        }
    }
    return true;
}

int32_t OutFile_hdf4::addProduct(productInfo_t* productInfo) {
    return addProductNonDisplay(productInfo);
}

//-------------------------------------------------------------------------------------------
// OutFile_netcdf4
//-------------------------------------------------------------------------------------------

OutFile_netcdf4::OutFile_netcdf4() : OutFile() {
    fileData = NULL;
    ncid = -1;
    qualityid = -1;
}

OutFile_netcdf4::~OutFile_netcdf4() {
    if(fileData)
        free(fileData);
}

void OutFile_netcdf4::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if(fileData)
        free(fileData);
}

bool OutFile_netcdf4::open() {
    const char* tmpStr;
    float tmpFloat;
    int status;
    char buf[2048];
    int varid1;

    currentLine = 0;

    status = nc_create(fileName.c_str(), NC_NETCDF4, &ncid);
    if (status) {
        printf("-E- Could not create netCDF4 file %s\n", fileName.c_str());
        exit(EXIT_FAILURE);
    }

    string prodName;
    size_t pos = fileName.find_last_of('/');
    if(pos == string::npos)
        prodName = fileName;
    else
        prodName = fileName.substr(pos+1);
    nc_put_att(ncid, NC_GLOBAL, "product_name", NC_CHAR, prodName.size()+1, prodName.c_str());
    nc_put_att(ncid, NC_GLOBAL, "instrument", NC_CHAR, strlen(metaData->sensor)+1, metaData->sensor);
    nc_put_att(ncid, NC_GLOBAL, "title", NC_CHAR, strlen(metaData->title)+1, metaData->title);

    strcpy(buf, "Ocean Biology Processing Group (NASA/GSFC/OBPG)");
    status = nc_put_att(ncid, NC_GLOBAL, "project", NC_CHAR, strlen(buf)+1, buf);

    if (strcmp(metaData->mission, "") != 0)
      nc_put_att(ncid, NC_GLOBAL, "platform",   NC_CHAR, strlen(metaData->mission)+1, metaData->mission);

    nc_put_att(ncid, NC_GLOBAL, "temporal_range",  NC_CHAR, strlen(metaData->prod_type)+1, metaData->prod_type);
    nc_put_att(ncid, NC_GLOBAL, "processing_version",  NC_CHAR, strlen(metaData->pversion)+1, metaData->pversion);

    strcpy(metaData->ptime, unix2isodate(time(NULL), 'G'));
    nc_put_att(ncid, NC_GLOBAL, "date_created", NC_CHAR, strlen(metaData->ptime)+1, metaData->ptime);
    nc_put_att(ncid, NC_GLOBAL, "history", NC_CHAR, strlen(metaData->proc_con)+1, metaData->proc_con);
    nc_put_att(ncid, NC_GLOBAL, "l2_flag_names",NC_CHAR, strlen(metaData->flag_names)+1, metaData->flag_names);

    strcpy(buf, unix2isodate(metaData->startTime,'G'));
    nc_put_att(ncid, NC_GLOBAL, "time_coverage_start", NC_CHAR, strlen(buf)+1, buf);

    strcpy(buf, unix2isodate(metaData->endTime,'G'));
    nc_put_att(ncid, NC_GLOBAL, "time_coverage_end", NC_CHAR, strlen(buf)+1, buf);
    nc_put_att(ncid, NC_GLOBAL, "start_orbit_number", NC_INT, 1, &metaData->start_orb);
    nc_put_att(ncid, NC_GLOBAL, "end_orbit_number", NC_INT, 1, &metaData->end_orb);
    nc_put_att(ncid, NC_GLOBAL, "map_projection", NC_CHAR, mapProjection.size()+1, mapProjection.c_str());

    nc_put_att(ncid, NC_GLOBAL, "latitude_units", NC_CHAR, strlen(metaData->lat_units)+1, metaData->lat_units);
    nc_put_att(ncid, NC_GLOBAL, "longitude_units", NC_CHAR, strlen(metaData->lon_units)+1, metaData->lon_units);

    nc_put_att(ncid, NC_GLOBAL, "northernmost_latitude", NC_FLOAT, 1, &metaData->north);
    nc_put_att(ncid, NC_GLOBAL, "southernmost_latitude", NC_FLOAT, 1, &metaData->south);
    nc_put_att(ncid, NC_GLOBAL, "westernmost_longitude", NC_FLOAT, 1, &metaData->west);
    nc_put_att(ncid, NC_GLOBAL, "easternmost_longitude", NC_FLOAT, 1, &metaData->east);
    nc_put_att(ncid, NC_GLOBAL, "geospatial_lat_max", NC_FLOAT, 1, &metaData->north);
    nc_put_att(ncid, NC_GLOBAL, "geospatial_lat_min", NC_FLOAT, 1, &metaData->south);
    nc_put_att(ncid, NC_GLOBAL, "geospatial_lon_max", NC_FLOAT, 1, &metaData->east);
    nc_put_att(ncid, NC_GLOBAL, "geospatial_lon_min", NC_FLOAT, 1, &metaData->west);

    tmpStr = "latitude_longitude";
    nc_put_att(ncid, NC_GLOBAL, "grid_mapping_name", NC_CHAR, strlen(tmpStr)+1, tmpStr);

    float latStep = (metaData->north - metaData->south) / (float)getHeight();
    float lonStep;
    if(metaData->east < metaData->west)
        lonStep = (360 + metaData->east - metaData->west) / (float)getWidth();
    else
        lonStep = (metaData->east - metaData->west) / (float)getWidth();
    nc_put_att(ncid, NC_GLOBAL, "latitude_step", NC_FLOAT, 1, &latStep);
    nc_put_att(ncid, NC_GLOBAL, "longitude_step", NC_FLOAT, 1, &lonStep);

    tmpFloat = (metaData->south + latStep/2.0);
    nc_put_att(ncid, NC_GLOBAL, "sw_point_latitude", NC_FLOAT, 1, &tmpFloat);
    tmpFloat = (metaData->west + lonStep/2.0);
    nc_put_att(ncid, NC_GLOBAL, "sw_point_longitude", NC_FLOAT, 1, &tmpFloat);

    const char* geoUnits;
    if(resolution > 1000.0) {
        tmpFloat = resolution / 1000.0;
        geoUnits = "km";
    } else {
        tmpFloat = resolution;
        geoUnits = "m";
    }
    nc_put_att(ncid, NC_GLOBAL, "geospatial_lon_resolution", NC_FLOAT, 1, &tmpFloat);
    nc_put_att(ncid, NC_GLOBAL, "geospatial_lat_resolution", NC_FLOAT, 1, &tmpFloat);
    nc_put_att(ncid, NC_GLOBAL, "geospatial_lat_units", NC_CHAR, strlen(metaData->lat_units)+1, metaData->lat_units);
    nc_put_att(ncid, NC_GLOBAL, "geospatial_lon_units", NC_CHAR, strlen(metaData->lon_units)+1, metaData->lon_units);
    sprintf(buf, "%3.2f %s", tmpFloat, geoUnits);
    nc_put_att(ncid, NC_GLOBAL, "spatialResolution", NC_CHAR, strlen(buf)+1, buf);

    nc_put_att(ncid, NC_GLOBAL, "number_of_lines", NC_INT, 1, &height);
    nc_put_att(ncid, NC_GLOBAL, "number_of_columns", NC_INT, 1, &width);
    nc_put_att(ncid, NC_GLOBAL, "measure", NC_CHAR, 5, "Mean");
    tmpFloat = productStuff[0]->productInfo->displayMin;
    nc_put_att(ncid, NC_GLOBAL, "suggested_image_scaling_minimum", NC_FLOAT, 1, &tmpFloat);
    tmpFloat = productStuff[0]->productInfo->displayMax;
    nc_put_att(ncid, NC_GLOBAL, "suggested_image_scaling_maximum", NC_FLOAT, 1, &tmpFloat);

    if(!strcmp(productStuff[0]->productInfo->displayScale, "log"))
        tmpStr = "LOG";
    else if(!strcmp(productStuff[0]->productInfo->displayScale, "arctan"))
        tmpStr = "ATAN";
    else
        tmpStr = "LINEAR";

    nc_put_att(ncid, NC_GLOBAL, "suggested_image_scaling_type", NC_CHAR, strlen(tmpStr)+1, tmpStr);
    nc_put_att(ncid, NC_GLOBAL, "suggested_image_scaling_applied", NC_CHAR, 3, "No");
    nc_put_att(ncid, NC_GLOBAL, "_lastModified", NC_CHAR, strlen(metaData->ptime)+1, (VOIDP)metaData->ptime);
    tmpStr = "CF-1.6";
    nc_put_att(ncid, NC_GLOBAL, "Conventions", NC_CHAR, strlen(tmpStr)+1, tmpStr);
    tmpStr = "NASA Goddard Space Flight Center, Ocean Ecology Laboratory, Ocean Biology Processing Group";
    nc_put_att(ncid, NC_GLOBAL, "institution", NC_CHAR, strlen(tmpStr)+1, tmpStr);
    tmpStr = "NetCDF Climate and Forecast (CF) Metadata Convention";
    nc_put_att(ncid, NC_GLOBAL, "standard_name_vocabulary", NC_CHAR, strlen(tmpStr)+1, tmpStr);
    tmpStr = "Unidata Dataset Discovery v1.0";
    nc_put_att(ncid, NC_GLOBAL, "Metadata_Conventions", NC_CHAR, strlen(tmpStr)+1, tmpStr);
    tmpStr = "gov.nasa.gsfc.sci.oceandata";
    nc_put_att(ncid, NC_GLOBAL, "naming_authority", NC_CHAR, strlen(tmpStr)+1, tmpStr);

    // create id
    strcpy(buf,metaData->pversion);
    if (strcmp(metaData->pversion,"Unspecified") != 0){
        strcpy(buf, metaData->product_name);
        strcat(buf, "/L3/");
    } else {
        strcpy(buf,"L3/");
    }
    strcat(buf,metaData->product_name);
    nc_put_att(ncid, NC_GLOBAL, "id", NC_CHAR, strlen(buf)+1, buf);

    tmpStr = "http://science.nasa.gov/earth-science/earth-science-data/data-information-policy/";
    nc_put_att(ncid, NC_GLOBAL, "license", NC_CHAR, strlen(tmpStr)+1, tmpStr);
    tmpStr = "NASA/GSFC/OBPG";
    nc_put_att(ncid, NC_GLOBAL, "creator_name", NC_CHAR, strlen(tmpStr)+1, tmpStr);
    nc_put_att(ncid, NC_GLOBAL, "publisher_name", NC_CHAR, strlen(tmpStr)+1, tmpStr);

    strcpy(buf, "data@oceancolor.gsfc.nasa.gov");
    nc_put_att(ncid, NC_GLOBAL, "creator_email", NC_CHAR, strlen(buf)+1, buf);
    nc_put_att(ncid, NC_GLOBAL, "publisher_email", NC_CHAR, strlen(buf)+1, buf);

    strcpy(buf, "http://oceandata.sci.gsfc.nasa.gov");
    nc_put_att(ncid, NC_GLOBAL, "creator_url", NC_CHAR, strlen(buf)+1, buf);
    nc_put_att(ncid, NC_GLOBAL, "publisher_url", NC_CHAR, strlen(buf)+1, buf);

    strcpy(buf, "L3 Mapped");
    status = nc_put_att(ncid, NC_GLOBAL, "processing_level", NC_CHAR, strlen(buf)+1, buf);

    strcpy(buf, "grid");
    status = nc_put_att(ncid, NC_GLOBAL, "cdm_data_type", NC_CHAR, strlen(buf)+1, buf);

    strcpy(buf, "http://dx.doi.org");
    status = nc_put_att(ncid, NC_GLOBAL, "identifier_product_doi_authority", NC_CHAR, strlen(buf)+1, buf);

    if (strcmp(metaData->sensor_name, "CZCS") == 0)
        strcpy(buf, "10.5067/NIMBUS-7/CZCS_OC.2014.0");
    if (strcmp(metaData->sensor_name, "OCTS") == 0)
        strcpy(buf, "10.5067/ADEOS/OCTS_OC.2014.0");
    if (strcmp(metaData->sensor_name, "SEAWIFS") == 0)
        strcpy(buf, "10.5067/ORBVIEW-2/SEAWIFS_OC.2014.0");
    if (strstr(metaData->sensor_name, "MODIS") != NULL){
        if (strstr(metaData->mission, "Aqua") != NULL)
            strcpy(buf, "10.5067/AQUA/MODIS_OC.2014.0");
        if (strstr(metaData->mission, "Terra") != NULL)
            strcpy(buf, "10.5067/TERRA/MODIS_OC.2014.0");
    }
    nc_put_att(ncid, NC_GLOBAL, "identifier_product_doi", NC_CHAR, strlen(buf)+1, buf);

    if ( strstr(productStuff[0]->productInfo->productName, "sst") != NULL)
        tmpStr = "Oceans > Ocean Temperature > Sea Surface Temperature";
    else
        tmpStr = "Oceans > Ocean Chemistry > Chlorophyll; Oceans > Ocean Optics > Ocean Color";
    nc_put_att(ncid, NC_GLOBAL, "keywords", NC_CHAR, strlen(tmpStr)+1, tmpStr);

    tmpStr = "NASA Global Change Master Directory (GCMD) Science Keywords";
    nc_put_att(ncid, NC_GLOBAL, "keywords_vocabulary", NC_CHAR, strlen(tmpStr)+1, tmpStr);

    int grpid;
    nc_def_grp(ncid, "processing_control", &grpid);

    status = nc_put_att(grpid, NC_GLOBAL, "software_name", NC_CHAR, strlen(metaData->soft_name)+1, metaData->soft_name);
    check_err(status,__LINE__,__FILE__);
    status = nc_put_att(grpid, NC_GLOBAL, "software_version", NC_CHAR, strlen(metaData->soft_ver)+1, metaData->soft_ver);
    check_err(status,__LINE__,__FILE__);

    if ((tmpStr = strrchr(metaData->infiles, '/')) != NULL)
        tmpStr++;
    else
        tmpStr = metaData->infiles;
    status = nc_put_att(grpid, NC_GLOBAL, "source", NC_CHAR, strlen(tmpStr)+1, tmpStr);
    check_err(status,__LINE__,__FILE__);
    status = nc_put_att(grpid, NC_GLOBAL, "l2_flag_names", NC_CHAR, strlen(metaData->flag_names)+1, metaData->flag_names);
    check_err(status,__LINE__,__FILE__);

    nc_def_grp(grpid, "input_parameters", &grpid);
    char *end_str;
    char *token = strtok_r(metaData->input_parms, "|", &end_str);
    while (token != NULL) {
        char *end_token;
        strcpy(buf, token);
        char *name = strtok_r(token, "=", &end_token);
        for (uint32_t i=0; i<strlen(name); i++) {
            if (name[i] == ' ') {
                name[i] = 0;
                break;
            }
        }
        strcpy(buf, strtok_r(NULL, "|", &end_token));
        nc_put_att(grpid, NC_GLOBAL, name, NC_CHAR, strlen(buf)+1, buf);
        token = strtok_r(NULL, "|", &end_str);
    }

    if(fileData)
        free(fileData);

    int dimids[3];
    int dim_sizes[3];

    dim_sizes[0] = height;
    dim_sizes[1] = width;

    // create the dimensions
    status = nc_def_dim(ncid, "lat", dim_sizes[0], &dimids[0]);
    check_err(status,__LINE__,__FILE__);

    status = nc_def_dim(ncid, "lon", dim_sizes[1], &dimids[1]);
    check_err(status,__LINE__,__FILE__);

    idDS ds_id;
    ds_id.deflate = deflate;
    ds_id.fftype = DS_NCDF;
    ds_id.fid = ncid;

    size_t dataSize = 1;
    
    for(size_t prodNum=0; prodNum<productStuff.size(); prodNum++) {

        int varDataType;
        switch(productStuff[prodNum]->dataStorage) {
        case ByteDS:
            varDataType = NC_BYTE;
            if(sizeof(signed char) > dataSize)
                dataSize = sizeof(signed char);
            break;
        case UByteDS:
            varDataType = NC_UBYTE;
            if(sizeof(unsigned char) > dataSize)
                dataSize = sizeof(unsigned char);
            break;
        case ShortDS:
            varDataType = NC_SHORT;
            if(sizeof(short) > dataSize)
                dataSize = sizeof(short);
            break;
        case UShortDS:
            varDataType = NC_USHORT;
            if(sizeof(unsigned short) > dataSize)
                dataSize = sizeof(unsigned short);
            break;
        case IntDS:
            varDataType = NC_INT;
            if(sizeof(int) > dataSize)
                dataSize = sizeof(int);
            break;
        case UIntDS:
            varDataType = NC_UINT;
            if(sizeof(unsigned int) > dataSize)
                dataSize = sizeof(unsigned int);
            break;
        case FloatDS:
            varDataType = NC_FLOAT;
            if(sizeof(float) > dataSize)
                dataSize = sizeof(float);
            break;
        case DoubleDS:
            varDataType = NC_DOUBLE;
            if(sizeof(double) > dataSize)
                dataSize = sizeof(double);
            break;
        default:
            printf("-E- OutFile_netcdf4::open - illegal data storage type = %d\n", productStuff[prodNum]->dataStorage);
            exit(EXIT_FAILURE);
        }

        char* productNameStr = getProductNameFull(productStuff[prodNum]->productInfo);

        status = CreateNCDF(ds_id,
                productNameStr,    /* short name */
                productStuff[prodNum]->productInfo->description,    /* long name */
                productStuff[prodNum]->productInfo->standardName, /* NetCDF standard name (not set if passed NULL or "") */
                productStuff[prodNum]->productInfo->reference,
                productStuff[prodNum]->productInfo->comment,
                productStuff[prodNum]->productInfo->units,    /* units (not set if passed NULL or "") */
                productStuff[prodNum]->productInfo->validMin,     /* low end of valid range */
                productStuff[prodNum]->productInfo->validMax,    /* high end of range (no range set if low >= high) */
                productStuff[prodNum]->productInfo->scaleFactor,    /* scale factor (not set if 1.0)  */
                productStuff[prodNum]->productInfo->addOffset,   /* scaling offset (not set if 0.0)  */
                productStuff[prodNum]->productInfo->fillValue,       /* fill value */
                varDataType,       /* NCDF number type */
                2,     /* number of dimensions (must be <= 3) */
                dimids /* dimension ids */
        );
        check_err(status,__LINE__,__FILE__);

        int tmpVarId;
        status = nc_inq_varid(ncid, productNameStr, &tmpVarId);
        check_err(status,__LINE__,__FILE__);

        varid.push_back(tmpVarId);
        
        // add display info
        status = nc_put_att_text(ncid, tmpVarId, "display_scale",
                      strlen(productStuff[prodNum]->productInfo->displayScale) + 1,
                      productStuff[prodNum]->productInfo->displayScale);
        check_err(status,__LINE__,__FILE__);
        tmpFloat = productStuff[prodNum]->productInfo->displayMin;
        status = nc_put_att_float(ncid, tmpVarId, "display_min", NC_FLOAT, 1, &tmpFloat);
        check_err(status,__LINE__,__FILE__);
        tmpFloat = productStuff[prodNum]->productInfo->displayMax;
        status = nc_put_att_float(ncid, tmpVarId, "display_max", NC_FLOAT, 1,&tmpFloat);
        check_err(status,__LINE__,__FILE__);
    } // for prodNum
    
    // allocate the data for the biggest size
    fileData = allocateMemory(width*dataSize, "OutFile_netCDF::open fileData");

    // add the quality variable
    if(qualityData) {
        productInfo_t* tmpInfo = allocateProductInfo();
        string qualityName = (string)"qual_" + getProductNameFull(productStuff[0]->productInfo);
        status = findProductInfo(qualityName.c_str(), metaData->sensorID, tmpInfo);
        if(!status) {
            printf("-E- OutFile_netcdf4::open - can not find product %s in product.xml\n", qualityName.c_str());
            exit(EXIT_FAILURE);
        }
        
        status = CreateNCDF(ds_id,
                qualityName.c_str(),    /* short name */
                tmpInfo->description,    /* long name */
                tmpInfo->standardName, /* NetCDF standard name (not set if passed NULL or "") */
                tmpInfo->reference,
                tmpInfo->comment,
                tmpInfo->units,    /* units (not set if passed NULL or "") */
                tmpInfo->validMin,     /* low end of valid range */
                tmpInfo->validMax,    /* high end of range (no range set if low >= high) */
                tmpInfo->scaleFactor,    /* scale factor (not set if 1.0)  */
                tmpInfo->addOffset,   /* scaling offset (not set if 0.0)  */
                tmpInfo->fillValue,       /* fill value */
                NC_BYTE,       /* NCDF number type */
                2,     /* number of dimensions (must be <= 3) */
                dimids /* dimension ids */
        );
        check_err(status,__LINE__,__FILE__);

        status = nc_inq_varid(ncid, qualityName.c_str(), &qualityid);
        check_err(status,__LINE__,__FILE__);

        freeProductInfo(tmpInfo);
    }
    
    // write lat/lon arrays if it is an SMI file
    if(mapProjection.compare("Equidistant Cylindrical") == 0) {
        float*  lonarray = (float*) allocateMemory(getWidth() * sizeof(float), "lonarray");
        float*  latarray = (float*) allocateMemory(getHeight() * sizeof(float), "latarray");
        for (int i=0;i<getWidth();i++){
            lonarray[i] = metaData->west + lonStep * i + lonStep/2.0;
        }
        for (int i=0;i<getHeight();i++){
            latarray[i] = metaData->north - latStep * i - latStep/2.0;
        }

        //Lat
        productInfo_t *p_info;
        p_info = allocateProductInfo();

        status = nc_def_var(ncid, "lat", NC_FLOAT, 1, &dimids[0], &varid1);
        check_err(status,__LINE__,__FILE__);
        if (!findProductInfo("lat", metaData->sensorID, p_info)) {
            printf("lat not found in XML product table\n");
            exit(EXIT_FAILURE);
        }
        status = nc_put_att_text(ncid, varid1, "long_name",
                strlen(p_info->description) + 1,
                p_info->description);

        status = nc_put_att_text(ncid, varid1, "units",
                strlen(p_info->units) + 1,
                p_info->units);

        if (p_info->standardName != NULL) {
            status = nc_put_att_text(ncid, varid1, "standard_name",
                    strlen(p_info->standardName) + 1,
                    p_info->standardName);
        }
        float fv_f32 = (float) p_info->fillValue;
        status = nc_put_att(ncid, varid1, "_FillValue", NC_FLOAT, 1, &fv_f32);
        float valid;
        valid = p_info->validMin;
        status = nc_put_att_float(ncid, varid1, "valid_min", NC_FLOAT, 1,
                &valid);
        valid = p_info->validMax;
        status = nc_put_att_float(ncid, varid1, "valid_max", NC_FLOAT, 1,
                &valid);
        status = nc_put_var(ncid, varid1, latarray);
        check_err(status,__LINE__,__FILE__);

        //Lon
        status = nc_def_var(ncid, "lon", NC_FLOAT, 1, &dimids[1], &varid1);
        check_err(status,__LINE__,__FILE__);
        if (!findProductInfo("lon", metaData->sensorID, p_info)) {
            printf("lon not found in XML product table\n");
            exit(EXIT_FAILURE);
        }
        status = nc_put_att_text(ncid, varid1, "long_name",
                strlen(p_info->description) + 1,
                p_info->description);

        status = nc_put_att_text(ncid, varid1, "units",
                strlen(p_info->units) + 1,
                p_info->units);

        if ( p_info->standardName != NULL) {
            status = nc_put_att_text(ncid, varid1, "standard_name",
                    strlen(p_info->standardName) + 1,
                    p_info->standardName);
        }
        fv_f32 = (float) p_info->fillValue;
        status = nc_put_att(ncid, varid1, "_FillValue", NC_FLOAT, 1, &fv_f32);
        valid = p_info->validMin;
        status = nc_put_att_float(ncid, varid1, "valid_min", NC_FLOAT, 1,
                &valid);
        valid = p_info->validMax;
        status = nc_put_att_float(ncid, varid1, "valid_max", NC_FLOAT, 1,
                &valid);
        status = nc_put_var(ncid, varid1, lonarray);
        check_err(status,__LINE__,__FILE__);

        free(latarray);
        free(lonarray);
        freeProductInfo(p_info);
    }

    /*-----------------------------------------------------------------------*/
    /*  write map_palette */

    // make sure a palette has been loaded
    if(red && green && blue) {
        uint8_t data[768];
        int j=0;
        for (int i = 0; i < 256; i++) {
            data[j++] = red[i];
            data[j++] = green[i];
            data[j++] = blue[i];
        }

        status = nc_def_dim(ncid, "rgb", 3, &dimids[0]);
        check_err(status,__LINE__,__FILE__);

        status = nc_def_dim(ncid, "eightbitcolor", 256, &dimids[1]);
        check_err(status,__LINE__,__FILE__);

        status = nc_def_var(ncid, "palette", NC_UBYTE, 2, dimids, &varid1);
        check_err(status,__LINE__,__FILE__);

        status = nc_put_var(ncid, varid1, data);
        check_err(status,__LINE__,__FILE__);
    }

    return true;
}

void OutFile_netcdf4::writeLine() {
    int status;

    size_t start[3];
    size_t count[3];

    start[0] = currentLine;
    start[1] = 0;
    count[0] = 1;
    count[1] = width;

    for(size_t prodNum=0; prodNum<productStuff.size(); prodNum++) {
        switch(productStuff[prodNum]->dataStorage) {
        case ByteDS:
            for(int i=0; i<width; i++)
                ((int8_t*)fileData)[i] = round(productStuff[prodNum]->calcOutputVal(productStuff[prodNum]->lineData[i]));
            break;
        case UByteDS:
            for(int i=0; i<width; i++)
                ((uint8_t*)fileData)[i] = round(productStuff[prodNum]->calcOutputVal(productStuff[prodNum]->lineData[i]));
            break;
        case ShortDS:
            for(int i=0; i<width; i++)
                ((int16_t*)fileData)[i] = round(productStuff[prodNum]->calcOutputVal(productStuff[prodNum]->lineData[i]));
            break;
        case UShortDS:
            for(int i=0; i<width; i++)
                ((uint16_t*)fileData)[i] = round(productStuff[prodNum]->calcOutputVal(productStuff[prodNum]->lineData[i]));
            break;
        case IntDS:
            for(int i=0; i<width; i++)
                ((int32_t*)fileData)[i] = round(productStuff[prodNum]->calcOutputVal(productStuff[prodNum]->lineData[i]));
            break;
        case UIntDS:
            for(int i=0; i<width; i++)
                ((uint32_t*)fileData)[i] = round(productStuff[prodNum]->calcOutputVal(productStuff[prodNum]->lineData[i]));
            break;
        case FloatDS:
            for(int i=0; i<width; i++)
                ((float*)fileData)[i] = productStuff[prodNum]->lineData[i];
            break;
        case DoubleDS:
            for(int i=0; i<width; i++)
                ((double*)fileData)[i] = productStuff[prodNum]->lineData[i];
            break;
        default:
            printf("-E- OutFile_hdf4::writeLine - unrecognized data type = %d\n", productStuff[prodNum]->dataStorage);
            exit(EXIT_FAILURE);
        }

        status = nc_put_vara(ncid, varid[prodNum], start, count, fileData);
        check_err(status,__LINE__,__FILE__);
    } // for prodNum

    if(qualityData) {
        status = nc_put_vara(ncid, qualityid, start, count, qualityData);
        check_err(status,__LINE__,__FILE__);
    }
        
    currentLine++;
}

bool OutFile_netcdf4::close() {
    int status;

    if(metaData) {
        status = nc_put_att(ncid, NC_GLOBAL, "data_bins", NC_INT, 1, &metaData->data_bins);
        check_err(status,__LINE__,__FILE__);
    }
    float tmpFloat = getFileMinVal();
    status = nc_put_att(ncid, NC_GLOBAL, "data_minimum", NC_FLOAT, 1, &tmpFloat);
    check_err(status,__LINE__,__FILE__);

    tmpFloat = getFileMaxVal();
    status = nc_put_att(ncid, NC_GLOBAL, "data_maximum", NC_FLOAT, 1, &tmpFloat);
    check_err(status,__LINE__,__FILE__);

    status = nc_close(ncid);
    check_err(status,__LINE__,__FILE__);
    ncid = -1;

    varid.clear();
    
    if(fileData) {
        free(fileData);
        fileData = NULL;
    }
    if(qualityData) {
        free(qualityData);
        qualityData = NULL;
    }
    
    return true;
}

int32_t OutFile_netcdf4::addProduct(productInfo_t* productInfo) {
    return addProductNonDisplay(productInfo);
}


