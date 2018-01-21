
#include <stdio.h>
#include <stdint.h>

#include <png.h>
#include <xtiffio.h>
#include <geotiffio.h>

#include <meta_l3b.h>
#include <productInfo.h>
#include <vector>
#include <string>

//earth radius in meters from WGS84 equatorial radius
#define EARTH_RADIUS 6378137.0
#define EARTH_CIRCUMFERENCE 40075016.6856



//---------------------------------------------------
class OutFile {
public:
    enum ColorType { Grayscale, ColorIndex, RGB };
    enum ScaleType { Linear, Log, ArcTan };
    enum DataStorage { ByteDS, UByteDS, ShortDS, UShortDS, IntDS, UIntDS, FloatDS, DoubleDS };
    
    class ProductStuff {
    public:
    
        int32_t width;
        productInfo_t* productInfo;
        DataStorage dataStorage; //< data type for the file storage
        ScaleType scaleType;  //< display Linear, Log, ATan scaling
        double scale;        //< display slope for scaling.  If log it is actually log10(slope)
        double offset;       //< display offset for scaling.
        double minOutputVal; //< display min output value
        double maxOutputVal; //< display max output value
        double minVal;       //< display min physical value
        double maxVal;       //< display max physical value
        double missingValue; //< missing value from product XML (val stored in file)
        double* lineData;

        ProductStuff(int32_t width, const productInfo_t* productInfo);
        ProductStuff(const OutFile::ProductStuff& pStuff);
        ~ProductStuff();
        
        int32_t getWidth() { return width; };
        void setScale(double min, double max, ScaleType scaleType);
        void setScale(double min, double max, ScaleType scaleType,
            double minOutput, double maxOutput);
        void setScaleOffset(double scale, double offset, ScaleType scaleType);
        void setScaleOffset(double scale, double offset, ScaleType scaleType,
            double minOutput, double maxOutput);
        double calcOutputVal(double val) const;
        double calcPhysicalVal(double val) const;
    };
    

protected:
    
    static const uint8_t qualityUnused = 255;
    static const double badPixelValue = -32767.0;

    std::string fileName;
    int32_t width;
    int32_t height;
    uint8_t* qualityData;
    uint32_t currentLine; //< current line number (0 based)

    ColorType colorType;

    double fileMinVal;      //< current min value of data written to file
    double fileMaxVal;      //< current max value of data written to file
    double resolution;   //< geospatial resolution of a pixel in the center of the scene (meters)
    int deflate;         //< compression setting for netCDF files

    uint8_t* red;
    uint8_t* green;
    uint8_t* blue;
    meta_l3bType* metaData;
    std::string mapProjection;
    
    std::vector<ProductStuff*> productStuff;

    OutFile();
    virtual ~OutFile();
    virtual int addProductNonDisplay(productInfo_t* productInfo);

public:
    virtual void setSize(int32_t width, int32_t height);
    virtual int32_t getWidth() const;
    virtual int32_t getHeight() const;
    virtual void setFileName(std::string fileName);
    virtual std::string getFileName() { return fileName; };
    
    /**
     * Open the output file. Make sure you have set these functions first:
     *      setSize, setScale, setPalette, setMetaData, setProductInfo, setQualityProcessing.
     * @param fileName
     * @return true if successful.
     */
    virtual bool open() = 0;
    virtual bool close() = 0;
    virtual double getMinValue(int32_t prod=0) { return productStuff[prod]->minVal; };
    virtual double getMaxValue(int32_t prod=0) { return productStuff[prod]->maxVal; };
    virtual void setPixel(int32_t x, double val, int32_t prod=0);
    virtual void setPixelRGB(int32_t x, float red, float green, float blue);
    virtual void setQuality(int32_t x, uint8_t val);
    virtual void fillPixel(int32_t x);
    virtual void missingPixel(int32_t x);
    virtual void writeLine() = 0;
    virtual bool setPalette(const char* paletteName);
    virtual void setMetaData(meta_l3bType* metaData);
    virtual meta_l3bType* getMetadata() { return metaData; };
    virtual int32_t addProduct(productInfo_t* productInfo);
    virtual int32_t getNumProducts() { return productStuff.size(); };
    virtual void setMapProjection(std::string projection);
    virtual void setNumFilledPixels(int32_t num);
    virtual int32_t getNumFilledPixels();
    virtual float getPercentFilledPixels();
    virtual void resetFileMinMax();
    virtual double getFileMinVal() { return fileMinVal; };
    virtual double getFileMaxVal() { return fileMaxVal; };
    virtual void setResolution(std::string resolutionStr);
    virtual void setResolution(double resolution) { this->resolution = resolution; };
    virtual double getResolution() { return resolution; };
    virtual void setQualityProcessing(bool val);
    virtual bool getQualityProcessing();
    virtual void setDeflate(int val) { deflate = val; };
    virtual int getDeflate() { return deflate; };
};

//---------------------------------------------------
class OutFile_pgm : public OutFile {
protected:
    FILE *outfp;
    uint8_t* fileData;

public:
    OutFile_pgm();
    virtual ~OutFile_pgm();
    virtual void setSize(int32_t width, int32_t height);
    virtual bool open();
    virtual bool close();
    virtual void writeLine();
};

//---------------------------------------------------
class OutFile_ppm : public OutFile_pgm {

public:
    virtual void setSize(int32_t width, int32_t height);
    virtual bool open();
    virtual void writeLine();
};

//---------------------------------------------------
class OutFile_ppm_rgb : public OutFile {
    FILE *outfp;
    uint8_t* fileData;

public:
    OutFile_ppm_rgb();
    virtual ~OutFile_ppm_rgb();
    virtual void setSize(int32_t width, int32_t height);
    virtual bool open();
    virtual bool close();
    virtual void setPixel(int32_t x, double val, int32_t prod=0);
    virtual void setPixelRGB(int32_t x, float red, float green, float blue);
    virtual void fillPixel(int32_t x);
    virtual void missingPixel(int32_t x);
    virtual void writeLine();
};

//---------------------------------------------------
class OutFile_png : public OutFile {
    FILE *outfp;
    uint8_t* fileData;
    bool isColor;
    png_structp png_ptr;
    png_infop info_ptr;

public:
    OutFile_png(bool color);
    virtual ~OutFile_png();
    virtual void setSize(int32_t width, int32_t height);
    virtual bool open();
    virtual bool close();
    virtual void writeLine();
};

//---------------------------------------------------
class OutFile_png_rgb : public OutFile {
    FILE *outfp;
    uint8_t* fileData;
    png_structp png_ptr;
    png_infop info_ptr;

public:
    OutFile_png_rgb();
    virtual ~OutFile_png_rgb();
    virtual void setSize(int32_t width, int32_t height);
    virtual bool open();
    virtual bool close();
    virtual void setPixel(int32_t x, double val, int32_t prod=0);
    virtual void setPixelRGB(int32_t x, float red, float green, float blue);
    virtual void fillPixel(int32_t x);
    virtual void missingPixel(int32_t x);
    virtual void writeLine();
};

//---------------------------------------------------
class OutFile_tiff_color : public OutFile {
    uint8_t* fileData;
    TIFF *tiff;
    GTIF *gtif;

public:
    OutFile_tiff_color();
    virtual ~OutFile_tiff_color();
    virtual void setSize(int32_t width, int32_t height);
    virtual bool open();
    virtual bool close();
    virtual void writeLine();
};

//---------------------------------------------------
class OutFile_tiff_gray : public OutFile {
    float* fileData;
    TIFF *tiff;
    GTIF *gtif;

public:
    OutFile_tiff_gray();
    virtual ~OutFile_tiff_gray();
    virtual void setSize(int32_t width, int32_t height);
    virtual bool open();
    virtual bool close();
    virtual void writeLine();
};

//---------------------------------------------------
class OutFile_tiff_rgb : public OutFile {
    uint8_t* fileData;
    TIFF *tiff;
    GTIF *gtif;

public:
    OutFile_tiff_rgb();
    virtual ~OutFile_tiff_rgb();
    virtual void setSize(int32_t width, int32_t height);
    virtual bool open();
    virtual bool close();
    virtual void setPixel(int32_t x, double val, int32_t prod=0);
    virtual void setPixelRGB(int32_t x, float red, float green, float blue);
    virtual void fillPixel(int32_t x);
    virtual void missingPixel(int32_t x);
    virtual void writeLine();
};

//---------------------------------------------------
class OutFile_hdf4 : public OutFile {
    void* fileData;
    int32 sdfid;  ///< HDF SD file ID
    int32 sdsid;  ///< HDF SDS dataset ID
    int32 quality_sdsid;  ///< HDF SDS dataset ID for SST quality factor
    int32 hdfDataType;

public:
    OutFile_hdf4();
    virtual ~OutFile_hdf4();
    virtual void setSize(int32_t width, int32_t height);
    virtual bool open();
    virtual bool close();
    virtual void writeLine();
    virtual int addProduct(productInfo_t* productInfo);
};

//---------------------------------------------------
class OutFile_netcdf4 : public OutFile {
    void* fileData;
    int ncid;   ///< NetCDF file ID
    std::vector<int> varid;  ///< NetCDF variable IDs
    int qualityid; ///< NetCDF variable ID for quality factor

public:
    OutFile_netcdf4();
    virtual ~OutFile_netcdf4();
    virtual void setSize(int32_t width, int32_t height);
    virtual bool open();
    virtual bool close();
    virtual void writeLine();
    virtual int addProduct(productInfo_t* productInfo);
};
