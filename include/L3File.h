#ifndef L3File_h
#define L3File_h

#include <hdf_bin.h>
#include <vector>
#include <deque>
#include <list>
#include <stdint.h>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/foreach.hpp>

#include <L3BinShape.h>

namespace l3 {
    const float missingPixelValue = -32767.0; 
    
    typedef boost::geometry::model::d2::point_xy<float> Point_t;
    typedef boost::geometry::model::polygon<Point_t> Polygon_t;
    typedef boost::geometry::model::box<Point_t> Box_t;


    template<typename Geometry>
    inline float calcIntersectionArea(Box_t &box, Geometry &geo) {
        std::deque<Geometry> outPolygons;
        boost::geometry::intersection(box, geo, outPolygons);
        float theArea = 0.0;
        BOOST_FOREACH(Geometry const& p, outPolygons)
        {
            theArea += boost::geometry::area(p);
        }
        return theArea;
    }

    template<>
    inline float calcIntersectionArea<Box_t>(Box_t &box, Box_t &geo) {
        Box_t outBox;
        boost::geometry::intersection(box, geo, outBox);
        return boost::geometry::area(outBox);
    }

    class L3Row;
    class L3File;

    /**
     * class used to return all of the information about a bin
     */
    class L3Bin {
        friend class L3Row;
        friend class L3File;
        friend class L3FileSMI;

    private:
        static const uint8_t qualityUnused = 255;

        void internalAllocate(int32_t numProducts);
        void internalCopy(const L3Bin &bin);
        void internalAdd(const L3Bin &bin);
        void internalInit(int32_t numProducts);
        void checkProdId(int32_t prodId) const;

    protected:
        int32_t numProducts;
        int32_t binNum;
        int32_t recordNum;
        int32_t nobs;
        int32_t nscenes;
        uint8_t quality;
        float timeRec; //< timeRec/nobs = observation time (sec since 1993 TAI)
        float weights;
        float *sums;
        float *sumSquares;

    public:
        L3Bin(int32_t numProducts=1);
        L3Bin(const L3Bin &bin);
        ~L3Bin();

        void clear();

        L3Bin& operator=(const L3Bin &bin);
        L3Bin& operator+=(const L3Bin &bin);

        void addWeighted(const L3Bin &bin, float weighting);

        template<typename Geometry1, typename Geometry2>
        inline void addWeighted(const L3Bin &bin, Geometry1 &box, Geometry2 &geo) {
            float theArea = calcIntersectionArea(box, geo);
            if(theArea > 0.0)
                addWeighted(bin, theArea);
        }

        void setNumProducts(int32_t numProducts);
        int32_t getNumProducts() const;
        int32_t getBinNum() const;
        int32_t getNobs() const;
        int32_t getNscenes() const;
        int32_t getRecordNum() const;
        float getObsTime() const;
        float getWeights() const;
        uint8_t getQuality() const;

        float getSum(int32_t prodId=0) const;
        float getSumSquares(int32_t prodId=0) const;

        float getMean(int32_t prodId=0) const;
        float getVariance(int32_t prodId=0) const;
        float getStdev(int32_t prodId=0) const;
        
    };

    /**
     * class to hold the information for one row of a L3 bin file
     */
    class L3Row {
        friend class L3File;
        friend class L3FileSMI;

    protected:
        int32_t row;
        int32_t numProducts;
        int32_t numBins;
        std::vector<L3Bin*> binArray;
        int32_t lastBin; //< bin index found by getBin last time

    public:
        L3Row(int32_t row, int32_t numBins, int32_t numProducts);
        ~L3Row();

        int32_t getRow() const;
        int32_t getNumProducts() const;
        int32_t getNumBins() const;
        void setNumBins(int32_t numBins);
        L3Bin* getBin(int32_t binNum);
        L3Bin* getBinByIndex(int32_t index);
    };

    /**
     * class to read the information out of a L3 bin file
     */
    class L3File {

    protected:
        L3BinShape* shape;
        Hdf::hdf_bin* binObj;
        std::list<L3Row*> rowList;
        int numCacheRows;
        int32_t *baseRecord; ///< first record index of each row (0 based)
        int32_t *extentbin; ///< save the number of records in each row
        float* sumBuffer; ///< buffer for the sums (sum, sumSq, 2*nrows, numProd)
        L3Bin outBin; ///< accumulation bin for lookups
        uint8_t* qualityBuffer; ///< buffer for the quality (2*nrows)
        size_t *prodMap; ///< index map used to map my product index into binObj's index
        std::vector<string> activeProdNameList; ///< array of active product names
        float fudgeFactor; ///< fudge factor used to change the L3 pixel size
        
        virtual L3Row* readRow(int32_t row);

        virtual int initRecordLookup();

        /**
         * Search for bins in this row that are inside the geometry.  Add the bins 
         * found into class member outBin
         * @param lat0 latitude to start at
         * @param lon0 longitude to start at
         * @param geo geometry to intersect with
         * @return true if a bin in this row intersected the geometry
         */
        template<typename Geometry>
        inline bool addBinsFromRow(float lat0, float lon0, Geometry geo, bool areaWeighted) {        
            bool foundIntersect = false;
            int32 row0, col0;
            shape->latlon2rowcol(lat0, lon0, row0, col0);

            int32 col = col0;
            float lat, lon;

            float deltaLat = 90.0 / shape->getNumRows() * fudgeFactor;
            float deltaLon = 180.0 / shape->getNumCols(row0) * fudgeFactor;

            // look right
            bool binIntersects = true;
            while(binIntersects) {
                if(col > shape->getNumCols(row0))
                    break;
                shape->rowcol2latlon(row0, col, lat, lon);
                Point_t pMin(lon-deltaLon, lat-deltaLat);
                Point_t pMax(lon+deltaLon, lat+deltaLat);
                Box_t box(pMin, pMax);
                binIntersects = boost::geometry::intersects(geo, box);
                if(binIntersects) {
                    foundIntersect = true;
                    L3Bin* tmpBin = getBin(row0, col);
                    if(tmpBin) {
                        if(areaWeighted) {
                            outBin.addWeighted(*tmpBin, box, geo);
                        } else
                            outBin += *tmpBin;
                    }
                    col++;
                }
            }

           // look left
            col = col0 - 1;
            binIntersects = true;
            while(binIntersects) {
                if(col < 0)
                    break;
                shape->rowcol2latlon(row0, col, lat, lon);
                Point_t pMin(lon-deltaLon, lat-deltaLat);
                Point_t pMax(lon+deltaLon, lat+deltaLat);
                Box_t box(pMin, pMax);
                binIntersects = boost::geometry::intersects(geo, box);
                if(binIntersects) {
                    foundIntersect = true;
                    L3Bin* tmpBin = getBin(row0, col);
                    if(tmpBin) {
                        if(areaWeighted) {
                            outBin.addWeighted(*tmpBin, box, geo);
                        } else
                            outBin += *tmpBin;
                    }
                    col--;
                }
            }

            return foundIntersect;
        };

    public:
        L3File();
        virtual ~L3File();

        virtual int32_t rowbin2record(int32_t row, int32_t bin);
        virtual int32_t rowcol2record(int32_t row, int32_t col);
        virtual int32_t latlon2record(float lat, float lon);
        virtual int32_t bin2record(int32_t bin);

        virtual void clearCache();
        virtual void setNumCacheRows(int32_t numRows);
        virtual bool open(const char* fileName);
        virtual void close();
        virtual meta_l3bType* getMetaData();
        virtual int32_t getNumProducts();
        virtual std::string getProductName(size_t index=0);
        virtual bool setActiveProductList(const char* prodStr);
        virtual int32_t getNumActiveProducts();
        virtual std::string getActiveProductName(size_t index=0);
        virtual int32_t getNumRows();
        virtual L3Row* getRow(int32_t row);
        virtual L3Bin* getBin(int32 row, int32 col);
        virtual L3Bin* getClosestBin(float lat, float lon);

        /**
         * return a bin with binned data for all of the bins that intersect
         * Geometry 
         * @param geo boost geometry to define an area
         * @return pointer to a bin containing all of the binned data found in this 
         * area, or NULL if no bins found. This pointer is good until the next 
         * call to this function.
         */
        template<typename Geometry>
        inline L3Bin* getBinsInside(Geometry geo, bool areaWeighted=false) {
            L3Bin* bin;

            outBin.clear();

            // get the center of the geometry
            Point_t geoCenter;
            boost::geometry::centroid(geo, geoCenter);

            // find r,c of closest bin
            float geoLat = constrainLat(geoCenter.y());
            float geoLon = constrainLon(geoCenter.x());
            geoCenter.y(geoLat);
            geoCenter.x(geoLon);

            int32 geoRow, geoCol; // note row,col are 0 based
            float tmpFloat;
            shape->latlon2rowcol(geoLat, geoLon, geoRow, geoCol);
            shape->rowcol2latlon(geoRow, geoCol, geoLat, tmpFloat); // adjust lat on the row center

            float deltaLat = 180.0 / shape->getNumRows();
            float lat = geoLat;
            bool binFound = true;

            // look up
            while(binFound) {
                if(lat > 90.0)
                    break;
                binFound = addBinsFromRow(lat, geoLon, geo, areaWeighted);
                lat = lat + deltaLat;
            }

            // look down
            binFound = true;
            lat = geoLat - deltaLat;
            while(binFound) {
                if(lat < -90.0)
                    break;
                binFound = addBinsFromRow(lat, geoLon, geo, areaWeighted);
                lat = lat - deltaLat;
            }

            if(outBin.nobs == 0)
                return NULL;
            else
                return &outBin;
        };

        virtual Hdf::hdf_bin* getHdfBinObject() const;
        virtual bool hasQuality();
        virtual void setQualityProcessing(bool val);
        virtual bool getQualityProcessing() const;
        virtual L3BinShape* getShape() const { return shape; };
        virtual void setFudgeFactor(float fudge) { fudgeFactor = fudge; };
        virtual float getFudgeFactor() { return fudgeFactor; };
        
    };

}; // namespace l3

#endif

