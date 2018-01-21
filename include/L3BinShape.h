/* 
 * File:   L3BinShape.h
 * Author: dshea
 *
 * Created on July 21, 2015, 8:45 AM
 */

#ifndef L3BINSHAPE_H
#define	L3BINSHAPE_H

#include <stdint.h>


namespace l3 {

    double constrainLat(double lat);
    double constrainLon(double lon);
    
    
    class L3BinShape {
    protected:
        int32_t totalBins; // total number of bins in the L3 bin shape
        int32_t totalRows; // total number of rows in the L3 bin shape
        float seamLon; // longitude of the start of the shape

    public:
        L3BinShape(int32_t numRows);
        virtual ~L3BinShape();

        virtual int32_t getNumRows() const;
        virtual void constrainRow(int32_t &row) const;
        virtual void constrainRowCol(int32_t &row, int32_t &col) const = 0;
        virtual int32_t getBaseBin(int32_t row) const = 0;
        virtual int32_t getNumCols(int32_t row) const = 0;

        virtual int32_t bin2row(int32_t bin) = 0;
        virtual void bin2rowcol(int32_t bin, int32_t &row, int32_t &col) = 0;
        virtual int32_t rowcol2bin(int32_t row, int32_t col) const = 0;
        virtual void rowcol2latlon(int32_t row, int32_t col, float &lat, float &lon) const = 0;
        virtual void bin2latlon(int32_t bin, float &lat, float &lon);
        virtual int32_t lat2row(float lat) const = 0;
        virtual void latlon2rowcol(float lat, float lon, int32_t &row, int32_t &col) const = 0;
        virtual int32_t latlon2bin(float lat, float lon) const = 0;
        virtual void rowcol2bounds(int32_t row, int32_t col, float &north, float &south, float &east, float &west) const = 0;
        virtual void bin2bounds(int32_t bin, float &north, float &south, float &east, float &west);
    };

    
    class L3BinShapeInst : public L3BinShape {
    protected:
        int32_t *numBin; // number of bins in each row
        int32_t *baseBin; // 1st bin of each row
        float *latBin; // center latitude of each row

        int32_t oldRow; // row that was found on last search

    public:
        L3BinShapeInst(int32_t numRows);
        virtual ~L3BinShapeInst();

        virtual void constrainRowCol(int32_t &row, int32_t &col) const;
        virtual int32_t getBaseBin(int32_t row) const;
        virtual int32_t getNumCols(int32_t row) const;

        virtual int32_t bin2row(int32_t bin);
        virtual void bin2rowcol(int32_t bin, int32_t &row, int32_t &col);
        virtual int32_t rowcol2bin(int32_t row, int32_t col) const;
        virtual void rowcol2latlon(int32_t row, int32_t col, float &lat, float &lon) const;
        virtual int32_t lat2row(float lat) const;
        virtual void latlon2rowcol(float lat, float lon, int32_t &row, int32_t &col) const;
        virtual int32_t latlon2bin(float lat, float lon) const;
        virtual void rowcol2bounds(int32_t row, int32_t col, float &north, float &south, float &east, float &west) const;
    };

}

#endif	/* L3BINSHAPE_H */

