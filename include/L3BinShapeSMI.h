/* 
 * File:   L3BinShape.h
 * Author: dshea
 *
 * Created on July 21, 2015, 8:45 AM
 */

#ifndef L3BINSHAPESMI_H
#define	L3BINSHAPESMI_H

#include <L3BinShape.h>

namespace l3 {

    class L3BinShapeSMI : public L3BinShape {
    protected:
        int32_t totalCols;
        float north, south, east, west;
        
    public:
        L3BinShapeSMI(int32_t numRows, int32_t numCols, 
                float north, float south, float east, float west);
        virtual ~L3BinShapeSMI();
        void constrainRowCol(int32_t &row, int32_t &col) const;

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

