#include <L3BinShapeSMI.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


namespace l3 {


    /**
     * constructor for a SMI shaped l3 bin file
     * @param numRows number of rows in the bin file
     */
    L3BinShapeSMI::L3BinShapeSMI(int32_t numRows, int32_t numCols, 
                float north, float south, float east, float west) 
    : L3BinShape(numRows) {
        totalCols = numCols;
        totalBins = totalRows * totalCols;
        
        this->north = north;
        this->south = south;
        this->east = east;
        this->west = west;
    }

    /**
     * destructor that frees the memory for this class
     */
    L3BinShapeSMI::~L3BinShapeSMI() {
    }

    /**
     * keep row and col in the valid range
     * @param row to modify if necessary
     * @param col to modify if necessary
     */
    void L3BinShapeSMI::constrainRowCol(int32_t &row, int32_t &col) const {
        if (row < 0)
            row = 0;
        else if (row >= totalRows)
            row = totalRows - 1;

        if (col < 0)
            col = 0;
        else if (col >= totalCols)
            col = totalCols - 1;
    }

    /**
     * get the bin number of the first column in the given row
     * @param row to look up
     * @return the first bin in this row
     */
    int32_t L3BinShapeSMI::getBaseBin(int32_t row) const {
        constrainRow(row);
        return row * totalCols + 1;
    }

    /**
     * get the number of columns in the given row
     * @param row to look up
     * @return the number of cols in this row
     */
    int32_t L3BinShapeSMI::getNumCols(int32_t row) const {
        return totalCols;
    }

    /**
     * find the row that contains the given bin number
     * @param bin to look up
     * @return row this bin is in
     */
    int32_t L3BinShapeSMI::bin2row(int32_t bin) {
        if (bin < 1)
            return 0; /* south pole */
        if (bin >= totalBins)
            return totalRows - 1; /* north pole */
        return (bin-1) / totalCols;
    }

    /**
     * get the row and col for the given bin
     * @param bin to look up
     * @param row that contains the bin
     * @param col for the bin
     */
    void L3BinShapeSMI::bin2rowcol(int32_t bin, int32_t &row, int32_t &col) {
        if (bin < 1) {
            row = col = 0; /* south pole */
            return;
        }
        if (bin >= totalBins) {
            row = totalRows - 1;
            col = totalCols - 1; /* north pole */
            return;
        }
        bin--;
        row = bin / totalCols;
        col = bin - row*totalCols;
    }

    /**
     * get the bin at the given row/col
     * @param row of the bin
     * @param col of the bin
     * @return bin number at row/col
     */
    int32_t L3BinShapeSMI::rowcol2bin(int32_t row, int32_t col) const {
        constrainRowCol(row, col);
        return row*totalCols + col + 1;
    }

    /**
     * get the center lat/lon for the given row/col
     * @param row of bin
     * @param col of bin
     * @param lat center lat
     * @param lon center lon
     */
    void L3BinShapeSMI::rowcol2latlon(int32_t row, int32_t col, float &lat, float &lon) const {
        constrainRowCol(row, col);
        lat = (north - south) * (row + 0.5) / totalRows + south;
        lon = (east - west) * (col + 0.5) / totalCols + west;
    }

    /**
     * get the row for the given lat
     * @param lat to look up
     * @return row of containing bin
     */
    int32_t L3BinShapeSMI::lat2row(float lat) const {
        lat = constrainLat(lat);
        if(lat>north || lat<south)
            return -1;

        int32_t row = (int32_t) ((lat - south) * (float) totalRows / (north - south));
        if (row >= totalRows)
            row = totalRows - 1;
        if(row < 0)
            row = 0;
        return row;
    }
    
    /**
     * get the row/col that contains the given lat/lon
     * @param lat to look up
     * @param lon to look up
     * @param row of containing bin
     * @param col of containing bin
     */
    void L3BinShapeSMI::latlon2rowcol(float lat, float lon, int32_t &row, int32_t &col) const {
        row = lat2row(lat);
        if(row == -1) {
            col = -1;
            return;
        }
        lon = constrainLon(lon);
        if(lon < west || lon > east) {
            row = col = -1;
            return;
        }
        col = (int32_t) ((float) totalCols * (lon - west) / (east - west));
        if (col >= totalCols || col < 0)
            row = col = -1;
    }

    /**
     * get the bin number containing the given lat/lon
     * @param lat to look up
     * @param lon to look up
     * @return  bin number containing lat/lon
     */
    int32_t L3BinShapeSMI::latlon2bin(float lat, float lon) const {
        int32_t row, col;

        latlon2rowcol(lat, lon, row, col);
        return row * totalCols + col;
    }

    /**
     * get the boundaries of the bin at row/col
     * @param row of bin
     * @param col of bin
     * @param north extent of bin
     * @param south extent of bin
     * @param east extent of bin
     * @param west extent of bin
     */
    void L3BinShapeSMI::rowcol2bounds(int32_t row, int32_t col, float &north, float &south, float &east, float &west) const {
        float lat, lon;
        rowcol2latlon(row, col, lat, lon);
        float tmp = (north - south) / totalRows / 2;
        north = lat + tmp;
        south = lat - tmp;
        
        tmp = (east - west) / totalCols / 2;
        east = lon + tmp;
        west = lon - tmp;
    }
    
}