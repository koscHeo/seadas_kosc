#include <L3BinShape.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// this should not be necessary, but M_PI seems to not be defined in C99
#define PI 3.141592653589793


namespace l3 {

    /**
     * limit latitude values to +-90
     * @param lat latitude to normalize
     * @return lat value clamped to -90 to 90
     */
    double constrainLat(double lat) {
        if (lat > 90.0)
            return 90.0;
        if (lat < -90.0)
            return -90.0;
        return lat;
    }

    /**
     * limit the longitude to +-180
     * @param lon longitude to normalize
     * @return lat value normalized to -180 to 180
     */
    double constrainLon(double lon) {
        int count = 0;
        while (lon < -180.0) {
            lon += 360.0;
            count++;
            if (count > 10) {
                fprintf(stderr, "-E- %s %d: constrainLon : longitude=%f way out of range.\n", __FILE__,
                        __LINE__, lon);
                exit(EXIT_FAILURE);
            }
        }

        count = 0;
        while (lon > 180.0) {
            lon -= 360.0;
            count++;
            if (count > 10) {
                fprintf(stderr, "-E- %s %d: constrainLon : longitude=%f way out of range.\n", __FILE__,
                        __LINE__, lon);
                exit(EXIT_FAILURE);
            }
        }
        return lon;
    }

    /**
     * constructor for a integerized sine shaped l3 bin file
     * @param numRows number of rows in the bin file
     */
    L3BinShape::L3BinShape(int32_t numRows) {
        seamLon = -180.0; /*  this value should be passed in  */
        totalRows = numRows;
        totalBins = 0;
    }

    /**
     * destructor that frees the memory for this class
     */
    L3BinShape::~L3BinShape() {
    }

    /**
     * get the number of rows in this shape
     * @return number of rows
     */
    int32_t L3BinShape::getNumRows() const {
        return totalRows;
    }
    
    /**
     * keep the row within 0 to number of rows
     * @param row to change if necessary
     */
    void L3BinShape::constrainRow(int32_t &row) const {
        if (row < 0)
            row = 0;
        else if (row >= totalRows)
            row = totalRows - 1;
    }

    /**
     * get the center lat/lon for the given bin number
     * @param bin to look up
     * @param lat center lat
     * @param lon center lon
     */
    void L3BinShape::bin2latlon(int32_t bin, float &lat, float &lon) {
        int32_t row, col;

        bin2rowcol(bin, row, col);
        rowcol2latlon(row, col, lat, lon);
    }

    /**
     * get the boundaries of the bin number
     * @param bin to look up
     * @param north extent of bin
     * @param south extent of bin
     * @param east extent of bin
     * @param west extent of bin
     */
    void L3BinShape::bin2bounds(int32_t bin, float &north, float &south, float &east, float &west) {
        int32_t row, col;
        bin2rowcol(bin, row, col);
        rowcol2bounds(row, col, north, south, east, west);
    }

    
    /**
     * constructor for a integerized sine shaped l3 bin file
     * @param numRows number of rows in the bin file
     */
    L3BinShapeInst::L3BinShapeInst(int32_t numRows) : L3BinShape(numRows) {
        int32_t i;
        double radfac = PI / 180.0;

        oldRow = 0;

        numBin = (int32_t *) calloc(totalRows, sizeof (int32_t));
        baseBin = (int32_t *) calloc(totalRows + 1, sizeof (int32_t));
        latBin = (float *) calloc(totalRows, sizeof (float));

        for (i = 0; i < totalRows; i++) {
            latBin[i] = (i + 0.5) * (180.0 / totalRows) - 90.0;
            numBin[i] = (int32_t) (cos(latBin[i] * radfac) * (2.0 * totalRows) + 0.5);
        }

        baseBin[0] = 1;

        for (i = 1; i <= totalRows; i++)
            baseBin[i] = baseBin[i - 1] + numBin[i - 1];
        totalBins = baseBin[totalRows] - 1;
    }

    /**
     * destructor that frees the memory for this class
     */
    L3BinShapeInst::~L3BinShapeInst() {
        free(numBin);
        free(baseBin);
        free(latBin);
    }

    /**
     * keep row and col in the valid range
     * @param row to modify if necessary
     * @param col to modify if necessary
     */
    void L3BinShapeInst::constrainRowCol(int32_t &row, int32_t &col) const {
        if (row < 0)
            row = 0;
        else if (row >= totalRows)
            row = totalRows - 1;

        if (col < 0)
            col = 0;
        else if (col >= numBin[row])
            col = numBin[row] - 1;
    }

    /**
     * get the bin number of the first column in the given row
     * @param row to look up
     * @return the first bin in this row
     */
    int32_t L3BinShapeInst::getBaseBin(int32_t row) const {
        constrainRow(row);
        return baseBin[row];
    }

    /**
     * get the number of columns in the given row
     * @param row to look up
     * @return the number of cols in this row
     */
    int32_t L3BinShapeInst::getNumCols(int32_t row) const {
        constrainRow(row);
        return numBin[row];
    }

    /**
     * find the row that contains the given bin number
     * @param bin to look up
     * @return row this bin is in
     */
    int32_t L3BinShapeInst::bin2row(int32_t bin) {
        int32_t rlow, rhi, rmid;

        if (baseBin[oldRow] <= bin && baseBin[oldRow + 1] > bin) {
            return oldRow;
        } else {
            if (bin < 1)
                return 0; /* south pole */
            if (bin >= totalBins)
                return totalRows - 1; /* north pole */

            /* binary search for row */
            rlow = 0;
            rhi = totalRows - 1;
            while (1) {
                rmid = (rlow + rhi + 1) / 2;
                if (baseBin[rmid] > bin)
                    rhi = rmid - 1;
                else
                    rlow = rmid;

                if (rlow == rhi) {
                    oldRow = rlow;
                    return rlow;
                }
            }
        }
        return 0;
    }

    /**
     * get the row and col for the given bin
     * @param bin to look up
     * @param row that contains the bin
     * @param col for the bin
     */
    void L3BinShapeInst::bin2rowcol(int32_t bin, int32_t &row, int32_t &col) {
        if (bin < 1) {
            row = col = 0; /* south pole */
            return;
        }
        if (bin >= totalBins) {
            row = totalRows - 1;
            col = numBin[row] - 1; /* north pole */
            return;
        }
        row = bin2row(bin);
        col = bin - baseBin[row];
    }

    /**
     * get the bin at the given row/col
     * @param row of the bin
     * @param col of the bin
     * @return bin number at row/col
     */
    int32_t L3BinShapeInst::rowcol2bin(int32_t row, int32_t col) const {
        constrainRowCol(row, col);
        return baseBin[row] + col;
    }

    /**
     * get the center lat/lon for the given row/col
     * @param row of bin
     * @param col of bin
     * @param lat center lat
     * @param lon center lon
     */
    void L3BinShapeInst::rowcol2latlon(int32_t row, int32_t col, float &lat, float &lon) const {
        constrainRowCol(row, col);
        lat = latBin[row];
        lon = 360.0 * (col + 0.5) / numBin[row] + seamLon;
    }

    /**
     * get the row for the given lat
     * @param lat to look up
     * @return row of containing bin
     */
    int32_t L3BinShapeInst::lat2row(float lat) const {
        lat = constrainLat(lat);
        int32_t row = (int32_t) ((90.0 + lat) * (float) totalRows / 180.0);
        if (row >= totalRows)
            row = totalRows - 1;
        return row;
    }
    
    /**
     * get the row/col that contains the given lat/lon
     * @param lat to look up
     * @param lon to look up
     * @param row of containing bin
     * @param col of containing bin
     */
    void L3BinShapeInst::latlon2rowcol(float lat, float lon, int32_t &row, int32_t &col) const {
        row = lat2row(lat);
        lon = constrainLon(lon);
        col = (int32_t) ((float) numBin[row] * (lon - seamLon) / 360.0);
        if (col >= numBin[row])
            col = numBin[row] - 1;
    }

    /**
     * get the bin number containing the given lat/lon
     * @param lat to look up
     * @param lon to look up
     * @return  bin number containing lat/lon
     */
    int32_t L3BinShapeInst::latlon2bin(float lat, float lon) const {
        int32_t row, col;

        latlon2rowcol(lat, lon, row, col);
        return baseBin[row] + col;
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
    void L3BinShapeInst::rowcol2bounds(int32_t row, int32_t col, float &north, float &south, float &east, float &west) const {
        float lat, lon;
        rowcol2latlon(row, col, lat, lon);
        north = lat + 90.0 / totalRows;
        south = lat - 90.0 / totalRows;
        east = lon + 180.0 / numBin[row];
        west = lon - 180.0 / numBin[row];
    }


    
    
    
    
    
}