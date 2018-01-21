/*
 * l2extract.h
 *
 *  Created on: May 6, 2014
 *      Author: dshea
 */

#ifndef L2EXTRACT_H_
#define L2EXTRACT_H_


int extractNetCDF(const char* infile, const char* outfile, int spix, int epix,
        int sscan, int escan, const char* prodlist);



#endif /* L2EXTRACT_H_ */
