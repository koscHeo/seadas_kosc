/* 
 * File:   setupflags.h
 * Author: dshea
 *
 * Created on January 31, 2012, 8:47 AM
 */

#ifndef SETUPFLAGS_H
#define	SETUPFLAGS_H

#include <hdf.h>


#ifdef	__cplusplus
extern "C" {
#endif

void setupflags (char *flagsdef, char *flaguse, uint32 *flagusemask, 
        uint32 *required, int *status);



#ifdef	__cplusplus
}
#endif

#endif	/* SETUPFLAGS_H */
