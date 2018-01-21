/*
 * mph_flags.h
 *
 *  Created on: Sep 8, 2015
 *      Author: rhealy
 */

#ifndef SRC_L2GEN_MPH_FLAGS_H_
#define SRC_L2GEN_MPH_FLAGS_H_

#define NMPHFLAGS     3

#define MPH_FLOAT    0x0001       /* floating aquatic vegetation.  */
#define MPH_CYANO    0x0002       /* Cyanobacteria              .  */
#define MPH_ADJ      0x0004       /* Adjacency effect (stray light)*/

#define NHABFLAGS     3
#define HABS_WATER     0x0000       /* Water */
#define HABS_CLOUD     0x0001       /* Cloud Flag */
#define HABS_NONWTR    0x0002       /* Not Water */

static const char *mph_flag_lname[NMPHFLAGS] = {"FLOAT",
                                      "CYANO",
                                      "ADJACENCY"};
static const char *habs_flag_lname[NHABFLAGS] = {"WATER","CLOUD",
                                      "NONWATER"};

#endif /* SRC_L2GEN_MPH_FLAGS_H_ */
