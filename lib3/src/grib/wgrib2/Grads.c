#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grb2.h"
#include "wgrib2.h"
#include "fnlist.h"
#define RINT(x)         ((int) (floor((x) + 0.5)))

/*
 * here are options that are used g2ctl
 * (making GrADS control files for grib2 files
 *
 * this file can be deleted if not needed
 *
 *   This file is in the public domain  2007 Wesley Ebisuzaki
 */

extern enum input_type input;
extern int header, dump_rec, dump_submsg;
extern int mode, decode, latlon, ndata;
extern float *lat, *lon;

extern int file_append;
extern char *item_deliminator;


/*
 * HEADER:100:ctl_inv:inv:0:ctl inventory dump (for g2ctl)
 */
int f_ctl_inv(ARG0) {

    char name[STRING_SIZE], desc[STRING_SIZE], unit[STRING_SIZE];
    float val1, val2;
    int level_type1, level_type2;
    int undef_val1, undef_val2, i;
    char *string;

    if (mode == -1) return 0;
    if (mode == -2) return 0;

    if (getName(sec, mode, NULL, name, desc, unit) != 0) {
	return -1;
    }
    level_type1 = sec[4][22];
    undef_val1 = undef_val2 = 1;
    val1 = val2 = -999;

    if (sec[4][23] != 255 || sec[4][24] != 255 || sec[4][25] != 255 || sec[4][26] != 255 ||
                sec[4][24] != 255) {
        undef_val1 = 0;
        val1 = scaled2flt(INT1(sec[4][23]), int4(sec[4] + 24));
    }

    level_type2 = sec[4][28];
    if (sec[4][29] != 255 || sec[4][30] != 255 || sec[4][31] != 255 || sec[4][32] != 255 ||
         sec[4][33] != 255) {
        undef_val2 = 0;
        val2 = scaled2flt(INT1(sec[4][29]), int4(sec[4] + 30));
    }

    sprintf(inv_out,"%s %d", name,code_table_4_5a(sec));
    inv_out += strlen(inv_out);

    /* i is number of fields to print */    
    i = 3;
    if (level_type2 == 255) i = 1;
    if (level_type1 == level_type2) i = 2;
    if (i == 2 && undef_val2 == 1) i = 1;
    if (i == 1 && undef_val1 == 1) i = 0;

    if (i >= 1) {
        if (undef_val1) sprintf(inv_out,",");
        else if (val1 == (int) val1)  sprintf(inv_out,",%d", (int) val1);
        else sprintf(inv_out,",%g", val1);
        inv_out += strlen(inv_out);
    }
    if (i >= 2) {
        if (undef_val2) sprintf(inv_out,",");
        else if (val2 == (int) val2)  sprintf(inv_out,",%d", (int) val2);
        else sprintf(inv_out,",%g", val2);
        inv_out += strlen(inv_out);
    }
    if (i >= 3) {
        sprintf(inv_out,",%d", level_type2);
        inv_out += strlen(inv_out);
    }

    sprintf(inv_out," %d,%d,%d", GB2_Discipline(sec), GB2_ParmCat(sec), GB2_ParmNum(sec));
    inv_out += strlen(inv_out);

    /* statistical processing */
    string = NULL;
    i = code_table_4_10(sec);
    if (i != -1 && i != 255) {
        sprintf(inv_out,",%d", i);
        inv_out += strlen(inv_out);
    }

    /* add stat processing term to description */
    if (i != -1) {
        string = NULL;
        switch(i) {
#include "CodeTable_4.10.dat"
        }
    }
    if (string == NULL) sprintf(inv_out," %s %s [%s]", "none", desc, unit);
    else sprintf(inv_out," %s %s [%s]", string, desc, unit);

	
    inv_out += strlen(inv_out);

    return 0;
}

/*
 * HEADER:200:lev0:inv:0:level (for g2ctl)
 */
int f_lev0(ARG0) {
    int prod_def, center, subcenter;
    int level_type1, level_type2;
    float val1, val2;

    if (mode < 0) return 0;
    prod_def = GB2_ProdDefTemplateNo(sec);
    center = GB2_Center(sec);
    subcenter = GB2_Subcenter(sec);

    level_type1 = sec[4][22];
    val1 = scaled2flt(INT1(sec[4][23]), int4(sec[4] + 24));

    level_type2 = sec[4][28];
    val2 = scaled2flt(INT1(sec[4][29]), int4(sec[4] + 30));
    inv_out[0] = 0;

    if (level_type2 == 255) {
        if (center == 7 && level_type1 >= 192 && level_type1 <= 254) {
           switch(level_type1) {
           case 200: strcpy(inv_out,"clm"); return 0;
           case 201: strcpy(inv_out,"ocn"); return 0;
           case 204: strcpy(inv_out,"top0C"); return 0;
           case 206: strcpy(inv_out,"gclb"); return 0;
           case 207: strcpy(inv_out,"gclt"); return 0;
           case 209: strcpy(inv_out,"blclb"); return 0;
           case 210: strcpy(inv_out,"blclt"); return 0;
           case 211: strcpy(inv_out,"blcll"); return 0;
           case 212: strcpy(inv_out,"lclb"); return 0;
           case 213: strcpy(inv_out,"lclt"); return 0;
           case 214: strcpy(inv_out,"lcll"); return 0;
           case 220: strcpy(inv_out,"pbl"); return 0;
           case 222: strcpy(inv_out,"mclb"); return 0;
           case 223: strcpy(inv_out,"mclt"); return 0;
           case 224: strcpy(inv_out,"mcll"); return 0;
           case 232: strcpy(inv_out,"hclb"); return 0;
           case 233: strcpy(inv_out,"hclt"); return 0;
           case 234: strcpy(inv_out,"hcll"); return 0;
           case 235: sprintf(inv_out,"ocn%gC",val1/10); return 0;
           case 237: sprintf(inv_out,"ocnml%gm",val1); return 0;
           case 238: sprintf(inv_out,"ocnil%gm",val1); return 0;
           case 242: strcpy(inv_out,"cclb"); return 0;
           case 243: strcpy(inv_out,"cclt"); return 0;
           case 244: strcpy(inv_out,"ccll"); return 0;
           case 248: strcpy(inv_out,"scclb"); return 0;
           case 249: strcpy(inv_out,"scclt"); return 0;
           case 250: strcpy(inv_out,"dcclb"); return 0;
           case 251: strcpy(inv_out,"scclt"); return 0;
	   }
	}

       switch(level_type1) {
       case 1: strcpy(inv_out,"sfc"); break;
       case 2: strcpy(inv_out,"clb"); break;
       case 3: strcpy(inv_out,"clt"); break;
       case 4: strcpy(inv_out,"0C"); break;
       case 6: strcpy(inv_out,"mwl"); break;
       case 7: strcpy(inv_out,"trop"); break;
       case 8: strcpy(inv_out,"toa"); break;
       case 100: sprintf(inv_out,"%dmb",RINT(val1/100)); break;
       case 101: strcpy(inv_out,"msl"); break;
       case 102: sprintf(inv_out,"_%dm",RINT(val1)); break;
       case 103: sprintf(inv_out,"%dm",RINT(val1)); break;
       case 104: sprintf(inv_out,"sig%d",RINT(1000*val1)); break;
       case 105: sprintf(inv_out,"hlev%d",RINT(val1)); break;
       case 109: if (val1 >= 0) { sprintf(inv_out,"%dpv",RINT(val1*1e6)); }
                 else { sprintf(inv_out,"neg%dpv",RINT(-val1*1e6)); }
                 break;
       case 160: sprintf(inv_out,"bsl%dm",RINT(val1)); break;
       default: sprintf(inv_out,"l%d",level_type1); break;
       }
       return 0;
    }
    if (level_type1 == level_type2) {
        switch(level_type1) {
        case 104: sprintf(inv_out,"sg%d_%d",RINT(1000*val1),RINT(1000*val2)); break;
	case 106: sprintf(inv_out, "%d_%dcm",RINT(100*val1),RINT(100*val2)); break;
	case 108: sprintf(inv_out, "%d_%dmb",RINT(val1/100),RINT(val2/100)); break;
        default: sprintf(inv_out,"l%d_%d",level_type1,level_type2); break;
        }
       return 0;
    }
   return 0;
}

