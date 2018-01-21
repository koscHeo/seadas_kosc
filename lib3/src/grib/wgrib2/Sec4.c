#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "grb2.h"
#include "wgrib2.h"
#include "fnlist.h"

/*
 * some routines that involve Sec4
 *
 * Public Domain 2006: Wesley Ebisuzaki
 * 1/2007  cleanup M Schwarb
 */

/*
 * HEADER:400:Sec4:inv:0:Sec 4 values (Product definition section)
 */
int f_Sec4(ARG0) {
    if (mode >= 0) {
	sprintf(inv_out,"Sec4 len=%d #vert coordinate=%d Product Defn Template=4.%d", 
          uint4(&(sec[4][0])), uint2(&(sec[4][5])),GB2_ProdDefTemplateNo(sec) );
    }
    return 0;
}

/*
 * HEADER:400:processid:inv:0:process id (locally defined)
 */
int f_processid(ARG0) {
    if (mode >= 0) {
	sprintf(inv_out,"process_id backgnd=%d anl/fcst=%d", sec[4][12], sec[4][13]);
    }
    return 0;
}



/*
 * HEADER:400:0xSec:inv:1:Hex dump of section X (0..8)
 */
int f_0xSec(ARG1) {
    int i;
    unsigned int j, len;
    unsigned char *s;

    if (mode >= 0) {
        i = atoi(arg1);
        if (i < 0  || i > 8) return 1;
        s = sec[i];
        if (s == NULL) {
            sprintf(inv_out,"sec%d is missing", i);
            return 0;
        }

        if (i == 0) { 
	    len = GB2_Sec0_size;
	}
	else if (i == 8){
	    len = GB2_Sec8_size;
        }
        else { 
            len = uint4(&(sec[i][0]));
        }

        if (mode == 0) sprintf(inv_out,"Sec%d(1..%d)=0x",i,len);
	else sprintf(inv_out,"Sec%d(1..%d)=",i,len);
	inv_out += strlen(inv_out);

	if (mode == 0) {
            while (len--) {
         	sprintf(inv_out,"%.2x", *s++);
		inv_out += strlen(inv_out);
            }
	}
	else if (mode == 1) {
            for (j = 1; j <= len; j++) {
		sprintf(inv_out," %.2x", *s++);
		inv_out += strlen(inv_out);
	    }
	}
	else if (mode == 2) {
            for (j = 1; j <= len; j++) {
		sprintf(inv_out,"%d:%.2x", j, *s++);
		inv_out += strlen(inv_out);
		sprintf(inv_out,j % 15 == 0 ? "\n  " : " ");
		inv_out += strlen(inv_out);
	    }
	}
    }
    return 0;
}

/*
 * HEADER:400:var:inv:0:short variable name
 */

int f_var(ARG0) {
    if (mode >= 0) {
        getName(sec, mode, inv_out, NULL, NULL, NULL);
	inv_out += strlen(inv_out);
    }
    return 0;
}

/*
 * HEADER:400:varX:inv:0:raw variable name - discipline mastertab localtab center parmcat parmnum
 */

int f_varX(ARG0) {

    int discipline, center,mastertab,localtab,parmcat,parmnum;

    if (mode >= 0) {
        discipline = GB2_Discipline(sec);
        center = GB2_Center(sec);
        mastertab = GB2_MasterTable(sec);
        localtab = GB2_LocalTable(sec);
        parmcat = GB2_ParmCat(sec);
        parmnum = GB2_ParmNum(sec);
	if (mode == 0) {
            if (parmnum == 255) {
                sprintf(inv_out,"missing");
                return 0;
	    }
	    sprintf(inv_out,"var%d_%d_%d_%d_%d_%d", discipline, mastertab, localtab, center, parmcat, parmnum);
        }
	else if (mode > 0) {
	    if (parmnum == 255) {
                sprintf(inv_out,"missing definition [?]");
                return 0;
            }
            if (parmnum < 192) {
                sprintf(inv_out,"var discipline=%d master_table=%d parmcat=%d parm=%d", 
                      discipline, mastertab, parmcat, parmnum);
            }
            else {
                sprintf(inv_out,"var discipline=%d local_table=%d center=%d parmcat=%d parm=%d",
                discipline, localtab, center, parmcat, parmnum);
            }
        }
    }
    return 0;
}

/*
 * HEADER:440:ftime:inv:0:forecast time
 */

/*
Code Table 4.4: Indicator of unit of time range
Code figure   Meaning
  0	Minute
  1	Hour
  2	Day
  3	Month
  4	Year
  5	Decade (10 years)
  6	Normal (30 years)
  7	Century (100 years)
  8-9	Reserved
 10	3 hours
 11	6 hours
 12	12 hours
 13	Second
 14-191	Reserved
192-254	Reserved for local use
255	Missing
*/


static void print_ftime (int unit1, int value1, int unit2, int value2, int format, char *inv_out) {
    char *units1, *units2;
    units1 = units2 = NULL;


    switch (unit1) {
    case 0:
      units1 = "min";
      break;
    case 1:
      units1 = "hour";
      break;
    case 2:
      units1 = "day";
      break;
    case 3:
      units1 = "month";
      break;
    case 4:
      units1 = "year";
      break;
    case 5:
      value1 *= 10;
      units1 = "year";
      unit1 = 4;
      break;
    case 6:
      value1 *= 30;
      units1 = "year";
      unit1 = 4;
      break;
    case 7:
      value1 *= 100;
      units1 = "year";
      unit1 = 4;
      break;
    case 10:
      value1 *= 3;
      units1 = "hour";
      unit1 = 1;
      break;
    case 11:
      value1 *= 6;
      units1 = "hour";
      unit1 = 1;
      break;
    case 12:
      value1 *= 12;
      units1 = "hour";
      unit1 = 1;
      break;
    case 13:
      units1 = "sec";
      break;
    case 255:
      units1 = "missing";
      break;
    default:
      units1 = "?";
      break;
    }

    if (format != 1) {
        switch (unit2) {
        case 0:
          units2 = "min";
          break;
        case 1:
          units2 = "hour";
          break;
        case 2:
          units2 = "day";
          break;
        case 3:
          units2 = "month";
          break;
        case 4:
          units2 = "year";
          break;
        case 5:
          value2 *= 10;
          units2 = "year";
          unit2 = 4;
          break;
        case 6:
          value2 *= 30;
          units2 = "year";
          unit2 = 4;
          break;
        case 7:
          value2 *= 100;
          units2 = "year";
          unit2 = 4;
          break;
        case 10:
          value2 *= 3;
          units2 = "hour";
          unit2 = 1;
          break;
        case 11:
          value2 *= 6;
          units2 = "hour";
          unit2 = 1;
          break;
        case 12:
          value2 *= 12;
          units2 = "hour";
          unit2 = 1;
          break;
        case 13:
          units2 = "sec";
          break;
        case 255:
          units1 = "missing";
          break;
        default:
          units2 = "?";
          break;
        }
    }

    if (format == 1) sprintf(inv_out,"%d %s",value1,units1);
    else if (format == 2) {
	if (unit1 == unit2) sprintf(inv_out,"%d-%d %s",value1,value1+value2,units1);
	else {
	    sprintf(inv_out,"%d %s-(%d %s+%d %s",value1,units1,value1,units1,value2,units2);
	}
    }
}


int f_ftime(ARG0) {
    int unit, unit2, value,value2,n,i,consistent,code_4_11;
    char *string;
    if (mode >= 0) {
	if (GB2_ProdDefTemplateNo(sec) != 8) {
	    if ((unit = code_table_4_4(sec)) < 0) return -1;
                value = GB2_ForecastTime(sec);
                if (value == 0) {
                    sprintf(inv_out,"anl");
                }
                else {
    	            print_ftime(unit, value, 0, 0, 1, inv_out);
		    inv_out += strlen(inv_out);
                    sprintf(inv_out," fcst");
                }
            return 0;
        }
        else {
	    if ((unit = code_table_4_4(sec)) < 0) return -1;
            n = (int) sec[4][41];

	    /* simple case n = 1 */
            if (n == 1) {
                value = GB2_ForecastTime(sec);
	        unit2 = sec[4][48];
                value2 = uint4(sec[4]+49);
	        print_ftime(unit,value,unit2,value2,2, inv_out);
	        inv_out += strlen(inv_out);

		/* average or accumulation or */
                string = "";
	        switch(sec[4][46]) {
#include           "CodeTable_4.10.dat"
                }
                sprintf(inv_out," %s",string);
	        inv_out += strlen(inv_out);

		/* analysis or forecast or */
                string = "";
	        switch(code_table_4_3(sec)) {
#include           "CodeTable_4.3.dat"
                }
                sprintf(inv_out," %s",string);
	        inv_out += strlen(inv_out);
                return 0;
            }
		/* check for consistent specs */
		consistent = 1;
		for (i = 12; i < n*12; i++) {
		    if (sec[4][46+(i % 12)] != sec[4][46+i]) {
			consistent = 0;
			break;
		    }
		}


		/* type of increment */
		code_4_11 = sec[4][47];

		/* climo type calculation */
		if (consistent && code_4_11 == 1) {
                    value = GB2_ForecastTime(sec);
                    unit2 = sec[4][48];
                    value2 = uint4(sec[4]+49);
                    print_ftime(unit,value,unit2,value2,2, inv_out);
                    inv_out += strlen(inv_out);

                    /* average or accumulation or */
                    string = "";
                    switch(sec[4][46]) {
#include               "CodeTable_4.10.dat"
                    }
                    sprintf(inv_out," %s",string);
                    inv_out += strlen(inv_out);

                    /* analysis or forecast or */
                    string = "";
                    switch(code_table_4_3(sec)) {
#include               "CodeTable_4.3.dat"
                    }
                    sprintf(inv_out," %s",string);
                    inv_out += strlen(inv_out);

		    /* time between intervals */
		    value = int4(sec[4]+54);
		    unit = sec[4][53];
	            sprintf(inv_out,"@");
	            inv_out += strlen(inv_out);
    	            print_ftime(unit, value, 0, 0, 1, inv_out);
	            inv_out += strlen(inv_out);
                    sprintf(inv_out,"*%d sample", n);
	            inv_out += strlen(inv_out);
                    return 0;
		}

                sprintf(inv_out," help -- need more code ftime");
	        inv_out += strlen(inv_out);

        }
    }
    return 0;
}
