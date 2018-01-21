#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "grb2.h"
#include "wgrib2.h"
#include "fnlist.h"

/*
 * CodeTable.c
 *
 * returns values of various code tables
 *  use this routine instead of .h files
 *
 * 12/2006 Public Domain Wesley Ebisuzaki
 * 1/2007 M. Schwarb
 * 1/2008 S. Varlamov fixed to code_table_5.5
 */

/*
 * HEADER:-1:code_table_0.0:inv:0:code table 0.0 discipline
 */

int f_code_table_0_0(ARG0) {
    int val;
    char *string;

    if (mode >= 0) {
        val = code_table_0_0(sec);
        string = NULL;
        switch(val) {
#include "CodeTable_0.0.dat"
        }
        if (string) sprintf(inv_out,"code table 0.0=%d %s", val, string);
        else sprintf(inv_out,"code table 0.0=%d", val);
    }
    return 0;
}
int code_table_0_0(unsigned char **sec) {
    return  (int) sec[0][9];
}

/*
 * HEADER:-1:code_table_1.0:inv:0:code table 1.0 master table version
 */

int f_code_table_1_0(ARG0) {

    if (mode >= 0) {
        sprintf(inv_out,"code table 1.0=%d", code_table_1_0(sec));
    }
    return 0;
}
int code_table_1_0(unsigned char **sec) {
    return  (int) sec[1][9];
}

/*
 * HEADER:-1:code_table_1.1:inv:0:code table 1.1 local table version
 */

int f_code_table_1_1(ARG0) {

    if (mode >= 0) {
        sprintf(inv_out,"code table 1.1=%d", code_table_1_1(sec));
    }
    return 0;
}
int code_table_1_1(unsigned char **sec) {
    return  (int) sec[1][10];
}


/*
 * HEADER:-1:code_table_1.2:inv:0:code table 1.2 significance of reference time
 */

int f_code_table_1_2(ARG0) {
    char *string;
    int val;
    if (mode >= 0) {
        val = code_table_1_2(sec);
        string = NULL;
        switch(val) {
#include "CodeTable_1.2.dat"
        }
	if (string == NULL) sprintf(inv_out,"code table 1.2=%d", val);
        else sprintf(inv_out,"code table 1.2=%d %s", val, string);
    }
    return 0;
}
int code_table_1_2(unsigned char **sec) {
    return  (int) sec[1][11];
}

/*
 * HEADER:-1:code_table_1.3:inv:0:code table 1.3 production status of processed data
 */

int f_code_table_1_3(ARG0) {
    char *string;
    int val;

    if (mode >= 0) {
        val = code_table_1_3(sec);
        string = NULL;
	switch(val) {
#include "CodeTable_1.3.dat"
        }
        if (string == NULL) sprintf(inv_out,"code table 1.3=%d", val);
        else sprintf(inv_out,"code table 1.3=%d %s", val, string);
    }
    return 0;
}
int code_table_1_3(unsigned char **sec) {
    return  (int) sec[1][19];
}


/*
 * HEADER:-1:code_table_1.4:inv:0:code table 1.4 type of processed data
 */

int f_code_table_1_4(ARG0) {
    char *string;
    int val;
    if (mode >= 0) {
        val = code_table_1_4(sec);
	string = NULL;
	switch(val) {
#include "CodeTable_1.4.dat"
	}
	if (string == NULL) sprintf(inv_out,"code table 1.4=%d", val);
	else sprintf(inv_out,"code table 1.4=%d %s", val, string);
    }
    return 0;
}
int code_table_1_4(unsigned char **sec) {
    return  (int) sec[1][20];
}

/*
 * HEADER:-1:code_table_3.0:inv:0:code table 3.0 Source of grid definition
 */

int f_code_table_3_0(ARG0) {

    if (mode >= 0) {
        sprintf(inv_out,"code table 3.0=%d", code_table_3_0(sec));
    }
    return 0;
}
int code_table_3_0(unsigned char **sec) {
    return  (int) sec[3][5];
}

/*
 * HEADER:-1:code_table_3.1:inv:0:code table 3.1 Grid definition template number
 */

int f_code_table_3_1(ARG0) {
    char *string;
    int val;
    if (mode >= 0) {
        val = code_table_3_1(sec);
	string = NULL;
	switch(val) {
#include "CodeTable_3.1.dat"
	}
	if (string == NULL) sprintf(inv_out,"code table 3.1=%d", val);
	else sprintf(inv_out,"code table 3.1=%d %s", val, string);
    }
    return 0;
}
int code_table_3_1(unsigned char **sec) {
    return  (int) uint2(sec[3]+12);
}

/*
 * HEADER:-1:code_table_3.2:inv:0:code table 3.2 Size and Shape of Earth
 */

int f_code_table_3_2(ARG0) {
    int val;
    char *string;
    if (mode >= 0) {
	val = code_table_3_2(sec);
        if (val >= 0) {
            string = NULL;
            switch(val) {
#include "CodeTable_3.2.dat"
	    }
	    if (string == NULL) sprintf(inv_out,"code table 3.2=%d", val);
	    else sprintf(inv_out,"code table 3.2=%d %s", val, string);
        }
    }
    return 0;
}
int code_table_3_2(unsigned char **sec) {
    int grid_def;
    grid_def = code_table_3_1(sec);
    if (grid_def < 50 || grid_def == 90 || grid_def == 110 || grid_def == 204 || grid_def == 1000 
               || grid_def == 1100) return (int) sec[3][14];
    return  -1;
}


/*
 * HEADER:-1:code_table_3.6:inv:0:code table 3.6 Spectral data representation type
 */

int f_code_table_3_6(ARG0) {
    int val;
    if (mode >= 0) {
	val = code_table_3_6(sec);
	if (val >= 0) sprintf(inv_out,"code table 3.6=%d", val);
    }
    return 0;
}
int code_table_3_6(unsigned char **sec) {
    int grid_def;
    grid_def = code_table_3_1(sec);
    if (grid_def >= 50 && grid_def <= 53)  return GDS_Harmonic_code_3_6(sec[3]);
    return -1;
}

/*
 * HEADER:-1:code_table_3.7:inv:0:code table 3.7 Spectral data representation mode
 */

int f_code_table_3_7(ARG0) {
    int val;
    if (mode >= 0) {
        val = code_table_3_7(sec);
        if (val >= 0) sprintf(inv_out,"code table 3.7=%d", val);
    }
    return 0;
}
int code_table_3_7(unsigned char **sec) {
    int grid_def;
    grid_def = code_table_3_1(sec);
    if (grid_def >= 50 && grid_def <= 53)  return GDS_Harmonic_code_3_7(sec[3]);
    return -1;
}

/*
 * HEADER:-1:code_table_3.8:inv:0:code table 3.8 Grid point position
 */

int f_code_table_3_8(ARG0) {
    int val;
    char *string;
    if (mode >= 0) {
        val = code_table_3_8(sec);
	if (val >= 0) {
            string = NULL;
            switch(val) {
#include "CodeTable_3.8.dat"
            }
            if (string == NULL) sprintf(inv_out,"code table 3.8=%d", val);
            else sprintf(inv_out,"code table 3.8=%d %s", val, string);
        }
    }
    return 0;
}

int code_table_3_8(unsigned char **sec) {
    int grid_def;
    grid_def = code_table_3_1(sec);
    if (grid_def == 100)  return sec[3][31];
    return -1;
}

/*
 * HEADER:-1:code_table_3.11:inv:0:code table 3.11 regional/global thinned/reduced grid

 */
int f_code_table_3_11(ARG0) {
    int val;
    char *string;
    if (mode >= 0) {
        if ((val = code_table_3_11(sec)) >= 0) {
            string = NULL;
            switch(val) {
		case 0: string = "not used";
			break;
		case 1: string = "thinned global grid";
			break;
		case 2: string = "thinned regional grid";
			break;
            }
            if (string == NULL) sprintf(inv_out,"code table 3.11=%d", val);
            else sprintf(inv_out,"code table 3.11=%d %s", val, string);
        }
    }
    return 0;
}
int code_table_3_11(unsigned char **sec) {
    return sec[3][11];
}

/*
 * HEADER:-1:code_table_3.15:inv:0:code table 3.15 Physical meaning of vertical coordinate
 */
int f_code_table_3_15(ARG0) {
    int val;
    char *string;
    if (mode >= 0) {
        if ((val = code_table_3_15(sec)) >= 0) {
            string = NULL;
            switch(val) {
#include "CodeTable_3.15.dat"
            }
            if (string == NULL) sprintf(inv_out,"code table 3.15=%d", val);
            else sprintf(inv_out,"code table 3.15=%d %s", val, string);
        }
    }
    return 0;
}
int code_table_3_15(unsigned char **sec) {
    int grid_def;
    grid_def = code_table_3_1(sec);
    if (grid_def == 1000)  return sec[3][62];
    if (grid_def == 1200)  return sec[3][38];
    return -1;
}

/*
 * HEADER:-1:code_table_3.21:inv:0:code table 3.21 Vertical Dimension coordinate values defn
 */
int f_code_table_3_21(ARG0) {
    int val;
    char *string;
    if (mode >= 0) {
        val = code_table_3_21(sec);
        if (val >= 0) {
            string = NULL;
            switch(val) {
#include "CodeTable_3.21.dat"
            }
            if (string == NULL) sprintf(inv_out,"code table 3.21=%d", val);
            else sprintf(inv_out,"code table 3.21=%d %s", val, string);
        }
    }
    return 0;
}
int code_table_3_21(unsigned char **sec) {
    int grid_def;
    grid_def = code_table_3_1(sec);
    if (grid_def == 1000)  return sec[3][63];
    if (grid_def == 1200)  return sec[3][39];
    return -1;
}

/*
 * HEADER:-1:code_table_4.0:inv:0:code table 4.0 Product Definition Template Number
 */

int f_code_table_4_0(ARG0) {
    int val;
    char *string;
    if (mode >= 0) {
        val = code_table_4_0(sec);
        string = NULL;
        switch(val) {
#include "CodeTable_4.0.dat"
        }
        if (string == NULL) sprintf(inv_out,"code table 4.0=%d", val);
        else sprintf(inv_out,"code table 4.0=%d %s", val, string);
    }
    return 0;
}
int code_table_4_0(unsigned char **sec) {
    return  GB2_ProdDefTemplateNo(sec);
}

/*
 * HEADER:-1:code_table_4.1:inv:0:code table 4.1
 */

int f_code_table_4_1(ARG0) {
    int p;
    if (mode >= 0) {
	p = code_table_4_1(sec);
	if (p >= 0) sprintf(inv_out,"code table 4.1=%d", p);
    }
    return 0;
}
int code_table_4_1(unsigned char **sec) {
    int p;
    p = GB2_ProdDefTemplateNo(sec);
	if (p <= 12 || p == 20 || p == 30 ||  p == 254 || p == 1000 || p == 1001 || 
		p == 1002 || p == 1100 || p == 1101) return (int) sec[4][9];
    return -1;
}

/*
 * HEADER:-1:code_table_4.2:inv:0:code table 4.2
 */

int f_code_table_4_2(ARG0) {
    int p;
    if (mode >= 0) {
	p = code_table_4_2(sec);
	if (p >= 0) sprintf(inv_out,"code table 4.2=%d", p);
    }
    return 0;
}
int code_table_4_2(unsigned char **sec) {
    int p;
    p = GB2_ProdDefTemplateNo(sec);
    if (p <= 12 || p == 20 || p == 30 ||  p == 1000 || p == 1001 || p == 1002 || p == 1100 || p == 1101) 
        return (int) sec[4][10];
    return -1;
}

/*
 * HEADER:-1:code_table_4.3:inv:0:code table 4.3 Type of Generating Process
 */
int f_code_table_4_3(ARG0) {
    int val;
    char *string;
    if (mode >= 0) {
        val = code_table_4_3(sec);
        string = NULL;
        switch(val) {
#include "CodeTable_4.3.dat"
        }
        if (string == NULL) sprintf(inv_out,"code table 4.3=%d", val);
        else sprintf(inv_out,"code table 4.3=%d %s", val, string);
    }
    return 0;
}
int code_table_4_3(unsigned char **sec) {
    int p;
    p = GB2_ProdDefTemplateNo(sec);
    if (p <= 12 || p == 20 || p == 30 ||  p == 1000 || p == 1001 || p == 1002 || p == 1100 || p == 1101) 
        return (int) sec[4][11];
    return -1;
}

/*
 * HEADER:-1:code_table_4.4:inv:0:code table 4.4
 */
int f_code_table_4_4(ARG0) {
    int p;
    if (mode >= 0) {
	p = code_table_4_4(sec);
	if (p >= 0) sprintf(inv_out,"code table 4.4=%d", p);
    }
    return 0;
}
int code_table_4_4(unsigned char **sec) {
    int p;
    p = GB2_ProdDefTemplateNo(sec);
    if (p <= 12 || p == 1000 || p == 1001 || p == 1002 || p == 1100 || p == 1101) 
		return (int) sec[4][17];
    return -1;
}

/*
 * HEADER:-1:code_table_4.5a:inv:0:code table 4.5 (1st value)
 */
int f_code_table_4_5a(ARG0) {
    int p;
    if (mode >= 0) {
	p = code_table_4_5a(sec);
	if (p >= 0) sprintf(inv_out,"code table 4.5a=%d", p);
    }
    return 0;
}
int code_table_4_5a(unsigned char **sec) {
    int p;
    p = GB2_ProdDefTemplateNo(sec);

    if (p <= 12 || p == 1100 || p == 1101) return (int) sec[4][22];
    return -1;
}
/*
 * HEADER:-1:code_table_4.5b:inv:0:code table 4.5 (2nd value)
 */
int f_code_table_4_5b(ARG0) {
    int p;
    if (mode >= 0) {
	p = code_table_4_5b(sec);
	if (p >= 0) sprintf(inv_out,"code table 4.5b=%d", p);
    }
    return 0;
}
int code_table_4_5b(unsigned char **sec) {
    int p;
    p = GB2_ProdDefTemplateNo(sec);

    if (p <= 12 || p == 1100 || p == 1101) return (int) sec[4][28];
    return -1;
}

/*
 * HEADER:-1:code_table_4.6:inv:0:code table 4.6 ensemble type
 */
int f_code_table_4_6(ARG0) {
    int p;
    if (mode >= 0) {
	p = code_table_4_6(sec);
	if (p >= 0) sprintf(inv_out,"code table 4.6=%d", p);
    }
    return 0;
}
int code_table_4_6(unsigned char **sec) {
    int p;
    p = GB2_ProdDefTemplateNo(sec);
    if (p == 1 || p == 11) return (int) sec[4][34];
    return -1;
}

/*
 * HEADER:-1:code_table_4.7:inv:0:code table 4.7 derived forecast
 */
int f_code_table_4_7(ARG0) {
    int val;
    if (mode >= 0) {
	val = code_table_4_7(sec);
	if (val >= 0) sprintf(inv_out,"code table 4.7=%d", val);
    }
    return 0;
}
int code_table_4_7(unsigned char **sec) {
    int val;
    val = GB2_ProdDefTemplateNo(sec);
    switch (val) {
	case 2:
	case 3:
	case 4:
	case 12:
		return (int) sec[4][34]; break;
    }
    return -1;
}
/*
 * HEADER:-1:code_table_4.9:inv:0:code table 4.9 Probability Type
 */
int f_code_table_4_9(ARG0) {
    int val;
    char *string;
    if (mode >= 0) {
        val = code_table_4_9(sec);
        if (val >= 0) {
            string = NULL;
            switch(val) {
#include "CodeTable_4.9.dat"
            }
            if (string == NULL) sprintf(inv_out,"code table 4.9=%d", val);
            else sprintf(inv_out,"code table 4.9=%d %s", val, string);
        }
    }
    return 0;
}
int code_table_4_9(unsigned char **sec) {
    int val;
    val = GB2_ProdDefTemplateNo(sec);
    if (val == 5) return sec[4][36];
    return -1;
}


/*
 * HEADER:-1:code_table_4.10:inv:0:code table 4.10 statistical processing
 */
int f_code_table_4_10(ARG0) {
    int val;
    char *string;
    if (mode >= 0) {
	val = code_table_4_10(sec);
	if (val >= 0) {
	    string = NULL;
	    switch(val) {
#include "CodeTable_4.10.dat"
	    }
	    if (string == NULL) sprintf(inv_out,"code table 4.10=%d", val);
	    else sprintf(inv_out,"code table 4.10=%d %s", val, string);
	}
    }
    return 0;
}
int code_table_4_10(unsigned char **sec) {
    int val, i; 
    val = GB2_ProdDefTemplateNo(sec);
    switch (val) {
	case 8: i = sec[4][46]; break;
	case 9: i = sec[4][59]; break;
	case 10: i = sec[4][47]; break;
	case 11: i = sec[4][49]; break;
	case 12: i = sec[4][48]; break;
	case 13: i = sec[4][80]; break;
	case 14: i = sec[4][76]; break;
	case 1001: i = sec[4][26]; break;
	case 1002: i = sec[4][24]; break;
	case 1101: i = sec[4][38]; break;
	default: i = -1; break;
    }
    return i;
}

/*
 * HEADER:-1:code_table_5.0:inv:0:code table 5.0 data representation number
 */
int f_code_table_5_0(ARG0) {
    int p;
    char *string;
    if (mode >= 0) {
	p = code_table_5_0(sec);
	if (p >= 0) {
	    string = NULL;
	    switch(p) {
#include "CodeTable_5.0.dat"
            }
	    if (string == NULL) sprintf(inv_out,"code table 5.0=%d", p);
	    else sprintf(inv_out,"code table 5.0=%d %s", p, string);
	}
    }
    return 0;
}

int code_table_5_0(unsigned char **sec) {
    return (int) uint2(sec[5]+9);
}

/*
 * HEADER:-1:code_table_5.1:inv:0:code table 5.1 type of original field values
 */
int f_code_table_5_1(ARG0) {
    int p;
    char *string;

    if (mode >= 0) {
        p = code_table_5_1(sec);
        if (p >= 0) {
	    string = NULL;
	    switch(p) {
#include "CodeTable_5.1.dat"
            }
	    if (string == NULL) sprintf(inv_out,"code table 5.1=%d", p);
	    else sprintf(inv_out,"code table 5.1=%d %s", p, string);
	}
    }
    return 0;
}

int code_table_5_1(unsigned char **sec) {

    switch(code_table_5_0(sec)) {
    case 0:
    case 1:
    case 2:
    case 3:
        return (int) (sec[5][20]);
        break;
    default:
	return -1;
	break;
    }
    return -1;
}


/*
 * HEADER:-1:code_table_5.5:inv:0:code table 5.5 missing value management for complex packing
 */
int f_code_table_5_5(ARG0) {
    int p;
    if (mode >= 0) {
        p = code_table_5_5(sec);
        if (p >= 0) sprintf(inv_out,"code table 5.5=%d", p);
    }
    return 0;
}

int code_table_5_5(unsigned char **sec) {
    if (code_table_5_0(sec) < 2 || code_table_5_0(sec) > 3) return -1;
    return (int) (sec[5][22]);
}

/*
 * HEADER:-1:code_table_6.0:inv:0:code table 6.0 Bitmap indicator
 */
int f_code_table_6_0(ARG0) {
    int p;
    if (mode >= 0) {
        p = code_table_6_0(sec);
        if (p >= 0) sprintf(inv_out,"code table 6.0=%d", p);
    }
    return 0;
}

int code_table_6_0(unsigned char **sec) {
    return (int) (sec[6][5]);
}

