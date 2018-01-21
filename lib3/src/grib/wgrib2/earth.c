#include <stdio.h>
#include "grb2.h"
#include "wgrib2.h"
#include "fnlist.h"

/* function to return the size of the earth */

double radius_earth(unsigned char **sec) {
    int table_3_2;
    double radius;

    table_3_2 = code_table_3_2(sec);

    /* set a default value .. not sure what to do with most values */
    radius = 6367.47;
    if (table_3_2 == 2)  radius = (6378.160 + 6356.775)*0.5;
    if (table_3_2 == 4)  radius = (6378.137 + 6356.772)*0.5;
    if (table_3_2 == 6)  radius = 6371.2290;

    return radius * 1000.0;
}

