#ifndef MAPLISTS_H /* avoid re-inclusion */
#define MAPLISTS_H

/*#define L3M_PARAMS 70*/
int32 L3M_PARAMS = 0;

char **parmname_list;
char **parmname_short;
char **unit_list;
char **scaling_list;
float32 *maximum_list;
float32 *minimum_list;
char **palette_list;
char **precision_list;

float32 base = 10.0;

const char *measure_list[] = {
        "Mean",
        "Variance",
        "Standard deviation",
        "Pixels per bin",
        "Scenes per bin"
};

#endif /* MAPLISTS_H */


