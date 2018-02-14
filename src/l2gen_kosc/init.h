#include "l12_proto.h"

typedef struct loadl1_struct {

    double      radeg;
    int32_t     sensorID;

    int32_t     *Lambda;
    float       *Fobar;
    float       *Tau_r;
    float       *k_oz;
    float       *aw;
    float       *bbw;

} loadl1str;

typedef struct get_f0_thuillier_ext_struct {

    float f0_table[2198];
   
    int32_t firstCall; 

} f0str;



typedef struct init_struct {

    loadl1str   *loadl1rec;
    f0str       *f0rec;

} initstr;




