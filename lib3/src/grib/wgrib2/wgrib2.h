/*
 * public domain 12/2006 wesley ebisuzaki
 *               1/2007 M. Schwarb
 */

#include <stdio.h>
#define VERSION "v0.1.6 1/2008 wesley ebisuzaki, Jaakko Hyvätti, Karl Pfeiffer, Manfred Schwarb, Kristian Nilssen, Sergey Varlamov"
/*  1/2007 M. Schwarb unsigned int ndata */

/* max number of -match options used */
#define MATCH_MAX 20

#define UNDEFINED       9.999e20
#define UNDEFINED_LOW   9.9989e20
#define UNDEFINED_HIGH  9.9991e20
#define UNDEFINED_VAL(x) ((x) >= UNDEFINED_LOW && (x) <= UNDEFINED_HIGH)


/* formatting length for function names for help screen */
#define HELP_NAME_LEN	15
#define N_ARGLIST	200
struct ARGLIST {int fn; int i_argc;};
#define STRING_SIZE	200

/* calling arguements for function API */

#define ARG0	int mode, unsigned char **sec, float *data, unsigned int ndata, char *inv_out, void **local
#define ARG1	int mode, unsigned char **sec, float *data, unsigned int ndata, char *inv_out, void **local, char *arg1
#define ARG2	int mode, unsigned char **sec, float *data, unsigned int ndata, char *inv_out, void **local, char *arg1, char *arg2
#define ARG3	int mode, unsigned char **sec, float *data, unsigned int ndata, char *inv_out, void **local, char *arg1, char *arg2, char *arg3
#define ARG4	int mode, unsigned char **sec, float *data, unsigned int ndata, char *inv_out, void **local, char *arg1, char *arg2, char *arg3, char *arg4

#define CALL_ARG0	mode, sec, data, ndata, inv_out, local

enum input_type {inv_mode, dump_mode, all_mode};
enum output_order_type {raw,wesn,wens};


#define INV_BUFFER	10000


struct gribtab_s {
  int disc;   /* Section 0 Discipline                                */
  int mtab;   /* Section 1 Master Tables Version Number              */
  int cntr;   /* Section 1 originating centre, used for local tables */
  int ltab;   /* Section 1 Local Tables Version Number               */
  int pcat;   /* Section 4 Template 4.0 Parameter category           */
  int pnum;   /* Section 4 Template 4.0 Parameter number             */
  const char *name;
  const char *desc;
  const char *unit;
};
extern struct gribtab_s gribtab[];

double Int_Power(double x, int y);

int int4(unsigned char *);
int int2(unsigned char *);
unsigned int uint4(unsigned char *);
int uint4_missing(unsigned char *);
unsigned int uint2(unsigned char *);
unsigned long int uint8(unsigned char *);
float scaled2flt(int scale_factor, int scale_value);
void uint8_char(unsigned long int i, unsigned char *p);
void uint_char(unsigned int i, unsigned char *p);
void int_char(int i, unsigned char *p);
void uint2_char(unsigned int i, unsigned char *p);
void int2_char(int i, unsigned char *p);


float ieee2flt(unsigned char *ieee);
 
unsigned char *seek_grib2(FILE *file, long int *pos, long int *len_grib,
        unsigned char *buffer, unsigned int buf_len, long int *n_bytes);

int read_grib2(FILE *file, long pos, long len_grib, unsigned char *buffer);

unsigned char *rd_grib2_msg(FILE *input, long int *pos, long int *len, int *num_submsgs);

int parse_1st_msg(unsigned char **sec);
int parse_next_msg(unsigned char **sec);
 
int missing_points(unsigned char *bitmap, int n);

unsigned int nval(unsigned char **sec);

void BDS_unpack(float *flt, unsigned char *bits, unsigned char *bitmap,
        int n_bits, int n, double ref, double scale);

int getName(unsigned char **sec, int mode, char *inv_out, char *name, char *desc, char *unit);

int rd_inventory(int *rec_num, int *submsg, long int *pos);
int get_nxny(unsigned char **sec, int *nx, int *ny, int *npnts, int *res, int *scan);

int wrtieee(float *array, int n, int header, FILE *output);
int flt2ieee(float x, unsigned char *ieee);
int add_time(int *year, int *month, int *day, int *hour, int *minute, int *second, unsigned int dtime, int unit);
int verftime(unsigned char **sec, int *year, int *month, int *day, int *hour, int *minute, int *second);

int is_match(char *string);
char *nc_strstr(char *s, char *t);

int code_table_0_0(unsigned char **sec);

int code_table_1_0(unsigned char **sec);
int code_table_1_1(unsigned char **sec);
int code_table_1_2(unsigned char **sec);
int code_table_1_3(unsigned char **sec);
int code_table_1_4(unsigned char **sec);

int code_table_3_0(unsigned char **sec);
int code_table_3_1(unsigned char **sec);
int code_table_3_2(unsigned char **sec);
int code_table_3_6(unsigned char **sec);
int code_table_3_7(unsigned char **sec);
int code_table_3_8(unsigned char **sec);
int code_table_3_11(unsigned char **sec);
int code_table_3_15(unsigned char **sec);
int code_table_3_21(unsigned char **sec);

int code_table_4_0(unsigned char **sec);
int code_table_4_1(unsigned char **sec);
int code_table_4_2(unsigned char **sec);
int code_table_4_3(unsigned char **sec);
int code_table_4_4(unsigned char **sec);
int code_table_4_5a(unsigned char **sec);
int code_table_4_5b(unsigned char **sec);
int code_table_4_6(unsigned char **sec);
int code_table_4_7(unsigned char **sec);
int code_table_4_9(unsigned char **sec);
int code_table_4_10(unsigned char **sec);
int code_table_5_0(unsigned char **sec);
int code_table_5_1(unsigned char **sec);
int code_table_5_5(unsigned char **sec);
int code_table_6_0(unsigned char **sec);

int flag_table_3_3(unsigned char **sec);
int flag_table_3_4(unsigned char **sec);
int set_flag_table_3_4(unsigned char **sec, unsigned int flag);
int flag_table_3_5(unsigned char **sec);
int flag_table_3_9(unsigned char **sec);
int flag_table_3_10(unsigned char **sec);

unsigned int pds_fcst_time(unsigned char **sec);
int ij2p(int i, int j, int scan_mode);
int to_we_ns_scan(float *data);
int to_we_sn_scan(float *data);
int get_latlon(unsigned char **sec);
void fatal_error(const char *fmt, const char *string);
void fatal_error_i(const char *fmt, const int i);
void set_mode(int new_mode);
int latlon_0(unsigned char **sec);
int new_gds(unsigned char **sec);
int closest_init(unsigned char **sec);
int closest( unsigned char **sec, float plat, float plon);
int regular2ll(unsigned char **sec, float **lat, float **lon);
int polar2ll(unsigned char **sec, float **lat, float **lon);
int gauss2ll(unsigned char **sec, float **lat, float **lon);
int lambert2ll(unsigned char **sec, float **lat, float **lon);
int mercator2ll(unsigned char **sec, float **lat, float **lon);

void flist2bitstream(float *list, unsigned char *bitstream, int ndata, int nbits);


void netcdf_command(int status);


double radius_earth(unsigned char **sec);

extern char *nl;
