/* @(#)vpftable.h	1.4 5/18/93 18:03:16 */

/* ================ SCCS_ID[] = "@(#) vpftable.h 2.2 10/29/91"  =================

   Environmental Systems Research Institute (ESRI) Applications Programming

       Project: 		Conversion from ARC/INFO to VPF
       Original Coding:		Barry Michaels	April 1991
       Modifications:		David Flinn	April 1991
                                Barry           June 1991
                                Dave            October 1991
   ======================================================================== */

#ifndef _VPF_TABLE_H_

#define _VPF_TABLE_H_

#include <stdio.h>

#if ( defined(CYGWIN) )
#include </usr/include/mingw/values.h>
#else
#include <values.h>
#endif

#include <math.h>

#include "vpfio.h"


#ifndef TRUE
#define TRUE  1
#endif
#ifndef FALSE
#define FALSE 0
#endif

/* Possible byte ordering of the data */

#define LEAST_SIGNIFICANT 0
#define MOST_SIGNIFICANT  1

#define VPF_DIR_SEPARATOR '\\'

/* 
*  The next definition is machine-specific and may need to be changed
*  before recompiling on a different machine
*/
#if defined( DEC_ALPHA ) || defined( LINUX ) || defined(LINUX32) || defined(LINUX64) || defined(IA64) || defined( SCO ) || defined(IA64) || defined(MACINTOSH) || defined(MACINTEL) || defined(CYGWIN) 
#define MACHINE_BYTE_ORDER LEAST_SIGNIFICANT
#else
#define MACHINE_BYTE_ORDER MOST_SIGNIFICANT
#endif
#define DIR_SEPARATOR '/'

/* This should be the ISO definition of date */

typedef char date_type[21] ;	/* Include null end of string */

/* NULL value type */

typedef union {
   char 	  *Char ;
   short int	  Short ;
   int       Int ;
   float	  Float ;
   double	  Double ;
   date_type	  Date ;
   char		  Other ;
} null_field;

/* The basic information carried for each field in the table */
typedef struct {
 char     *name;            /* Name of the field */
 char      description[81]; /* Field description */
 char      keytype;         /* Type of key - (P)rimary, (F)oreign, (N)onkey */
 char      vdt[13];         /* Value description table name */
 char     *tdx;	 	    /* Thematic index file name */
 char      type;            /* Data type - T,I,F,K,D */
 int  count;           /* Number of items in this column (-1 =>variable) */
 char     *narrative;       /* Name of a narrative table describing the field */
 null_field nullval;        /* This is used for the converter */
} header_cell, *header_type;

typedef enum { ram, disk, either, compute } storage_type;
#define RAM 0
#define DISK 1
#define EITHER 2
#define COMPUTE 3

typedef enum { Read, Write, Update } file_mode ;

#define CLOSED 0
#define OPENED 1

/* 
*  Each column in a table row has a count and a pointer to the data
*  and a null value default
*/
typedef struct {
   int   count;
   void *ptr;
} column_type;

/* A table row is an array of columns */
typedef column_type *row_type;

/* 
*  Index for variable width tables.
*  One index cell for each row in the table. 
*/
typedef struct {
   unsigned int pos;
   unsigned int length;
} index_cell, *index_type;

/* VPF table structure: */
typedef struct {
   char           name[13];        /* Name of the VPF table */
   char           *path;           /* Directory path to the table */
   int       nfields;         /* Number of fields */
   char           description[81]; /* Table description */
   char           narrative[13];   /* Table narrative file name */
   header_type    header;          /* Table header structure */
   FILE           *xfp;            /* Index file pointer */
   index_type     index;           /* Index structure */
   storage_type   xstorage;        /* Flag indicating where index stored */
   FILE           *fp;             /* Table file pointer */
   int       nrows;           /* Number of rows in the table */
   row_type       *row;            /* Array of table rows */
   int       reclen;          /* Table record length (-1 => variable */
   int       ddlen;           /* Data definition string length */
   char           *defstr ;        /* rdf, definition string */
   storage_type   storage;         /* Flag indicating table storage method */
   file_mode	  mode ;	   /* Table is either reading or writing */
   unsigned char  status;          /* VPF table status - OPENED or CLOSED */
   unsigned char  byte_order;      /* Byte order of the table's data */
   int       size;            /* Size of the table in bytes */
} vpf_table_type;

typedef struct {
   float x,y;
} coordinate_type;

typedef struct {
   double x,y;
} double_coordinate_type;

typedef struct {
   float x,y,z;
} tri_coordinate_type;

typedef struct {
   double x,y, z;
} double_tri_coordinate_type;

typedef enum { key_id, key_tile, key_exid } key_field_type;

/* These macros help determine the type in the key datatype */

#define TYPE0(cell) ((cell>>6)&(3))
#define TYPE1(cell) ((cell>>4)&(3))
#define TYPE2(cell) ((cell>>2)&(3))
#define TYPE3(cell) ((cell)&(3))

/* These macros set the value in the key datatype */

#define SETTYPE0(cell,value) cell = (((cell)&(~(3<<6)))|(((3)&(value))<<6))
#define SETTYPE1(cell,value) cell = (((cell)&(~(3<<4)))|(((3)&(value))<<4))
#define SETTYPE2(cell,value) cell = (((cell)&(~(3<<2)))|(((3)&(value))<<2))
#define SETTYPE3(cell,value) cell = (((cell)&(~(3)))|(((3)&(value))))

/* This macro helps to write out a key */

#define ASSIGN_KEY(tYPE,kEY,loc,val) { \
    if (val < 1) { \
      tYPE(kEY.type,0); \
    } else if (val < (1<<8)) { \
      tYPE(kEY.type,1); \
      kEY.loc = val ; \
    } else if ( val < (1<<16)) { \
      tYPE(kEY.type,2); \
      kEY.loc = val; \
    } else { \
      tYPE(kEY.type,3); \
      kEY.loc = val; }}

/* define NULL values */


#if ( defined(CYGWIN) )
#include </usr/include/mingw/values.h>
#else
#include <values.h>
#endif

#include <math.h>

#define		VARIABLE_STRING_NULL_LENGTH	10
#define 	NULLCHAR	' '
#define		NULLTEXT	" "
#define		NULLSHORT	-MAXSHORT
#define         NULLINT         -MAXLONG
                               /*12345678901234567890*/
#define		NULLDATE	"                    "
#if UNIX
#define		NULLFLOAT	((float)  quiet_notan (0))
#define		NULLDOUBLE	((double) quiet_notan (0))
#else
#define		NULLFLOAT	((float) MAXFLOAT)
#define		NULLDOUBLE	((double) MAXFLOAT)
#endif

typedef union {
   unsigned char      f1;
   unsigned short int f2;
   unsigned int       f3;
} key_field;

/* id triplet internal storage type */
typedef struct {
   unsigned char type;
   int  id, tile, exid;
} id_triplet_type;

typedef enum { idle_state, name_state, type_state,
	       tuple_state, count_state } ddef_state_type;

/* Functions: */

char *parse_get_string(int *ind,char *src,char delimit);
char  parse_get_char  (int *ind,char *src);
 int parse_get_number(int *ind,char *src,char delemit);

int parse_data_def( vpf_table_type *table );

char *read_text_defstr( FILE *infile, FILE *outerr );

int index_length( int row_number, vpf_table_type table );

int index_pos( int row_number, vpf_table_type table );

static int row_offset( int field, row_type row, vpf_table_type table);

row_type  read_next_row( vpf_table_type table );

row_type  read_row( int row_number, vpf_table_type table );

vpf_table_type vpf_open_table( char *tablename, storage_type storage,
			       char *mode, char *defstr );   /* rdf added */

row_type get_row( int row_number, vpf_table_type table );

void free_row( row_type row, vpf_table_type table );

int  table_pos( char *field_name, vpf_table_type table );

void *get_table_element( int field_number, row_type row,
			 vpf_table_type table, void *value, int *count );

void *named_table_element( char *field_name, row_type row,
			   vpf_table_type table, void *value, int *count);

void *table_element( int field_number, int row_number,
		     vpf_table_type table, void *value, int *count );
 
int get_table_int( int field_number,
                        row_type row,
                        vpf_table_type table,
                        key_field_type keyfield );
 
double get_table_number( int field_number,
                         row_type row,
                         vpf_table_type table,
                         key_field_type keyfield );
 
void vpf_close_table( vpf_table_type *table );

int  is_vpf_table( char *fname );

/* Write functions */

int write_next_row( row_type row, vpf_table_type *table );

row_type create_row( vpf_table_type table );

void destroy_table_element(int field, row_type row, vpf_table_type table);

void nullify_table_element(int field, row_type row, vpf_table_type table);

int put_table_element( int field, row_type row, vpf_table_type table,
			    void *value, int count );
 
int put_table_int( int field_number,
                        row_type row,
                        vpf_table_type table,
                        int value,
                        key_field_type keyfield );
 
int put_table_number( int field_number,
                           row_type row,
                           vpf_table_type table,
                           double value,
                           key_field_type keyfield );
 
void swap_two ( char *in, char *out );

void swap_four ( char *in, char *out );

void swap_eight ( char *in, char *out );

int calculate_row_size(row_type row, vpf_table_type *table);
 
int calculate_key_size(id_triplet_type key);
 
void adjust_table_size(int rownum, vpf_table_type *table, int offset);
 
void increase_table_size(int rownum, vpf_table_type *table,
                         int offset);
 
void decrease_table_size(int rownum, vpf_table_type *table,
                         int offset);
 
void adjust_index(int startrow, vpf_table_type *table, int offset);


#ifndef Max
#define Max(a,b)     ((a > b) ? a : b)
#endif

#ifndef Min
#define Min(a,b)     ((a < b) ? a : b)
#endif

#endif     /* #ifndef _VPF_TABLE_H_  */
