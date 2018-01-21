/*
*******************************************************************************
* hico_hs_stats.c
*
* Read a HICO health and status file and output stats in x-minute intervals.
*
* Call:
*   hico_hs_stats <input_file> <interval_minutes> <output_file>
*
*   <input_file> is the input health and status csv filename
*   <interval_minutes> is the minute interval to summarize over
*   <output_file> is the output filename
*
*******************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <ctype.h>

void usage_call(void);
int split_str(char *first_str, char *second_str, char *str_split);

int main (int argc, char *argv[])
{
  FILE *infile, *outfile
  char *in_filename, *out_filename;
  char tmp_str[100], this_line[5000];
  float interval_min;
  int num_fields=238;
  int i, j, year, month, day, hour, minute;
  float fmin[num_fields], fmean[num_fields], fmax[num_fields];
  float fstdev[num_fields], fcount[num_fields];

  /* only run if there are three arguments */
  if (argc == 4)
    {
    in_filename = argv[1];
    interval_min = atof(argv[2]);
    out_filename = argv[3];
    }
  else
    usage_call();

  /* open files */
  infile = fopen(in_filename, "rt");
  if (infile == NULL)
    {
    fprintf(stderr, "Input file %s could not be opened.\n",in_filename);
    exit(-1);
    }
  /* read header line */
  this_line[0] = '\0';
  fgets(this_line,sizeof(this_line),infile);

  /* read in existing output file category breakpoints and data */
  outfile = fopen(out_filename, "wt");
  if (outfile == NULL)
    {
    fprintf(stderr, "Output file %s could not be opened.\n",out_filename);
    exit(-1);
    }

  /* echo header line to output file with extra stat columns */
  for (i=0; i<num_fields-1; i++)
    {
    split_str(&tmp_str[0], &this_line[0], ",");
    fprintf(outfile,"MIN_%s,MAX_%s,MEAN_%s,STDEV_%s,COUNT_%s,",tmp_str,tmp_str,tmp_str,tmp_str,tmp_str);
    }
  fprintf(outfile,"MIN_%s,MAX_%s,MEAN_%s,STDEV_%s,COUNT_%s\n",this_line,this_line,this_line,this_line,this_line);

  /* process the rest of the file */
  while (!feof(infile))
    {
    /* the following is a fail-safe in case the last line of the file is blank */
    this_line[0] = '\0';
    fgets(this_line,sizeof(this_line),infile);
    if (strlen(this_line) > 2) /* not EOF */
      {
      
      }
    }

  /* close files and exit */
  fclose(infile);
  fclose(outfile);

  fprintf(stdout,"Output successfully written to %s\n",out_filename);

  return(0);
  }

/* =============================================================
   split_str - splits a string into two parts
     - second_str should start off containing the string to be
       split
     - when done:
       - first_str will contain everything before the first
         character of str_split
       - second_str will contain the rest of the string,
         starting with the first character after str_split
   ============================================================= */
int split_str(char *first_str, char *second_str, char *str_split)
  {
  char *line_ptr, *tmp_str;
  int i;

  tmp_str = second_str;
  line_ptr = strstr(second_str,str_split);
  for(i=0; i<(strlen(second_str)-strlen(line_ptr)); i++)
    first_str[i] = second_str[i];
  first_str[strlen(second_str)-strlen(line_ptr)] = '\0';
  for(i=strlen(str_split); i<=strlen(line_ptr); i++)
    tmp_str[i-strlen(str_split)] = *(line_ptr+i);
  strcpy(second_str, tmp_str);

  return(0);
  }

/* ============================================================
        usage_call()

   ============================================================ */
void usage_call(void)
  {
  printf("\n");
  printf("hico_hs_stats <input_file> <interval_minutes> <output_file>\n");
  printf("\n");
  printf("  <input_file> = HREP health and status CSV file\n");
  printf("  <interval_minutes> = time interval in minutes over which to calculate stats\n");
  printf("  <output_file> = output filename\n");
  printf("\n");
  exit(0);
  }
