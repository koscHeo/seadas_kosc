#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <netcdf.h>
#include <mfhdf.h>

#include "l2bin_input.h"
#include "genutils.h"
#include "passthebuck.h"
#include "sensorInfo.h"


int get_item(char *arg, instr *input);
int par_file(char *pfile, instr *input);



int input_init(instr *input_str)
{
  input_str->infile[0] = '\0';
  input_str->ofile[0] = '\0';
  input_str->pfile[0] = '\0';

  input_str->fileuse[0] = '\0';
  input_str->qual_prod[0] = '\0';

  input_str->oformat[0] = '\0';

  strcpy(input_str->pversion, "Unspecified");
  strcpy(input_str->prodtype, "day");
  strcpy(input_str->average, "standard");

  strcpy(input_str->l3bprod, ":ALL:");

  /*  strcpy(input_str->flaguse, DEF_FLAG);*/

  input_str->sday = 1970001;
  input_str->eday = 2038018;

  input_str->resolve[0] = '\0';

  /*  input_str->rowgroup = 2160/12;*/
  input_str->rowgroup = -1;
  input_str->interp = -1;
  input_str->noext = -1;

  input_str->night = 0;
  input_str->verbose = 0;
  input_str->minobs = 0;
  input_str->healpix = 0;

  input_str->meminfo = 0;
  input_str->dcinfo = 0;

  input_str->latsouth = -90.0;
  input_str->latnorth = +90.0;
  input_str->lonwest = 0.0;
  input_str->loneast = 0.0;

  input_str->qual_max = 255;

  input_str->deflate = 0;

  strcpy( input_str->suite, "");
  return 0;
}



/*-----------------------------------------------------------------------------
    Function:  l2bin_input

    Returns:   int (status)
        The return code is a negative value if any error occurs, otherwise,
        returns 0.

    Description:
        Convert the arguments from the command line into a structure input
        variable.

    Parameters: (in calling order)
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int       argc          I   number of arguments
        char      **argv        I   list of arguments
        instr     input         O   structure variable for inputs

----------------------------------------------------------------------------*/

int l2bin_input(int argc, char **argv, instr *input)
{
  int  i, status;
  char str_buf[4096];
  char small_buf[256];
  int32_t sd_id;

  char *tmp_str;
  char *char_ptr;
  int sensorID;
  char instrument[256];
  char platform[256];

  FILE *fp=NULL;

  if (argc == 1) return(-1);

  /*                                                                  */
  /* Set input values to defaults                                     */
  /*                                                                  */
  if ( input_init( input ) != 0 ) {
      printf("-E- %s: Error initializing input structure.\n",__FILE__);
      return(-1);
  }

  /*                                                                  */
  /* Loop through command arguments and update input struct           */
  /*                                                                  */
  for (i=1; i<argc; i++)
    if (get_item(argv[i], input) != 0)
      return -1;


  /* Single HDF input */
  /* ---------------- */
  if (Hishdf(input->infile) == TRUE || nc_open(input->infile, 0, &sd_id) == 0) {
    if (Hishdf(input->infile) == TRUE) {
      sd_id = SDstart(input->infile, DFACC_RDONLY);
      SDreadattr(sd_id, SDfindattr(sd_id, "Sensor Name"), small_buf);
      sensorID = sensorName2SensorId(small_buf);
      SDend(sd_id);
    } else {
        DPTB(nc_get_att( sd_id, NC_GLOBAL, "instrument", instrument));
        DPTB(nc_get_att( sd_id, NC_GLOBAL, "platform", platform));
        sensorID = instrumentPlatform2SensorID(instrument,platform);
      nc_close( sd_id);
    }
  }
  else {
    /* Open L2 input files */
    /* ------------------- */
    fp = fopen(input->infile, "r");
    if (fp == NULL) {
      printf("Input listing file: \"%s\" not found.\n", input->infile);
      return -1;
    }
    fgets(str_buf, 256, fp);
    str_buf[strlen(str_buf)-1] = 0;

    if (Hishdf(str_buf) == TRUE) {
      sd_id = SDstart(str_buf, DFACC_RDONLY);
      SDreadattr(sd_id, SDfindattr(sd_id, "Sensor Name"), small_buf);
      SDend(sd_id);
      sensorID = sensorName2SensorId(small_buf);
    } else {
      nc_close( sd_id);
      status = nc_open(str_buf, NC_NOWRITE, &sd_id);
      if (status == 0) {

          DPTB(nc_get_att( sd_id, NC_GLOBAL, "instrument", instrument));
          DPTB(nc_get_att( sd_id, NC_GLOBAL, "platform", platform));
          sensorID = instrumentPlatform2SensorID(instrument,platform);

	nc_close( sd_id);
      } else {
          // Assume HDF5 VIIRS file (JMG  10/18/13)
          strcpy(small_buf, "viirsn");
          sensorID = sensorName2SensorId(small_buf);
      }
    }
    fclose(fp);
  }

  if(sensorID == -1) {
      printf("-E- Can not look up sensor ID for %s.\n", small_buf);
      return(1);
  }

  if ((tmp_str = getenv("OCDATAROOT")) == NULL) {
    printf("OCDATAROOT environment variable is not defined.\n");
    return(1);
  }

  strcpy(str_buf,tmp_str);
  strcat(str_buf,"/");
  strcat(str_buf,sensorDir[sensorID]);
  strcat(str_buf,"/");
  strcat(str_buf,"l2bin_defaults.par");
  par_file(str_buf, input);

  /*                                                                  */
  /* Loop through command arguments and update input struct           */
  /*                                                                  */
  for (i=1; i<argc; i++)
    if (get_item(argv[i], input) != 0)
      return -1;

  // Check for suite entry
  if ( input->suite[0] != 0) {
    strcpy(str_buf,tmp_str);

    strcat(str_buf,"/");
    strcat(str_buf,sensorDir[sensorID]);
    strcat(str_buf,"/");
    strcat(str_buf,"l2bin_defaults_");
    strcat(str_buf, input->suite);
    strcat(str_buf, ".par");
    par_file(str_buf, input);

    for (i=1; i<argc; i++)
      if (get_item(argv[i], input) != 0)
	return -1;
  }

  /*                                                                  */
  /* Build string of parameters for metadata                          */
  /*                                                                  */
  sprintf(str_buf, "infile = %s\n",input->infile);
  strcat(input->parms, str_buf);
  sprintf(str_buf, "ofile = %s\n",input->ofile);
  strcat(input->parms, str_buf);
  sprintf(str_buf, "oformat = %s\n",input->oformat);
  strcat(input->parms, str_buf);
  sprintf(str_buf, "fileuse = %s\n",input->fileuse);
  strcat(input->parms, str_buf);
  /*
  sprintf(str_buf, "PFILE = %s\n",input->pfile);
  strcat(input->parms, str_buf);
  */
  sprintf(str_buf, "sday = %d\n",input->sday);
  strcat(input->parms, str_buf);
  sprintf(str_buf, "eday = %d\n",input->eday);
  strcat(input->parms, str_buf);

  sprintf(str_buf, "resolve = %s\n",input->resolve);
  strcat(input->parms, str_buf);

  sprintf(str_buf, "rowgroup = %d\n",input->rowgroup);
  strcat(input->parms, str_buf);

  sprintf(str_buf, "flaguse = %s\n",input->flaguse);
  strcat(input->parms, str_buf);

  sprintf(str_buf, "l3bprod = %s",&input->l3bprod[1]);
  str_buf[strlen(str_buf)-1] = 0;
  strcat(input->parms, str_buf);
  strcat(input->parms, "\n");

  sprintf(str_buf, "prodtype = %s\n",input->prodtype);
  strcat(input->parms, str_buf);

  sprintf(str_buf, "pversion = %s\n",input->pversion);
  strcat(input->parms, str_buf);

  sprintf(str_buf, "suite = %s\n",input->suite);
  strcat(input->parms, str_buf);

  sprintf(str_buf, "average = %s\n",input->average);
  strcat(input->parms, str_buf);

  /*
    Disable interpolation
  sprintf(str_buf, "INTERP = %d\n",input->interp);
  strcat(input->parms, str_buf);
  */

  sprintf(str_buf, "night = %d\n",input->night);
  strcat(input->parms, str_buf);

  sprintf(str_buf, "verbose = %d\n",input->verbose);
  strcat(input->parms, str_buf);

  sprintf(str_buf, "minobs = %d\n",input->minobs);
  strcat(input->parms, str_buf);

  sprintf(str_buf, "deflate = %d\n",input->deflate);
  strcat(input->parms, str_buf);

  sprintf(str_buf, "qual_prod = %s\n",input->qual_prod);
  strcat(input->parms, str_buf);

  sprintf(str_buf, "qual_max = %d\n",input->qual_max);
  strcat(input->parms, str_buf);

  sprintf(str_buf, "healpix = %d\n",input->healpix);
  strcat(input->parms, str_buf);
  return 0;
}

int get_item(char *arg, instr *input)
{
  char *tmp_str; 
  char keyword [20]; 
  char parm_str[4096];
  char tmp_file[FILENAME_MAX];
  int  ilen1, iret;
  char *ptr, *null;
  char delim[2];


  if ((tmp_str = strchr(arg, '=')) == NULL) {
    printf("Invalid argument \"%s\"\n", arg);
    exit(1);
  }

  ilen1 = tmp_str-arg;
  strncpy(keyword, arg, ilen1);
  keyword[ilen1] = '\0';
  strcpy(parm_str, tmp_str+1);

  /* change keyword to lower case */
  tmp_str = keyword;
  while (*tmp_str != '\0') {
    if (isupper(*tmp_str)) *tmp_str = tolower(*tmp_str);
    tmp_str++;
  }

  if (strncmp(keyword, "parfile",3) == 0) {
    iret = par_file(parm_str, input);
    return iret;
  
  } else if (strcmp(keyword, "infile") == 0) {
    parse_file_name(parm_str, tmp_file);
    strcpy(input->infile, tmp_file);

  } else if (strcmp(keyword, "ofile") == 0) {
    parse_file_name(parm_str, tmp_file);
    strcpy(input->ofile, tmp_file);

  } else if (strcmp(keyword, "fileuse") == 0) {
    parse_file_name(parm_str, tmp_file);
    strcpy(input->fileuse, tmp_file);

  } else if (strcmp(keyword, "sday") == 0) {
    input->sday = atoi(parm_str);

  } else if (strcmp(keyword, "eday") == 0) {
    input->eday = atoi(parm_str);

  } else if (strcmp(keyword, "resolve") == 0) {
    parse_file_name(parm_str, tmp_file);
    strcpy(input->resolve, tmp_file);

  } else if (strcmp(keyword, "rowgroup") == 0) {
    input->rowgroup = atoi(parm_str);

  } else if (strcmp(keyword, "flaguse") == 0) {

    input->flaguse[0] = '\0';

    null = 0;
    
    ptr = strtok(parm_str, ",");
    while ( ptr ) {
      if ( ptr ) {
	if ( strstr("default", ptr) ) {
	  if (input->flaguse[0] == '\0') {
	    strcpy(input->flaguse, DEF_FLAG);
	  } else {
	    strcat(input->flaguse, ",");
	    strcat(input->flaguse, DEF_FLAG);
	  }
	} else {
	  if (input->flaguse[0] == '\0') {
	    strcpy(input->flaguse, ptr);
	  } else {
	    strcat(input->flaguse, ",");
	    strcat(input->flaguse, ptr);
	  }
	}

      }
      ptr = strtok(null,",");
    }

  } else if (strcmp(keyword, "l3bprod") == 0) {

    strcpy(delim, ":");
    if (strchr(parm_str, ',') != NULL) strcpy(delim, ",");
    if (strchr(parm_str, ' ') != NULL) strcpy(delim, " ");

    strcpy(input->l3bprod, delim);
    strcat(input->l3bprod, parm_str);
    strcat(input->l3bprod, delim);

  } else if (strcmp(keyword, "prodtype") == 0) {
    strcpy(input->prodtype, parm_str);

  } else if (strcmp(keyword, "pversion") == 0) {
    strcpy(input->pversion, parm_str);

  } else if (strcmp(keyword, "suite") == 0) {
    strcpy(input->suite, parm_str);

  } else if (strcmp(keyword, "average") == 0) {
    strcpy(input->average, parm_str);

  } else if (strcmp(keyword, "oformat") == 0) {
    strcpy(input->oformat, parm_str);

  } else if (strcmp(keyword, "interp") == 0) {
    input->interp = atoi(parm_str);
    input->interp = 0; /* Disable interp */
    printf("INTERP parameter disabled.\n");

  } else if (strcmp(keyword, "latsouth") == 0) {
    input->latsouth = atoi(parm_str);

  } else if (strcmp(keyword, "latnorth") == 0) {
    input->latnorth = atoi(parm_str);

  } else if (strcmp(keyword, "lonwest") == 0) {
    input->lonwest = atoi(parm_str);

  } else if (strcmp(keyword, "loneast") == 0) {
    input->loneast = atoi(parm_str);

  } else if (strcmp(keyword, "meminfo") == 0) {
    input->meminfo = atoi(parm_str);

  } else if (strcmp(keyword, "dcinfo") == 0) {
    input->dcinfo = atoi(parm_str);

  } else if (strcmp(keyword, "noext") == 0) {
    input->noext = atoi(parm_str);

  } else if (strcmp(keyword, "night") == 0) {
    input->night = atoi(parm_str);

  } else if (strcmp(keyword, "verbose") == 0) {
    input->verbose = atoi(parm_str);

  } else if (strcmp(keyword, "healpix") == 0) {
    input->healpix = atoi(parm_str);

  } else if (strcmp(keyword, "minobs") == 0) {
    input->minobs = atoi(parm_str);

  } else if (strcmp(keyword, "deflate") == 0) {
    input->deflate = atoi(parm_str);

  } else if (strcmp(keyword, "qual_max") == 0) {
    input->qual_max = (uint8) atoi(parm_str);

  } else if (strcmp(keyword, "qual_prod") == 0) {
    parse_file_name(parm_str, tmp_file);
    strcpy(input->qual_prod, tmp_file);

  } else {
    goto Invalid_return;

  }

  return 0;


Invalid_return:
  printf("Invalid argument \"%s\"\n", arg);
  exit(1);
}


int par_file(char *pfile, instr *input)
{
  FILE *fp;
  char arg[2048];
  int  iret;

  if ((fp = fopen(pfile, "r")) == NULL) {
    printf("Error on opening the parameter file - %s\n", pfile);
    exit(1);
  }

  while ((fgets(arg, 2047, fp)) != NULL) {
    /* skip the comment or blank line */
    if (arg[0] == '#' || arg[0] == ';' || arg[0] == ' ' || arg[0] == '\0' ||
        arg[0] == '\n')   
      continue;

    arg[strlen(arg)-1] = '\0';    /* replace the last char new_line to NULL */
    iret = get_item(arg, input);
    if (iret != 0) {
      fclose(fp);
      return -1;
    }
  }

  fclose(fp);
  return 0;
}
