/*-----------------------------------------------------------------------------
    File : get_cal_misc.c

    Contents:
	get_index	-  reads time vdata and returns appropriate index to 
			   access data
	read_parm_data  -  reads parameter data

    Other relevant files:
	cal.h		-  various #defined constants, TDI table, and also
				includes hdf.h
	getcal_proto.h  -  prototypes for get_cal functions
	get_cal.c       -  a higher layer of calibration input functions

    Notes:

        Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      03/11/94    Original development
        Lakshmi Kumar    Hughes STX      06/07/94    Updated to reflect v3.1
                                                        interface spcifications
	Lakshmi Kumar    Hughes STX	 03/21/96    Corrected non-prototype
						     declarations
	Lakshmi Kumar	 Hughes STX	 03/17/97    Removed non-ANSI proto
						     declarations.  In-detector
						     offsets are redefined as
						     idoffs[8][16]. 
------------------------------------------------------------------------------*/

#include "get_cal_osmi.h"
#include "getcal_proto_osmi.h"
#include <hdf4utils.h>

/*-----------------------------------------------------------------------------
    Function: get_index 

    Returns: int32 (index)
 	On successful it returns the index of the given time entry and if
	time entry not found, returns -3.

    Description:
        The function get_index reads time vdata and searches for the
	given time entry.  If the given time found, it rerurns the entry
	number to access related information from slopes and parameter vdatas.

    Arguments: (in calling order)
      Type       Name        I/O     Description
      ----       ----        ---     -----------
      int32      fid          I      file ID
      int16      syear        I      year of data start time 
      int16      sday         I      day-of-year for data start time 
      int16      eday         I      day-of-year for data end time
      int32      smsec        I      milliseconds-of-day for data start time
      int16      *cal_year    O      year the cal entry was made
      int16      *cal_day     O      day of the year the cal entry was made

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      03/11/94    Original development
	Lakshmi Kumar    Hughes STX      06/07/94    Updated to reflect v3.1
                                                        interface spcifications
   	Lakshmi Kumar	 Hughes STX	 02/07/94    Added code to return 
						     cal entry year and day
						     (ref to I/O specs v4.2)
------------------------------------------------------------------------------*/

int32 get_index_osmi(int32 fid, int16 syear, int16 sday, int16 eday, int32 msec,
		     int16 *cal_year, int16 *cal_day, int32 *cal_msec)
{

  int16 dyear, dday, *cal_syear, *cal_sday, *cal_eyear, *cal_eday;
  int16 *entry_year, *entry_day;
  int32 i, *cal_smsec, *cal_emsec, vsid, elts;


  if ((vsid = attach_vdata(fid, TIME)) < 0)
     return RDERR;

  if ((elts = VSelts(vsid)) < 0)
     return RDERR;
  if((cal_syear = (int16 *) malloc(elts * sizeof(int16)))==NULL)
          return BUFERR;

  if((cal_eyear = (int16 *) malloc(elts * sizeof(int16)))==NULL)
          return BUFERR;

  if((cal_sday  = (int16 *) malloc(elts * sizeof(int16)))==NULL)
          return BUFERR;

  if((cal_eday  = (int16 *) malloc(elts * sizeof(int16)))==NULL)
          return BUFERR;

  if((cal_smsec = (int32 *) malloc(elts * sizeof(int32)))==NULL)
          return BUFERR;

  if((cal_emsec = (int32 *) malloc(elts * sizeof(int32)))==NULL)
          return BUFERR;

  if((entry_year = (int16 *) malloc(elts * sizeof(int16)))==NULL)
          return BUFERR;

  if((entry_day = (int16 *) malloc(elts * sizeof(int16)))==NULL)
          return BUFERR;

  rdvdata(vsid, SYEAR, 0, elts, (unsigned char *)cal_syear);
  rdvdata(vsid, SDAY,  0, elts, (unsigned char *)cal_sday);
  rdvdata(vsid, SMSEC, 0, elts, (unsigned char *)cal_smsec);

  rdvdata(vsid, EYEAR, 0, elts, (unsigned char *)cal_eyear);
  rdvdata(vsid, EDAY,  0, elts, (unsigned char *)cal_eday);
  rdvdata(vsid, EMSEC, 0, elts, (unsigned char *)cal_emsec);

  rdvdata(vsid, ENTRY_YEAR, 0, elts, (unsigned char *)entry_year);
  rdvdata(vsid, ENTRY_DAY,  0, elts, (unsigned char *)entry_day);

  dyear = syear;
  dday = sday;
  if (sday != eday && msec < 43200000)
     dday = eday;
  if (dday < sday)
     dyear += 1;

  for(i = elts-1; i >= 0; i--) {
    if (cal_eyear[i] == 0){ /* onwards rec */
       if (dyear > cal_syear[i])
          break;
       if (dyear == cal_syear[i] && dday > cal_sday[i])
          break;
       if (dyear == cal_syear[i] && dday == cal_sday[i] &&
                msec >= cal_smsec[i])
          break;
     }
    else {  /* not an onwards rec */
       if (dyear > cal_syear[i]) {
          if (dyear < cal_eyear[i])
             break;
          if (dyear == cal_eyear[i]) {
             if (dday < cal_eday[i])
                break;
             if (dday == cal_eday[i] && msec <= cal_emsec[i])
                break;
           }
        }
       else
         if (dyear == cal_syear[i]) {
            if (dyear < cal_eyear[i])
               break;
            if (dyear == cal_eyear[i]) {
               if (dday >= cal_sday[i] && dday < cal_eday[i])  
                  break;
               if (dday >= cal_sday[i] && dday == cal_eday[i] && 
			msec <= cal_emsec[i])
                  break;
             }
          }
      }
   }

  *cal_year = cal_syear[i];
  *cal_day  = cal_sday[i];
  *cal_msec  = cal_smsec[i];

  VSdetach(vsid);
  free(cal_syear);
  free(cal_sday);
  free(cal_smsec);
  free(cal_eyear);
  free(cal_eday);
  free(cal_emsec);
  free(entry_year);
  free(entry_day);

  if (i < 0)
     return TMERR;
  else
     return i;
}

/*-----------------------------------------------------------------------------
    Function: read_parm_data

    Returns: int32 (Status)
	On success returns 0, otherwise returns -2 indicating read error. 
    Description:
        The function  attaches to the requested vdata

    Arguments: (in calling order)
      Type       Name        I/O     Description
      ----       ----        ---     -----------
      int32      fid          I      HDF file ID
      int32      sdfid        I      HDF SD file ID
      int32      index        I      element number
      int32      idoffs[8][16] O      detector zero offset counts        
      float32    gains[8][16] O      slopes (band*detector*gains)
      float32    temps[256][2]O      temp and ref temp correction factors
      float32    scan_mod[2][1285]O  scan-modulation buffer
      float32    mirror[8][2] O      mirror side1 and side2 for all bands
      int16      tdi_list[256][4] O  TDI values for all 8 bands

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      03/11/94    Original development
	Lakshmi Kumar    Hughes STX      06/07/94    Updated to reflect v3.1
                                                        interface spcifications
------------------------------------------------------------------------------*/
int32 read_parm_data_osmi(int32 fid, int32 sdfid, int32 index, float32 *eoffset, 
			  float32 *egain, float32 *temp, float32 *mirror,  
			  float32 *t_const, float32 *t_linear, float32 *t_quadratic, 
			  float32 *slope, float32 *dc, float32 *sm)
{

  int32 i, j, slpid, parmid;
  float32  egain_buf[9];
     if ((slpid = attach_vdata(fid, "Electronics")) < 0)
        return RDERR; 
     if ((rdvdata(slpid, gain_flds, index, 1, (unsigned char *)egain_buf)) < 0)
		return RDERR;
     VSdetach(slpid);
     for(i = 0; i < 8; i++)
     	egain[i] = egain_buf[i];
     *eoffset = egain_buf[8];
     if ((slpid = attach_vdata(fid, "Mirror")) < 0)
        return RDERR; 
     if ((rdvdata(slpid, mirror_flds, index, 1, (unsigned char *)mirror)) < 0)
		return RDERR;
     VSdetach(slpid);

     if ((slpid = attach_vdata(fid, "Longterm")) < 0)
        return RDERR; 
     if ((rdvdata(slpid, t_const_flds, index, 1, (unsigned char *)t_const)) < 0)
		return RDERR;
     if ((rdvdata(slpid, t_linear_flds, index, 1, (unsigned char *)t_linear)) < 0)
		return RDERR;
     if ((rdvdata(slpid, t_quadratic_flds, index, 1, (unsigned char *)t_quadratic)) < 0)
		return RDERR;
     VSdetach(slpid);
     if ((slpid = attach_vdata(fid, "TempCorr")) < 0)
        return RDERR; 
     if ((rdvdata(slpid, temp_flds, index, 1, (unsigned char *)temp)) < 0)
		return RDERR;
     VSdetach(slpid);
     if ((slpid = attach_vdata(fid, "Slopes")) < 0)
        return RDERR; 
     if ((rdvdata(slpid, slope_flds, index, 1, (unsigned char *)slope)) < 0)
		return RDERR;
     VSdetach(slpid);
     if ((slpid = attach_vdata(fid, "DCs")) < 0)
        return RDERR; 
     if ((rdvdata(slpid, dc_flds, index, 1, (unsigned char *)dc)) < 0)
		return RDERR;
     VSdetach(slpid);
     if ((slpid = attach_vdata(fid, "ScanMod")) < 0)
        return RDERR; 
     if ((rdvdata(slpid, sm_flds, index, 1, (unsigned char *)sm)) < 0)
		return RDERR;
     VSdetach(slpid);

  return SUCCEED;
}

