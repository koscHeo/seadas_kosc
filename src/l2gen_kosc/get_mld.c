#include "l12_proto.h"
#include <sys/types.h>
#include <unistd.h>

// module getmld.c -retrieve mixed layer depth from World Ocean Atlas Climatology

static float mldbad = BAD_FLT;

#define MLDXANC 360
#define MLDYANC 180
/* This program opens the netcdf/hdf mld_climatology_woa1994 file located in ocssw//run/data/common. It extracts the mixed layer depth information from the hdf file. */
/* 	:Title = "WOA Mixed Layer Depth Monthly Climatology"
 :Description = "World Ocean Atlas 1994, Mixed Layer Depth, Monthly"
 :Reference = "http://www.nodc.noaa.gov/OC5/WOA94/mix.html"
 */

float get_mld(float lon, float lat, int day) {
	static int firstCall = 1;
	static int mldx = MLDXANC;
	static int mldy = MLDYANC;
	static float dx = 360.0 / MLDXANC;
	static float dy = 180.0 / MLDYANC;
	static float mldref[MLDYANC + 2][MLDXANC + 2];
	static float mld_sav[MLDYANC + 2][MLDXANC + 2];

	float mld = mldbad;
	int i, j, ii;
	float xx, yy;
	float t, u;

	if (firstCall) {

		float mldtmp[MLDYANC][MLDXANC];
		char name[H4_MAX_NC_NAME];
		char sdsname[H4_MAX_NC_NAME];
		int32 sd_id;
		int32 sds_id;
		int32 rank;
		int32 status;
		int32 sds_index;
		int32 nt;
		int32 dims[H4_MAX_VAR_DIMS];
		int32 nattrs;
		int32 start[2] = { 0, 0 };
		int32 ix, iy, ct;
		int32 lymin, lxmin, lymax, lxmax;
		float sum;

		char *OCSSW_root, mldclimatology_file[300];
		char mldmonth[7];

		firstCall = 0;

		if ((OCSSW_root = getenv("OCDATAROOT")) == NULL) {
			printf(
					"-E- %s:  Error looking up environmental variable OCDATAROOT\n",
					__FILE__);
			exit(1);
		}
		strcpy(mldclimatology_file, OCSSW_root);
		strcat(mldclimatology_file, "/common/mld_climatology_woa1994.hdf");

		sd_id = SDstart(mldclimatology_file, DFACC_RDONLY);
		if (sd_id == -1) {
			printf("-E- %s:  Error opening mld climatology file %s.\n",
					__FILE__, mldclimatology_file);
			exit(1);
		}

		printf("\nOpening mld climatology file %s\n\n", mldclimatology_file);

		/* calculate the correct date and day of the  month using yd2md utility. */

		int16 mon = (int) day / 31 + 1; // month of year (no need for perfection..at least according to the sea surface salinity reference algorithm)

		
		
		sprintf(mldmonth, "mld%02i", mon);
		

		/* Get the SDS index */
		strcpy(sdsname, mldmonth);
		sds_index = SDnametoindex(sd_id, sdsname);
		if (sds_index == -1) {
			printf("-E- %s:  Error seeking %s SDS from %s.\n", __FILE__,
					sdsname, mldclimatology_file);
			exit(1);
		}

		

		sds_id = SDselect(sd_id, sds_index);

		status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
		if (status != 0) {
			printf("-E- %s:  Error reading SDS info for %s from %s.\n",
					__FILE__, sdsname, mldclimatology_file);
			exit(1);
		}
		status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) &mldtmp[0][0]);
		if (status != 0) {
			printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__,
					sdsname, mldclimatology_file);
			exit(1);
		}

		status = SDendaccess(sds_id);
		status = SDend(sd_id);

		/* rotate 180-deg and add wrapping border to simplify interpolation */
		/* new grid is -180.5,180.5 in i=0,361 and -90.5,90.5 in j=0,181    */

		

		for (j = 0; j < mldy; j++) {
			for (i = 0; i < mldx; i++) {
				ii = (i < mldx / 2) ? i + mldx / 2 : i - mldx / 2;
				if (mldtmp[j][i] > 0)
					mldref[j + 1][ii + 1] = mldtmp[j][i];
				else
					mldref[j + 1][ii + 1] = mldbad;
			}
			mldref[j + 1][0] = mldref[j + 1][mldx];
			mldref[j + 1][mldx + 1] = mldref[j + 1][1];
		}
		for (i = 0; i < mldx + 2; i++) {
			mldref[0][i] = mldref[1][i];
			mldref[mldy + 1][i] = mldref[mldy][i];
		}

		/*
		 *  Extend the grid by 1 point (use only full grid boxes later)
		 */

		

		memcpy(mld_sav, mldref,
				( MLDYANC + 2) * ( MLDXANC + 2) * sizeof(float));
		for (iy = 0; iy < mldy + 2; iy++) {
			lymin = (iy == 0) ? 0 : iy - 1;
			lymax = (iy == mldy + 1) ? mldy + 1 : iy + 1;

			for (ix = 0; ix < mldx + 2; ix++) {
				if (mldref[iy][ix] < mldbad + 1) {
					lxmin = (ix == 0) ? 0 : ix - 1;
					lxmax = (ix == mldx + 1) ? mldx + 1 : ix + 1;

					sum = 0.;
					ct = 0;
					for (j = lymin; j < lymax + 1; j++)
						for (i = lxmin; i < lxmax + 1; i++) {
							if (mld_sav[j][i] > mldbad + 1) {
								sum += mld_sav[j][i];
								ct++;
							}
						}
					if (ct > 0)
						mldref[iy][ix] = sum / ct;
				}
			}
		}
		
	} /* end 1-time grid set up */


	/* locate LL position within reference grid */
	i = MAX(MIN((int) ((lon+180.0+dx/2)/dx),MLDXANC+1),0);
	j = MAX(MIN((int) ((lat+ 90.0+dy/2)/dy),MLDYANC+1),0);

	/* compute longitude and latitude of that grid element */
	xx = i * dx - 180.0 - dx / 2;
	yy = j * dy - 90.0 - dy / 2;

	/* Bilinearly interpolate, only using full grid boxes */
	t = (lon - xx) / dx;
	u = (lat - yy) / dy;

	if ((mldref[j][i] > mldbad + 1) && (mldref[j][i + 1] > mldbad + 1)
			&& (mldref[j + 1][i] > mldbad + 1)
			&& (mldref[j + 1][i + 1] > mldbad + 1)) {

		mld = (1 - t) * (1 - u) * mldref[j][i] + t * (1 - u) * mldref[j][i + 1]
				+ t * u * mldref[j + 1][i + 1] + (1 - t) * u * mldref[j + 1][i];

	} else
		mld = mldbad;

	return (mld);
}

