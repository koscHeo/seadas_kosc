/* These constants are specified in Fred Patt's cdata.f module. */

#define PI	3.1415926535897932384626433832795029L
#define RADEG	(180.0/PI)		/* For radians to degrees conversion */
#define RE	6378.137		/* Earth equatorial radius (km) */
#define REM	6371.0			/* Earth mean radius (km) */
#define FF	(1.0/298.257)		/* Earth flattening factor */
#define OMF2	((1.0 - FF)*(1.0 - FF))	/* One minus flattening, squared */
#define OMEGAE	7.29211585494e-5L	/* Earth rotation rate (radians/sec) */
