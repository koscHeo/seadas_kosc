README for determine_hico_boresight
2013/03/08 
Marcos Montes, marcos.montes@nrl.navy.mil, 202-767-7308

This folder contains the code that is used to make the executable for
determining the HICO IGMs, as well as solar zenith and azimuth, and view
zenith and azimuth - all for each pixel in the image. 

In order to make compile it on your computer(s), you need to first edit the
makefile in order to point to the corrct paths. The paths is this makefile
are for my Macbook Pro running OSX Lion. 

I compiled and ran it using gfortran. No fancy flags were used. 

Some of the code comes from
1) my Tafkaa library (default_character_length, auxtools,
generic_info,chartools.f90, SOLCOR.f90, SUNCOR.f90, solar_geometry.f90,
HAZEL.f90, file_io_new.f90, JULIAN.f90, check_range.f90, envi_header.f90,
errors.f90, write_envi_header.f90);
2) translations of IDL code I've used for a variety of projects (well, HICO)
(ecef2latlon.f90, q_to_r.f90, intercept.f90); 
3) Vallado's code that is available online (ASTMATH.for, ASTMATH.CMN,
ASTREDUC.FOR, ASTREDUC.CMN, ASTTIME.f)
4) And some code that I wrote just for this project (proc_hico.f90,
read_gcp.f90, bore_site.f90, get_astro_data.f90)

If you have gfortran installed, and if you change the paths so that they
work for your system, and if you have a good make installed, then:
prompt> make 
should compile and produce an executable.  

INPUTS NEEDED:
Inputs needed  (examples are in example_inputs) are 
a HICO pos-vel-quat file

Bore-sight parameter offsets, in degrees, determined using, say,
determine_hico_boresight. 

ADDITIONALLY NEEDED (to be placed in astro_data_dir):
Some time and earth orientation parameters from USNO these probably need to
be obtained once a day for long term processing. The files are:
finals_daily_all.txt
tai-utc.dat.txt

finals_daily can be found at:
http://www.usno.navy.mil/USNO/earth-orientation/eo-products/daily

tai-utc.dat.txt can be found at http://maia.usno.navy.mil/ser7/tai-utc.dat
(one link to this was at
http://www.usno.navy.mil/USNO/earth-orientation/eo-info/general/date-time-def)


CALLING/USAGE
 proc_hico_mjm /path/to/iss*pos_vel_quat.csv /path/to/usno_data_dir/ 
 assumes bore sight offsets are all 0.d0
  OR
 proc_hico_mjm /path/to/iss*pos_vel_quat.csv /path/to/usno_data_dir/ -0.900000 0.1 -0.01
 where the last three are XYZ bore sight offsets in degrees (if one is present, all three must be)
XYZ are about the XYZ axes, roughly roll-like, pitch-like, heading-like.

ASSUMPTIONS:
ORBIT LOCATIONS are perfect.
Attitudes are perfect and represent attitudes at HICO. 

In reality, the attitudes are the attitudes at the ISS attitude center, as determined by an array of GPS sensors. 

Part of the assumptions are that we understand everything perfectly about all the clocks involved. In order to be good enough to a frame, we really need better than ~13ms accuracy.

Assumes the expression for the lens pointing is perfect. 

OUTPUTS
a *.hico_LonLatViewAngles.bil file and its header. 



