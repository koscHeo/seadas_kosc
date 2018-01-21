# 
# Clear all conditional flags
unset sgi_mode
unset pgs_daac
unset pgs_f90_comp
unset pgs_nag_flag
set pgs_absoft_flag=0
setenv fc_path /usr/bin
# set the toolkit home directory environment variable
# SDP installation done on Sun Jan 21 18:28:32 KST 2018 
# 
 
setenv PGSHOME /home/jaemoo/seadas-7.4/ocssw/build/lib3/src/sdptk/TOOLKIT       # set the toolkit home directory
set sgi_mode=32                   # set SGI for standard mode
set LINUX_BRAND=linux64		# set LINUX for -64 mode
setenv AA_install 0    #set option for AA tool instalation
 
#-----------------------------------------------------------------------------
# file:	
# 	pgs-dev-env.csh
#
# description:
# 	This file defines environment variables used by SDP (PGS) Toolkit.
# 	This version is for use under the C shell (csh).
#
# usage:
# 	This file should be called from your .cshrc file with the line:
#
# 	    source <toolkit-home-dir>/bin/pgs-dev-env.csh
#	
# 	where <toolkit-home-dir> is the full path of the toolkit home directory.
#
# author:
# 	Mike Sucher / A.R.C.
#       Guru Tej S. Khalsa / Applied Research Corporation
#	Abe Taaheri / Emergent Information Technologies, Inc.
#
# notes:
# 	1) This file is compatible with the following platforms:
# 	   Sun, SGI, HP-9000, IBM RS-6000, DEC Alpha and Cray.
# 	   It automatically figures out which platform you are on,
# 	   and sets environment variables accordingly.
# 	2) (OPTIONAL) if HDF is installed on your system, you may 
# 	   wish to edit the line near the end of this file that sets
# 	   up the environment variable HDFHOME and HDF5HOME.
# 	3) This file adds the PGS bin (and optionally the HDF bin) 
#	   directory to the 'path'  environment variable.
# 	4) This file defines a variable called pgs_path which contains
# 	   all the directories likely to be needed on a given machine
# 	   type, including the PGS (and HDF) bin directories.  Users 
# 	   may choose to set their path variable with:
# 	   
# 	           set path=(<user-path-additions> $pgs_path )
#
# 	   where <user-path-additions> is an optional list of other
# 	   directories added to the search path.
#
# history:
#	20-Apr-1994 MES  Initial version
#	31-Aug-1994 MES  Update for Release 3
# 	13-Sep-1994 MES  Added variables for Process Control tools
# 	14-Sep-1994 MES  Updated path for sun4 to include 
# 			   /usr/local/lang and /usr/lang
# 	26-Sep-1994 MES  Updated variable PGS_PC_INFO_FILE to point to the file
# 			 $PGSRUN/PCF.v3
# 	09-Feb-1995 MES  Updated variable PGS_PC_INFO_FILE to point to the file
# 			 $PGSRUN/PCF.v4
# 	03-May-1995 MES  Updated variable PGS_PC_INFO_FILE to point to the file
# 			 $PGSRUN/PCF.v5
# 			 Removed referenced to obsolete environment variables
# 			 no longer referenced by the Process Control tools.
# 	09-May-1995 MES  Updated compiler flags for Release 5
# 			   hp:                        ibm:
# 			   C_CFH:   remove -Dextname  C_CFH:   remove -Dextname
# 			   F77_CFH: remove +ppu       F77_CFH: remove -qextname
# 			 Note: this change was made to facilitate binding with
# 			 heritage code libraries.  Caution:  the cfortran.h 
# 			 bindings can ONLY be used with mixed-case C function
# 			 names.  Otherwise the bindings will not work.
#	15-Jun-1995 GTSK Added support for IRIX64, the sgi 64 bit O/S.
#			 This script will install the toolkit in 32 bit mode.
#	10-Jul-1995 MES  Removed flags "-h stdc" from Cray CFLAGS.
#	18-Aug-1995 MES  Revised environment for formal track directory
# 			 structure used by the DAAC version Toolkit for IR-1.
# 	31-Aug-1995 MES  Updated IRIX64 OS support to recognize both 32 and 64
# 			 bit modes.  Default is is 32-bit.  New BRAND setting
# 			 sgi64 supports the SGI 64 bit mode.
# 	30-Oct-1995 MES  Updates for IR-1 release:
# 			 - add support for SGI compiled -n32 ($BRAND = sgi32)
# 			 - formal directory structure is now set by pgs_formal
# 			 - revise PGSRUN: no longer architecture-specific
#  			 - set PGS_PC_INFO_FILE to $PGSRUN/$BRAND/PCF.ir1
# 	21-Nov-1995 MES  Handle the DAAC toolkit compilation.
# 	22-Nov-1995 MES  Add code to handle FORTRAN-90 compiler.
# 	15-Apr-1996 MES  Updates for Release A
#  			 - clear DAAC version env variables, if not used
#  			 - set PGS_PC_INFO_FILE to $PGSRUN/$BRAND/PCF.relA
# 	26-Aug-1996 MES  Updates for Release A
#  			 - update DAAC version logic
#  			 - update HDF environment variables setup
#       30-Mar-2000 AT   -update for IRIX 6.5
#       28-Feb-2001 AT   changed f77 to g77 for linux
#       30-Mar-2001 AT   Modified for HDF5 support 
#       16-Aug-2001 AA   Modified for solaris8
#       10-Ari-2003 PTN  Modified for pgf77 compiler 
#       12-Oct-2004 MP   Changed $HDFHOME/zlib to $ZLIBHOME/zlib,
#                        $HDFHOME/jpeg to $JPEGHOME/jpeg, $HDFHOME/szip to
#                        $SZIPHOME/szip (NCR_41266)
#-----------------------------------------------------------------------------

if (! $?sgi_mode ) then		# by default, SGI mode is standard 32-bit
    set sgi_mode=32
endif
if (! $?pgs_formal ) then 	# by default, DAAC toolkit version disabled
    set pgs_formal=0
endif
if (! $?pgs_daac ) then 	# by default, DAAC toolkit version disabled
    set pgs_daac=0
endif
if (! $?pgs_f90_comp ) then 	# by default, no FORTRAN-90 compiler
    set pgs_f90_comp=""
endif
if (! $?pgs_nag_flag ) then 	# by default, not using NAG FORTRAN-90
    set pgs_nag_flag=0
endif

if (! $?use_flavor ) then 	# by default, do not use "flavors"
    set use_flavor=0
endif
if (! $?use_flavor_r ) then 	# by default, do not use "flavors"
    set use_flavor_r=0
endif
if (! $?pgs_c_rlib ) then       # by default, not using NAG FORTRAN-90
    set pgs_c_rlib=0
endif

set user_path = ( $path )	# save user path

# set path to a base subset of directories, allowing startup on unknown host
# note: once the host has been determined the path is appropriately customized

set path=(/usr/local/bin /bin /usr/bin /etc /usr/etc /usr/ucb /usr/bin/X11)

if ($?LD_LIBRARY_PATH ) then
setenv LD_LIBRARY_PATH "/home/jaemoo/seadas-7.4/ocssw/build/lib3/lib:/home/jaemoo/seadas-7.4/ocssw/build/lib3/lib:$LD_LIBRARY_PATH"
else
setenv LD_LIBRARY_PATH "/home/jaemoo/seadas-7.4/ocssw/build/lib3/lib:/home/jaemoo/seadas-7.4/ocssw/build/lib3/lib"
endif

# get operating system type, login name
# special cases: SCO and Cray  - uname works differently,

setenv MACHINE `uname -m | awk '{print $1}'`	# needed on Cray & SCO
setenv temp_ostype `uname`

switch ( "$MACHINE" )

    case "i386":        			# SCO box
	setenv OSTYPE sco386
        breaksw

    case "CRAY":    				# CRAY
	setenv OSTYPE UNICOS
        breaksw

    default:					# everybody else
	setenv OSTYPE `uname`
        breaksw

endsw

setenv CYGPL `uname | awk -F_ '{print $1}'`

# Intel Macintosh is also i386 or i686 (?) machine

    if ( "$MACHINE" == "i386" ) then
	if ( "$temp_ostype" == "Darwin" ) then 
	    setenv OSTYPE DarwinIntel
	endif
	if("$CYGPL" == "CYGWIN") then 
	    setenv OSTYPE Cygwin
	endif
    endif
    if ( "$MACHINE" == "i686" ) then
	if ( "$temp_ostype" == "Darwin" ) then 
	    setenv OSTYPE DarwinIntel
	endif
	if("$CYGPL" == "CYGWIN") then 
	    setenv OSTYPE Cygwin
	endif
    endif


set user=`id | cut -d\( -f2 | cut -d\) -f1`
if ($?LOGNAME == 0) setenv LOGNAME $user	# make sure $LOGNAME is defined
setenv USER $LOGNAME				# make sure $USER is defined


# set machine-dependent environment variables:
# 	HOST  the host name of this machine
# 	BRAND used by other achitecture-specific code
# 	path  the execution search path exported to PATH

switch ( "$OSTYPE" )

    case "AIX": 
	set path=(/usr/local/bin /bin /usr/bin /etc /usr/etc /usr/ucb /usr/bin/X11 /usr/ccs/bin /usr/sbin)
	setenv HOST `hostname`
	setenv BRAND ibm
        breaksw

    case "HP-UX": 
	set path=(/usr/local/bin /bin /usr/bin /etc /usr/etc /usr/bin/X11)
	setenv HOST `hostname`
	setenv BRAND hp 
        breaksw

    case "IRIX":  
	set path=(/usr/local/bin /bin /usr/bin /etc /usr/etc /usr/bsd /usr/bin/X11 /usr/sbin)
	setenv HOST `hostname`
	setenv BRAND sgi 
        breaksw

    case "IRIX64":  
	set path=(/usr/local/bin /bin /usr/bin /etc /usr/etc /usr/bsd /usr/bin/X11 /usr/sbin)
	setenv HOST `hostname`
	if ("$sgi_mode" == 32) then
	    setenv BRAND sgi
        else if ("$sgi_mode" == 64) then
	    setenv BRAND sgi64
        else if ("$sgi_mode" == n32) then
	    setenv BRAND sgi32
	else if ("$sgi_mode" == 65) then
	    setenv BRAND irix65
        else
	    setenv BRAND sgi
        endif
        breaksw

    case "Linux": 
	set path=(/usr/local/bin /bin /usr/bin /etc /usr/etc /usr/X11/bin)
	setenv HOST `hostname`
	setenv BRAND linux
	if("$LINUX_BRAND" == "linux32") then
		 setenv BRAND linux32
	else if("$LINUX_BRAND" == "linux") then
		 setenv BRAND linux
	else if("$LINUX_BRAND" == "linux64") then
		 setenv BRAND linux64
	else if( "`uname -m | awk '{print $1}'`" == "x86_64" || "`uname -m | awk '{print $1}'`" == "ia64" ) then
		 setenv BRAND linux64
	else
		 setenv BRAND linux
	endif
        breaksw

    case "Darwin":
        set path=(/bin /sbin /usr/bin /usr/sbin)
        setenv HOST `hostname -s`
        setenv BRAND macintosh
        breaksw

    case "DarwinIntel":
        set path=(/bin /sbin /usr/bin /usr/sbin)
        setenv HOST `hostname -s`
        setenv BRAND macintel32
	if("$MAC_BRAND" == "macintel32") then 
		setenv BRAND macintel32
	else if("$MAC_BRAND" == "macintel") then 
		setenv BRAND macintel
	else if("$MAC_BRAND" == "macintel64") then
		setenv BRAND macintel64
	else if( "`uname -m | awk '{print $1}'`" == "x86_64" ||  "`uname -m | awk '{print $1}'`" == "ia64"   || "`uname -m | awk '{print $1}'`" == "i386"   || "`uname -m | awk '{print $1}'`" == "i686" ) then
		setenv BRAND macintel64
		setenv MAC_BRAND macintel64
	else
		setenv BRAND macintel32

	endif
        breaksw

    case "Cygwin":
        set path=(/bin /sbin /usr/bin /usr/sbin)
        setenv HOST `hostname`
        setenv BRAND cygwin
        breaksw

    case "OSF1":  
	set path=(/usr/local/bin /bin /usr/bin /etc /usr/etc /usr/ucb /usr/bin/X11 /usr/ccs/bin /usr/sbin)
	setenv HOST `hostname -s`
	setenv BRAND dec 
        breaksw

    case "sco386": 
	set path=(/usr/local/bin /bin /usr/bin /etc /usr/etc /usr/bin/X11)
	setenv HOST `hostname -s`
	setenv BRAND sco 
        breaksw

    case "SunOS": 
	# distinguish between SunOS 5.x and 4.x versions
	set cbrand=`uname -r | awk -F. '{print $1, $2}'`
	if ("$cbrand" == "5 10") then 
	    setenv BRAND sun5.10		# release V5.10 SunOS
	    set path=(/usr/local/bin /opt/SUNWspro/bin /bin /usr/bin /etc /usr/etc /usr/ucb /usr/openwin/bin /usr/openwin/demo /usr/ccs/bin /usr/sbin)
        else if ("$cbrand" == "5 9") then 
	    setenv BRAND sun5.9			# release V5.9 SunOS
	    set path=(/usr/local/bin /opt/SUNWspro/bin /bin /usr/bin /etc /usr/etc /usr/ucb /usr/openwin/bin /usr/openwin/demo /usr/ccs/bin /usr/sbin)
        else if ("$cbrand" == "5 8") then 
	    setenv BRAND sun5.8			# release V5.8 SunOS
	    set path=(/usr/local/bin /opt/SUNWspro/bin /bin /usr/bin /etc /usr/etc /usr/ucb /usr/openwin/bin /usr/openwin/demo /usr/ccs/bin /usr/sbin)
        else if ("$cbrand" == "5 5") then
            setenv BRAND sun5        # release V5.x SunOS
            set path=(/usr/local/bin /opt/SUNWspro/bin /bin /usr/bin /etc /usr/etc /usr/ucb /usr/openwin/bin /usr/openwin/demo /usr/ccs/bin /usr/sbin)
         else
            setenv BRAND sun4			# release V4.x SunOS
	    set path=(/usr/local/bin /usr/local/lang /usr/lang /bin /usr/bin /etc /usr/etc /usr/ucb /usr/openwin/bin /usr/openwin/demo) 
          endif

	setenv HOST `hostname`
        breaksw

    case "UNICOS": 
	set path=(/usr/local/bin /bin /usr/bin /etc /usr/ucb /usr/bin/X11)
	setenv HOST `hostname`
	setenv BRAND cray 
        breaksw

    default:
	echo "Operating system: $OSTYPE not supported"
	echo "This release of the PGS Toolkit supports: "
	echo "   Sun, SGI HP-9000 IBM-6000 DEC-Alpha and Cray/Unicos "
        breaksw

endsw




# set machine-dependent compilers and compilation switches:
#
#

switch ($BRAND)

    case cray:
	setenv CC cc 			# C compiler
	setenv CFLAGS "-O2 -DH5_USE_16_API " 		# default C flags (optimize, ansi)
	setenv C_CFH "-DCRAYFortran"    # C/cfortran.h called from FORTRAN
	setenv CFHFLAGS "$CFLAGS $C_CFH" # CFLAGS + C_CFH
	setenv F77 cf77			# FORTRAN compiler
	setenv F77FLAGS "-O2 -DH5_USE_16_API " 		# common FORTRAN flags
        setenv F77_CFH ""               # FORTRAN callable from C w/ cfortran.h
	setenv F77_C_CFH ""          	# calling C w/ cfortran.h
	setenv F77_C_LIB ""		# C lib called by FORTRAN main 
	setenv HDFSYS UNICOS		# system type as defined by HDF
	breaksw

    case dec:
	setenv CC cc 			# C compiler
	setenv CFLAGS "-O2 -DH5_USE_16_API  -std"		# default C flags (optimize, ansi)
	setenv C_CFH "-DDECFortran"	# C w/ cfortran.h callable from FORTRAN
	setenv CFHFLAGS "$CFLAGS $C_CFH" # CFLAGS + C_CFH
	setenv F77 f77 			# FORTRAN compiler
	setenv F77FLAGS "-O2 -DH5_USE_16_API " 		# common FORTRAN flags
        setenv F77_CFH ""               # FORTRAN callable from C w/ cfortran.h
	setenv F77_C_CFH ""             # calling C w/ cfortran.h
	setenv F77_C_LIB ""		# C lib called by FORTRAN main
	setenv HDFSYS DEC_ALPHA		# system type as defined by HDF
	breaksw

    case hp:
	setenv CC cc 			# C compiler
	setenv CFLAGS "-O2 -DH5_USE_16_API  -Ae" 		# default C flags (optimize, ansi)
	setenv C_CFH "" 		# C w/ cfortran.h callable from FORTRAN
	setenv CFHFLAGS "$CFLAGS $C_CFH" # CFLAGS + C_CFH
	setenv F77 fort77		# FORTRAN compiler
	setenv F77FLAGS "-O2 -DH5_USE_16_API " 		# common FORTRAN flags
        setenv F77_CFH ""               # FORTRAN callable from C w/ cfortran.h
	setenv F77_C_CFH ""          	# calling C w/ cfortran.h
	setenv F77_C_LIB ""		# C lib called by FORTRAN main
	setenv HDFSYS HP9000		# system type as defined by HDF
	breaksw

    case ibm:
	setenv CC cc 			# C compiler
	setenv CFLAGS "-O2 -DH5_USE_16_API  -qlanglvl=ansi" # default C flags (optimize, ansi)
	setenv C_CFH "" 		# C w/ cfortran.h callable from FORTRAN
	setenv CFHFLAGS "$CFLAGS $C_CFH" # CFLAGS + C_CFH
	setenv F77 xlf 			# FORTRAN compiler
	setenv F77FLAGS "-O2 -DH5_USE_16_API " 		# common FORTRAN flags
        setenv F77_CFH ""               # FORTRAN callable from C w/ cfortran.h
	setenv F77_C_CFH ""          	# calling C w/ cfortran.h
	setenv F77_C_LIB ""		# C lib called by FORTRAN main
	setenv HDFSYS IBM6000		# system type as defined by HDF
	breaksw

    case linux:
    case linux32:
    case linux64:
	setenv CC "gcc " 			# C compiler
	setenv CFLAGS "-O2 -DH5_USE_16_API  -ansi"	# default C flags (optimize, ansi)
	setenv C_CFH "-Df2cFortran"	# C w/ cfortran.h callable from FORTRAN
	setenv CFHFLAGS "$CFLAGS $C_CFH" # CFLAGS + C_CFH
        if ( "`fgrep pgf77 $PGSHOME/bin/$BRAND/pgs-env.csh`" != "" ) then
	    setenv F77 pgf77                    # pgf77 Fortran compiler
        else
	    setenv F77 "g77 "			# FORTRAN compiler - default
        endif

	if ( "`echo $pgs_f90_comp | grep pgf90`" != "" || "$F77" == "pgf77" ) then
		setenv F77FLAGS "-O2 -DH5_USE_16_API " # common FORTRAN flags
        	setenv F77_CFH ""         # FORTRAN callable from C w/ cfortran.h
	else
		setenv F77FLAGS "-O2 -DH5_USE_16_API  -fno-second-underscore"
        	setenv F77_CFH "-fno-second-underscore"
	endif

	setenv F77_C_CFH ""          	# calling C w/ cfortran.h
	setenv F77_C_LIB ""		# C lib called by FORTRAN main
	if( "`uname -m | awk '{print $1}'`" == "x86_64" ) then
	    if( "$BRAND" == "linux64" ) then
		setenv HDFSYS LINUX64      # system type as defined by HDF
	    else if( "$BRAND" == "linux32" || "$BRAND" == "linux" ) then
		setenv HDFSYS LINUX
	    endif
	else if( "`uname -m | awk '{print $1}'`" == "ia64" ) then
	    if( "$BRAND" == "linux64" ) then
		setenv HDFSYS IA64      # system type as defined by HDF
	    else if( "$BRAND" == "linux32" || "$BRAND" == "linux" ) then
		setenv HDFSYS LINUX
	    endif
	else
	    setenv HDFSYS LINUX		   # system type as defined by HDF
	endif
	breaksw

    case macintosh:
        setenv CC gcc                   # C compiler
        setenv CFLAGS "-O2 -DH5_USE_16_API  -DMACINTOSH" # default C flags (optimize, ansi)
        setenv C_CFH "-Df2cFortran"     # C w/ cfortran.h callable from FORTRAN
        setenv CFHFLAGS "$CFLAGS $C_CFH"# CFLAGS + C_CFH
        setenv F77 "g77"                # FORTRAN compiler - default
        setenv F77FLAGS "-O2 -DH5_USE_16_API  -fno-second-underscore" # common FORTRAN flags  
        setenv F77_CFH "-fno-second-underscore"      # FORTRAN callable from C w/cfortran.h
        setenv F77_C_CFH ""             # calling C w/ cfortran.h
        setenv F77_C_LIB ""             # C lib called by FORTRAN main
        setenv HDFSYS MACINTOSH         # system type as defined by HDF
        breaksw

    case macintel:
    case macintel32:
    case macintel64:
        setenv CC "gcc MACINTEL_CMP_FLAG"               # C compiler
        setenv CFLAGS "-O2 -DH5_USE_16_API  -DMACINTEL"  # default C flags (optimize, ansi)
        setenv C_CFH "-Df2cFortran"     # C w/ cfortran.h callable from FORTRAN
        setenv CFHFLAGS "$CFLAGS $C_CFH"# CFLAGS + C_CFH
        setenv F77 "gfortran MACINTEL_CMP_FLAG"         # FORTRAN compiler - default
        setenv F77FLAGS "-O2 -DH5_USE_16_API  -fno-second-underscore" # common FORTRAN flags  
        setenv F77_CFH "-fno-second-underscore"      # FORTRAN callable from C w/cfortran.h
        setenv F77_C_CFH ""             # calling C w/ cfortran.h
        setenv F77_C_LIB ""             # C lib called by FORTRAN main
        setenv HDFSYS MACINTEL          # system type as defined by HDF
        breaksw

    case cygwin:
        setenv CPP g++                  # C++ compiler
        setenv CC gcc                   # C compiler
        setenv CFLAGS "-O2 -DH5_USE_16_API  -DCYGWIN"     # default C flags (optimize, ansi)
        setenv C_CFH "-Df2cFortran"     # C w/ cfortran.h callable from FORTRAN
        setenv CFHFLAGS "$CFLAGS $C_CFH"  # CFLAGS + C_CFH

        setenv CPPFLAGS "-O2 -DH5_USE_16_API  -DCYGWIN"   # default C ++ flags (optimize)
        setenv CPP_CFH "-Df2cFortran"   # C ++ w/ cfortran.h callable from FORTRAN
        setenv CPPFHFLAGS "$CPPFLAGS $CPP_CFH"      # CFLAGS + C_CFH for C++
        setenv F77 "g77"                # FORTRAN compiler - default
        setenv F77FLAGS "-O2 -DH5_USE_16_API  -fno-second-underscore" # common FORTRAN flags  
        setenv F77_CFH "-fno-second-underscore"     # FORTRAN callable from C w/cfortran.h
        setenv F77_C_CFH ""             # calling C w/ cfortran.h
        setenv F77_C_LIB ""             # C lib called by FORTRAN main
        setenv HDFSYS CYGWIN            # system type as defined by HDF
        breaksw

    case sco:
	setenv CC cc 			# C compiler
	setenv CFLAGS "-O2 -DH5_USE_16_API  -posix"	# default C flags (optimize, ansi)
	setenv C_CFH "-Df2cFortran"	# C w/ cfortran.h callable from FORTRAN
	setenv CFHFLAGS "$CFLAGS $C_CFH" # CFLAGS + C_CFH
	setenv F77 ""			# FORTRAN compiler
	setenv F77FLAGS "-O2 -DH5_USE_16_API " 		# common FORTRAN flags
        setenv F77_CFH ""               # FORTRAN callable from C w/ cfortran.h
	setenv F77_C_CFH ""          	# calling C w/ cfortran.h
	setenv F77_C_LIB ""		# C lib called by FORTRAN main
	setenv HDFSYS SCO		# system type as defined by HDF
	breaksw

    case sgi:
	if ($OSTYPE == "IRIX64") then
		setenv CC "cc -32"	# C compiler (32 bit)
		setenv F77 "f77 -32"	# FORTRAN compiler (32 bit)
	else
                setenv CC cc		# C compiler
                setenv F77 f77		# FORTRAN compiler
	endif
	setenv CFLAGS "-O2 -DH5_USE_16_API  -xansi -D_POSIX_SOURCE"	# default C flags (optimize, ansi)
	setenv C_CFH ""	 		# C w/ cfortran.h callable from FORTRAN
	setenv CFHFLAGS "$CFLAGS $C_CFH" # CFLAGS + C_CFH
	setenv F77FLAGS "-O2 -DH5_USE_16_API " 		# common FORTRAN flags
        setenv F77_CFH ""               # FORTRAN callable from C w/ cfortran.h
	setenv F77_C_CFH ""          	# calling C w/ cfortran.h
	setenv F77_C_LIB ""		# C lib called by FORTRAN main
	setenv HDFSYS IRIX		# system type as defined by HDF
	breaksw

    case sgi32:
	setenv CC "cc -n32"		# C compiler (new-style 32 bit)
	setenv F77 "f77 -n32"		# FORTRAN compiler (new-style 32 bit)
	setenv CFLAGS "-O2 -DH5_USE_16_API  -xansi -D_POSIX_SOURCE"	# default C flags (optimize, ansi)
	setenv C_CFH ""	 		# C w/ cfortran.h callable from FORTRAN
	setenv CFHFLAGS "$CFLAGS $C_CFH" # CFLAGS + C_CFH
	setenv F77FLAGS "-O2 -DH5_USE_16_API " 		# common FORTRAN flags
        setenv F77_CFH ""               # FORTRAN callable from C w/ cfortran.h
	setenv F77_C_CFH ""          	# calling C w/ cfortran.h
	setenv F77_C_LIB ""		# C lib called by FORTRAN main
	setenv HDFSYS IRIX		# system type as defined by HDF
	breaksw

    case sgi64:
	set cpu_type=`hinv | fgrep CPU | head -1 | cut -d' ' -f3 | cut -b2`
	if ("$cpu_type" == "4") then
                setenv CC "cc -64 -mips3"       # C compiler (R4?00 chip)
                setenv F77 "f77 -64 -mips3"     # FORTRAN compiler (R4?00 chip)
        else
                setenv CC "cc -64"      # C compiler
                setenv F77 "f77 -64"    # FORTRAN compiler
        endif
	setenv CFLAGS "-O2 -DH5_USE_16_API  -xansi -D_POSIX_SOURCE"	# default C flags (optimize, ansi)
	setenv C_CFH ""	 		# C w/ cfortran.h callable from FORTRAN
	setenv CFHFLAGS "$CFLAGS $C_CFH" # CFLAGS + C_CFH
	setenv F77FLAGS "-O2 -DH5_USE_16_API " 		# common FORTRAN flags
        setenv F77_CFH ""               # FORTRAN callable from C w/ cfortran.h
	setenv F77_C_CFH ""          	# calling C w/ cfortran.h
	setenv F77_C_LIB ""		# C lib called by FORTRAN main
	setenv HDFSYS IRIX		# system type as defined by HDF
	breaksw

    case irix65:
	setenv CC "cc -n32"		# C compiler (new-style 32 bit)
	setenv F77 "f77 -n32"		# FORTRAN compiler (new-style 32 bit)
	setenv CFLAGS "-O2 -DH5_USE_16_API  -xansi -D_POSIX_SOURCE"	# default C flags (optimize, ansi)
	setenv C_CFH ""	 		# C w/ cfortran.h callable from FORTRAN
	setenv CFHFLAGS "$CFLAGS $C_CFH" # CFLAGS + C_CFH
	setenv F77FLAGS "-O2 -DH5_USE_16_API " 		# common FORTRAN flags
        setenv F77_CFH ""               # FORTRAN callable from C w/ cfortran.h
	setenv F77_C_CFH ""          	# calling C w/ cfortran.h
	setenv F77_C_LIB ""		# C lib called by FORTRAN main
	setenv HDFSYS IRIX		# system type as defined by HDF
	breaksw

	
    case sun4:
	setenv CC acc			# C compiler
	setenv CFLAGS "-O2 -DH5_USE_16_API  -Xa" 		# default C flags (optimize, ansi)
	setenv C_CFH "-DsunFortran"	# C w/ cfortran.h callable from FORTRAN
	setenv CFHFLAGS "$CFLAGS $C_CFH" # CFLAGS + C_CFH
	setenv F77 f77 			# FORTRAN compiler
	setenv F77FLAGS "-O2 -DH5_USE_16_API " 		# common FORTRAN flags
        setenv F77_CFH ""               # FORTRAN callable from C w/ cfortran.h
	setenv F77_C_CFH ""           	# calling C w/ cfortran.h
	setenv F77_C_LIB "-lm" 		# C lib called by FORTRAN main
	setenv HDFSYS SUN		# system type as defined by HDF
	breaksw

    case sun5:
	setenv CC cc			# C compiler
	setenv CFLAGS "-O2 -DH5_USE_16_API  -Xa" 		# default C flags (optimize, ansi)
	setenv C_CFH "-DsunFortran"	# C w/ cfortran.h callable from FORTRAN
	setenv CFHFLAGS "$CFLAGS $C_CFH" # CFLAGS + C_CFH
	setenv F77 f77 			# FORTRAN compiler
	setenv F77FLAGS "-O2 -DH5_USE_16_API " 		# common FORTRAN flags
        setenv F77_CFH ""               # FORTRAN callable from C w/ cfortran.h
	setenv F77_C_CFH ""           	# calling C w/ cfortran.h
	setenv F77_C_LIB "-lm" 		# C lib called by FORTRAN main
	setenv HDFSYS SUN		# system type as defined by HDF
	breaksw

     case sun5.8:
	setenv CC cc			# C compiler
	setenv CFLAGS "-O2 -DH5_USE_16_API  -Xa" 		# default C flags (optimize, ansi)
	setenv C_CFH "-DsunFortran -DSUNOS58X"	# C w/ cfortran.h callable from FORTRAN
	setenv CFHFLAGS "$CFLAGS $C_CFH" # CFLAGS + C_CFH
	setenv F77 f77 			# FORTRAN compiler
	setenv F77FLAGS "-O2 -DH5_USE_16_API " 		# common FORTRAN flags
        setenv F77_CFH ""               # FORTRAN callable from C w/ cfortran.h
	setenv F77_C_CFH ""           	# calling C w/ cfortran.h
	setenv F77_C_LIB "-lm" 		# C lib called by FORTRAN main
	setenv HDFSYS SUN		# system type as defined by HDF
	breaksw

     case sun5.9:
	setenv CC cc			# C compiler
	setenv CFLAGS "-O2 -DH5_USE_16_API  -Xa" 		# default C flags (optimize, ansi)
	setenv C_CFH "-DsunFortran -DSUNOS59X"	# C w/ cfortran.h callable from FORTRAN
	setenv CFHFLAGS "$CFLAGS $C_CFH" # CFLAGS + C_CFH
	setenv F77 f77 			# FORTRAN compiler
	setenv F77FLAGS "-O2 -DH5_USE_16_API " 		# common FORTRAN flags
        setenv F77_CFH ""               # FORTRAN callable from C w/ cfortran.h
	setenv F77_C_CFH ""           	# calling C w/ cfortran.h
	setenv F77_C_LIB "-lm" 		# C lib called by FORTRAN main
	setenv HDFSYS SUN		# system type as defined by HDF
	breaksw

     case sun5.10:
	setenv CC cc			# C compiler
	setenv CFLAGS "-O2 -DH5_USE_16_API  -Xa" 		# default C flags (optimize, ansi)
	setenv C_CFH "-DsunFortran -DSUNOS510X"	# C w/ cfortran.h callable from FORTRAN
	setenv CFHFLAGS "$CFLAGS $C_CFH" # CFLAGS + C_CFH
	setenv F77 f77 			# FORTRAN compiler
	setenv F77FLAGS "-O2 -DH5_USE_16_API " 		# common FORTRAN flags
        setenv F77_CFH ""               # FORTRAN callable from C w/ cfortran.h
	setenv F77_C_CFH ""           	# calling C w/ cfortran.h
	setenv F77_C_LIB "-lm" 		# C lib called by FORTRAN main
	setenv HDFSYS SUN		# system type as defined by HDF
	breaksw

    default:
	setenv CC cc			# C compiler
	setenv CFLAGS "-O2 -DH5_USE_16_API " 		# default C flags (optimize)
	setenv C_CFH ""	        	# C w/ cfortran.h callable from FORTRAN
	setenv CFHFLAGS "$CFLAGS $C_CFH" # CFLAGS + C_CFH
	setenv F77 f77 			# FORTRAN compiler
	setenv F77FLAGS "-O2 -DH5_USE_16_API " 		# common FORTRAN flags
        setenv F77_CFH ""               # FORTRAN callable from C w/ cfortran.h
	setenv F77_C_CFH ""           	# calling C w/ cfortran.h
	setenv F77_C_LIB "-lm" 		# C lib called by FORTRAN main
	setenv HDFSYS unknown		# system type as defined by HDF
	breaksw
endsw


# 
# set up environment to handle FORTRAN-90 compiler
#

if ("$pgs_f90_comp" != "") then		# using FORTRAN-90

    setenv F77 "$pgs_f90_comp"

    if ("$pgs_nag_flag" == "1") then	# using NAG f90
        setenv C_CFH "$C_CFH -DNAGf90F"
        setenv CFHFLAGS "$CFLAGS $C_CFH"
    endif

endif


if ("$pgs_absoft_flag" == "1") then	# using ABSOFT f77
     setenv F77 "$fc_path/f77 "
     setenv F77FLAGS "-W -f -s -N15 -N26 OPT_LVL"
endif

# copy the machine-specific path to variable pgs_path

set pgs_path = ($path)

# set PGS-related environment variables
# these may be referred to in makefiles and on compiler command lines

if ( $?PGSHOME ) then

# set up base set of PGS Toolkit directory variables.

    setenv PGSBIN 	$PGSHOME/bin/$BRAND		# exectuable files
    setenv PGSDAT 	$PGSHOME/database/$BRAND	# database files
    setenv PGSINC 	$PGSHOME/include		# include (header) files
    setenv PGSLIB 	$PGSHOME/lib/$BRAND		# library files
    setenv PGSMSG 	$PGSHOME/message		# SMF message files
    setenv PGSOBJ 	$PGSHOME/obj/$BRAND		# object files
    setenv PGSRUN 	$PGSHOME/runtime		# runtime work files
    setenv PGSSRC 	$PGSHOME/src			# toolkit source files
    setenv PGSTST 	$PGSHOME/test			# test source files


    if ("$use_flavor" == 1) then
	if ("$use_flavor_r" == 1) then
	    setenv PGSBIN 	${PGSBIN}_${flavor_r}	# exectuable files
	    setenv PGSLIB 	${PGSLIB}_${flavor_r} 	# library files
	    setenv PGSOBJ 	${PGSOBJ}_${flavor_r}	# object files
	else 
	    setenv PGSBIN 	${PGSBIN}_${flavor}	# exectuable files
	    setenv PGSLIB 	${PGSLIB}_${flavor} 	# library files
	    setenv PGSOBJ 	${PGSOBJ}_${flavor}	# object files
	endif
    endif

# update path variables

    set path  = ($path $PGSBIN)			# add PGSBIN to pgs  path
    set pgs_path  = ($pgs_path $PGSBIN)		# add PGSBIN to pgs  path
    set user_path = ($user_path $PGSBIN)	# add PGSBIN to user path

# set up variables needed by Process Control (PC) tools. 
# user may customize these in their private .cshrc file

    setenv PGS_PC_INFO_FILE    	$PGSRUN/$BRAND/PCF.relB0 
    if ( $use_flavor == 1 ) then
	if ( $flavor =~ *daac* ) then
	    setenv PGS_PC_INFO_FILE    	$PGSRUN/$BRAND/PCF.relB0_daac
	else
	    setenv PGS_PC_INFO_FILE    	$PGSRUN/$BRAND/PCF.relB0_scf
	endif
    endif



else

    echo "You must first set the environment variable PGSHOME"

endif


# set HDF-related environment variables
# these may be referred to in makefiles and on compiler command lines
# use the command 'sethdf <hdf-home-directory> to override the default

####
#
#
# if your site has HDF installed:
#    uncomment the next line and edit it to properly 
#    set your HDF home directory
#
####

# we will define a compilation DFLAG if hdf4 netCDF is disabled, which
# means "sd_" prefix is added to the netCDF function names

setenv hdf4_netcdf enabled

setenv HDFHOME /home/jaemoo/seadas-7.4/ocssw/build/lib3	# set HDF home directory
setenv HDF5HOME /home/jaemoo/seadas-7.4/ocssw/build/lib3        # set HDF5 home directory

if ( $use_flavor == 1 ) then
    if ( $flavor =~ *debug ) then
	set hdfhome=`echo $HDFHOME | sed "s/${BRAND}/${BRAND}_debug/"`
        set hdf5home=`echo $HDF5HOME | sed "s/${BRAND}/${BRAND}_debug/"`
        if ( -d $hdfhome ) then
	    setenv HDFHOME $hdfhome
        endif
	unset hdfhome
        if ( -d $hdf5home ) then
            setenv HDF5HOME $hdf5home
        endif
        unset hdf5home
    endif
endif

#
# Set HDF directories:
#
# - first look in $HDFHOME
# - if not found, default to $HDFHOME/hdf
#

if ($?HDFHOME) then

    if ( -d $HDFHOME/bin ) then		# HDF utilities
        setenv HDFBIN $HDFHOME/bin
    else
        setenv HDFBIN $HDFHOME/hdf/bin
    endif

    if ( -d $HDFHOME/include ) then	# HDF header files
        setenv HDFINC $HDFHOME/include
    else
        setenv HDFINC $HDFHOME/hdf/include
    endif

    if ( -d $HDFHOME/lib ) then		# HDF libraries
        setenv HDFLIB $HDFHOME/lib
    else
        setenv HDFLIB $HDFHOME/hdf/lib
    endif



    set pgs_path  = ($pgs_path  $HDFBIN)	# add HDFBIN to pgs  path
    set user_path = ($user_path $HDFBIN)	# add HDFBIN to user path

endif
 
# check to see if HDF4 is installed with --disable-netcdf
# if that is the case then netcdf funtions in libmfhdf.a will hav "sd_" prefix
# set a a TOOLKIT compilation flag to call netCDF function with or without
# the "sd" prefix 

    if (-f $HDFLIB/libmfhdf.a) then
	set check_netcdf = `strings $HDFLIB/libmfhdf.a|grep sd_ncopen`
	if ("$check_netcdf" != "") then
	    setenv hdf4_netcdf disabled
	else
	    setenv hdf4_netcdf enabled
	endif
   endif

setenv ZLIBHOME /home/jaemoo/seadas-7.4/ocssw/build/lib3  # set ZLIB home directory

#
# Set ZLIB directories:
#
# - first look in $ZLIBHOME
# - if not found, default to $ZLIBHOME/zlib
#
 
if ($?ZLIBHOME) then
 
if ( -d $ZLIBHOME/bin ) then         # ZLIB utilities
        setenv ZLIBBIN $ZLIBHOME/bin
    else
        setenv ZLIBBIN $ZLIBHOME/zlib/bin
    endif
 
if ( -d $ZLIBHOME/include ) then     # ZLIB header files
        setenv ZLIBINC $ZLIBHOME/include
    else
        setenv ZLIBINC $ZLIBHOME/zlib/include
    endif
 
    if ( -d $ZLIBHOME/lib ) then         # ZLIB libraries
        setenv ZLIBLIB $ZLIBHOME/lib
    else
        setenv ZLIBLIB $ZLIBHOME/zlib/lib
    endif
 
 
    set pgs_path  = ($pgs_path  $ZLIBBIN)        # add ZLIBBIN to pgs  path
    set user_path = ($user_path $ZLIBBIN)        # add ZLIBBIN to user path
 
endif

setenv JPEGHOME /home/jaemoo/seadas-7.4/ocssw/build/lib3  # set JPEG home directory

#
# Set JPEG directories:
#
# - first look in $JPEGHOME
# - if not found, default to $JPEGHOME/jpeg
#
 
if ($?JPEGHOME) then
 
if ( -d $JPEGHOME/bin ) then         # JPEG  utilities
        setenv JPEGBIN $JPEGHOME/bin
    else
        setenv JPEGBIN $JPEGHOME/jpeg/bin
    endif
 
if ( -d $JPEGHOME/include ) then     # JPEG header files
        setenv JPEGINC $JPEGHOME/include
    else
        setenv JPEGINC $JPEGHOME/jpeg/include
    endif
 
    if ( -d $JPEGHOME/lib ) then         # JPEG libraries
        setenv JPEGLIB $JPEGHOME/lib
    else
        setenv JPEGLIB $JPEGHOME/jpeg/lib
    endif
 
 
    set pgs_path  = ($pgs_path  $JPEGBIN)        # add JPEGBIN to pgs  path
    set user_path = ($user_path $JPEGBIN)        # add JPEGBIN to user path
 
endif

#
# Set HDF5 directories:
#
# - first look in $HDF5HOME
# - if not found, default to $HDF5HOME/hdf5
#

if ($?HDF5HOME) then

    if ( -d $HDF5HOME/bin ) then                # HDF5 utilities
        setenv HDF5BIN $HDF5HOME/bin
    else
        setenv HDF5BIN $HDF5HOME/hdf5/bin
    endif

    if ( -d $HDF5HOME/include ) then    	# HDF5 header files
        setenv HDF5INC $HDF5HOME/include
    else
        setenv HDF5INC $HDF5HOME/hdf5/include
    endif

    if ( -d $HDF5HOME/lib ) then                # HDF5 libraries
        setenv HDF5LIB $HDF5HOME/lib
    else
        setenv HDF5LIB $HDF5HOME/hdf5/lib
    endif



    set pgs_path  = ($pgs_path  $HDF5BIN)       # add HDF5BIN to pgs  path
    set user_path = ($user_path $HDF5BIN)       # add HDF5BIN to user path

endif

setenv SZIPHOME /home/jaemoo/seadas-7.4/ocssw/build/lib3  # set SZIP home directory
#
# Set SZIP directories:
#
# - first look in $SZIPHOME
# - if not found, default to $SZIPHOME/szip
#
 
if ($?SZIPHOME) then
 
if ( -d $SZIPHOME/bin ) then         # SZIP utilities
        setenv SZIPBIN $SZIPHOME/bin
    else
        setenv SZIPBIN $SZIPHOME/szip/bin
    endif
 
if ( -d $SZIPHOME/include ) then     # SZIP header files
        setenv SZIPINC $SZIPHOME/include
    else
        setenv SZIPINC $SZIPHOME/szip/include
    endif
 
    if ( -d $SZIPHOME/lib ) then         # SZIP libraries
        setenv SZIPLIB $SZIPHOME/lib
    else
        setenv SZIPLIB $SZIPHOME/szip/lib
    endif
 
 
    set pgs_path  = ($pgs_path  $SZIPBIN)        # add SZIPBIN to pgs  path
    set user_path = ($user_path $SZIPBIN)        # add SZIPBIN to user path
 
endif
 

####
#
#
# if your site has HDF-EOS installed:
#    uncomment the next line and edit it to properly 
#    set your HDF-EOS home directory
#
####

setenv HDFEOS_HOME /home/jaemoo/seadas-7.4/ocssw/build/lib3/EOS	# set HDF-EOS home directory

#
# Set HDF-EOS directories:
#

if ($?HDFEOS_HOME) then

    setenv HDFEOS_INC $HDFEOS_HOME/include
    setenv HDFEOS_LIB $HDFEOS_HOME/lib/$BRAND
    setenv HDFEOS_BIN $HDFEOS_HOME/bin/$BRAND

    if ( $use_flavor == 1 ) then
	if ( $flavor =~ *debug ) then
	    if ( -d ${HDFEOS_LIB}_debug ) then
		setenv HDFEOS_LIB ${HDFEOS_LIB}_debug
		setenv HDFEOS_BIN ${HDFEOS_BIN}_debug
	    endif
	endif
    endif

endif

####
#
#
# if your site has HDF-EOS5 installed:
#    uncomment the next line and edit it to properly 
#    set your HDF-EOS home directory
#
####

setenv HDFEOS5_HOME /home/jaemoo/seadas-7.4/ocssw/build/lib3/EOS        # set HDF-EOS5 home directory

#
# Set HDF-EOS5 directories:
#

if ($?HDFEOS5_HOME) then

    setenv HDFEOS5_INC $HDFEOS5_HOME/include
    setenv HDFEOS5_LIB $HDFEOS5_HOME/lib/$BRAND
    setenv HDFEOS5_BIN $HDFEOS5_HOME/bin/$BRAND
    if ( $use_flavor == 1 ) then
        if ( $flavor =~ *debug ) then
            if ( -d ${HDFEOS5_LIB}_debug ) then
                setenv HDFEOS5_LIB ${HDFEOS5_LIB}_debug
                setenv HDFEOS5_BIN ${HDFEOS5_BIN}_debug
            endif
        endif
    endif

endif

#
# check to see if DAAC version of the toolkit was installed
#

if ( $?ADD_IFLAGS ) setenv ADD_IFLAGS "" # make sure old values are gone
if ( $?ADD_LFLAGS ) setenv ADD_LFLAGS "" # make sure old values are gone
if ( $?ADD_LIBS )   setenv ADD_LIBS ""   # make sure old values are gone

if ("$pgs_daac" == "1") then

    setenv CFLAGS "$CFLAGS -DPGS_DAAC"
    setenv CFHFLAGS "$CFHFLAGS -DPGS_DAAC"

endif

if ( $?ECS_VERSION ) then

    setenv VFLAG "'-DPGSd_ECS_VERSION="'"${ECS_VERSION}"'"'"

endif

if ("$pgs_c_rlib" == "1") then

   switch ($BRAND)
	case sun5:
    	   setenv CFLAGS "$CFLAGS -D_PGS_THREADSAFE -D_POSIX_PTHREAD_SEMANTICS -D_REENTRANT"
           setenv CFHFLAGS "$CFHFLAGS -D_PGS_THREADSAFE -D_POSIX_PTHREAD_SEMANTICS -D_REENTRANT"
           setenv PGSRLIB "YES"		# Threadsafe library files
	   breaksw
       	case sun5.8:
    	   setenv CFLAGS "$CFLAGS -D_PGS_THREADSAFE -D_POSIX_PTHREAD_SEMANTICS -D_REENTRANT"
           setenv CFHFLAGS "$CFHFLAGS -D_PGS_THREADSAFE -D_POSIX_PTHREAD_SEMANTICS -D_REENTRANT"
           setenv PGSRLIB "YES"		# Threadsafe library files
	   breaksw
       	case sun5.9:
    	   setenv CFLAGS "$CFLAGS -D_PGS_THREADSAFE -D_POSIX_PTHREAD_SEMANTICS -D_REENTRANT"
           setenv CFHFLAGS "$CFHFLAGS -D_PGS_THREADSAFE -D_POSIX_PTHREAD_SEMANTICS -D_REENTRANT"
           setenv PGSRLIB "YES"		# Threadsafe library files
	   breaksw
       	case sun5.10:
    	   setenv CFLAGS "$CFLAGS -D_PGS_THREADSAFE -D_POSIX_PTHREAD_SEMANTICS -D_REENTRANT"
           setenv CFHFLAGS "$CFHFLAGS -D_PGS_THREADSAFE -D_POSIX_PTHREAD_SEMANTICS -D_REENTRANT"
           setenv PGSRLIB "YES"		# Threadsafe library files
	   breaksw
	case sgi:
     	   setenv CFLAGS "$CFLAGS"
           setenv CFHFLAGS "$CFHFLAGS"
	   setenv  PGSRLIB "YES"
	breaksw
	case sgi32:
     	   setenv CFLAGS "$CFLAGS -D_PGS_THREADSAFE -D_REENTRANT -D_POSIX_C_SOURCE=199506L"
           setenv CFHFLAGS "CFHFLAGS -D_PGS_THREADSAFE -D_REENTRANT -D_POSIX_C_SOURCE=199506L"
	   setenv  PGSRLIB "YES"
	breaksw
	case sgi64:
     	   setenv CFLAGS "$CFLAGS -D_PGS_THREADSAFE -D_REENTRANT -D_POSIX_C_SOURCE=199506L"
           setenv CFHFLAGS "$CFHFLAGS -D_PGS_THREADSAFE -D_REENTRANT -D_POSIX_C_SOURCE=199506L"
	   setenv  PGSRLIB "YES"
	breaksw
	default:
	  echo ""
	  breaksw
   endsw
endif

if ("$hdf4_netcdf" == "disabled") then
    setenv CFLAGS "$CFLAGS -DHDF4_NETCDF_HAVE_SD"
    setenv CFHFLAGS "$CFHFLAGS -DHDF4_NETCDF_HAVE_SD"
endif

#
# restore augmented user path
#
set path = ( $user_path )


# done
