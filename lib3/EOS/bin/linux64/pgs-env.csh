# 
# set the toolkit home directory environment variable
# SDP installation done on Sun Jan 21 18:28:27 KST 2018 
# 
 
setenv PGSHOME /home/jaemoo/seadas-7.4/ocssw/build/lib3/src/sdptk/TOOLKIT	# set the toolkit home directory
set pgs_formal=1	# formal directory structure (1=yes/0=no)
set sgi_mode=32			# set SGI for standard mode
set LINUX_BRAND=linux64		# set LINUX for -64 mode
setenv AA_install 0    #set option for AA tool instalation
 
# file: pgs-env.csh (C-shell version)
# set up the PGS directory environment variables
# this version is for Release B0 of the toolkit, DAAC and SCF versions

if (! $?sgi_mode ) then		# by default, SGI mode is standard 32-bit
    set sgi_mode=32
endif

if (! $?use_flavor ) then 	# by default, do not use "flavors"
    set use_flavor=0
endif

# save user path, then set path to a base subset of directories for unknown host

set user_path = ( $path )	# save user path
set path = (/usr/local/bin /bin /usr/bin /etc /usr/etc /usr/ucb /usr/bin/X11)

# get operating system type, login name

switch ( `uname -m | awk '{print $1}'` )
  case "i386": 
    setenv OSTYPE sco386 ; breaksw # SCO box
  case "CRAY": 
    setenv OSTYPE UNICOS ; breaksw # CRAY
  default:    
    setenv OSTYPE `uname` ; breaksw # everybody else
endsw

setenv CYGPL `uname | awk -F_ '{print $1}'`

    # Intel Macintosh is also i386 or i686 (?) machine

    if("`uname -m | awk '{print $1}'`" == "i386") then
	if("$OSTYPE" == "Darwin") then 
	    setenv OSTYPE DarwinIntel
	endif
	if("$CYGPL" == "CYGWIN") then 
	    setenv OSTYPE Cygwin
	endif
    endif
    if("`uname -m | awk '{print $1}'`" == "i686") then
	if("$OSTYPE" == "Darwin") then 
	    setenv OSTYPE DarwinIntel
	endif
	if("$CYGPL" == "CYGWIN") then 
	    setenv OSTYPE Cygwin
	endif
    endif
    if( "`uname -m | awk '{print $1}'`" == "x86_64" || "`uname -m | awk '{print $1}'`" == "ia64" ) then
	if("$OSTYPE" == "Darwin") then 
	    setenv OSTYPE DarwinIntel
	endif
	if("$CYGPL" == "CYGWIN") then 
	    setenv OSTYPE Cygwin
	endif
     endif
# set variable BRAND,  used in achitecture-specific path names

switch ( "$OSTYPE" )
    case "AIX": 
	setenv BRAND ibm ; breaksw 
    case "HP-UX": 
	setenv BRAND hp  ; breaksw 
    case "IRIX":  
	setenv BRAND sgi  ; breaksw 
    case "IRIX64":  
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
    case "OSF1":  
	setenv BRAND dec  ; breaksw 
    case "sco386": 
	setenv BRAND sco  ; breaksw 
    case "SunOS": 
	# distinguish between SunOS 5.x and 4.x versions
	set cbrand=`uname -r | awk -F. '{print $1, $2}'`
	if ("$cbrand" == "5 10") then 
	    setenv BRAND sun5.10	
        else if ("$cbrand" == "5 9") then 
	    setenv BRAND sun5.9	
        else if ("$cbrand" == "5 8") then 
	    setenv BRAND sun5.8	
        else if ("$cbrand" == "5 5") then
            setenv BRAND sun5        
         else
            setenv BRAND sun4
         endif
        breaksw
    case "UNICOS": 
	setenv BRAND cray ; breaksw
    case "Linux":
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
        setenv BRAND macintosh ; breaksw
    case "DarwinIntel":
        setenv BRAND macintel32
	if("$MAC_BRAND" == "macintel32") then 
		setenv BRAND macintel32
	else if("$MAC_BRAND" == "macintel") then 
		setenv BRAND macintel
	else if("$MAC_BRAND" == "macintel64") then
		setenv BRAND macintel64
	else if( "`uname -m | awk '{print $1}'`" == "x86_64" ||  "`uname -m | awk '{print $1}'`" == "ia64"   || "`uname -m | awk '{print $1}'`" == "i386"   || "`uname -m | awk '{print $1}'`" == "i686" ) then
		setenv BRAND macintel64
	else
		setenv BRAND macintel32
	endif
	breaksw
    case "Cygwin":
        setenv BRAND cygwin ; breaksw
    default:
	echo "Operating system: $OSTYPE not supported"
        breaksw
endsw


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


    if ( $use_flavor == 1 ) then

        setenv PGSBIN 	${PGSBIN}_${flavor}		# exectuable files
        setenv PGSLIB 	${PGSLIB}_${flavor} 		# library files
        setenv PGSOBJ 	${PGSOBJ}_${flavor}		# object files

    endif

#    update path variables

    set path  = ($path $PGSBIN)			# add PGSBIN to pgs  path
    set pgs_path  = ($pgs_path $PGSBIN)		# add PGSBIN to pgs  path
    set user_path = ($user_path $PGSBIN)	# add PGSBIN to user path

#   Set up variables needed by Process Control (PC) tools. 
#   User may customize these in their private .cshrc file

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

# restore user path

set path = ( $user_path )

# install option:  -batch -fc_path /usr/bin/gfortran -cc_path /usr/bin/gcc -hdfhome /home/jaemoo/seadas-7.4/ocssw/build/lib3 -hdf5home /home/jaemoo/seadas-7.4/ocssw/build/lib3 -hdfeos_home /home/jaemoo/seadas-7.4/ocssw/build/lib3/EOS -hdfeos5_home /home/jaemoo/seadas-7.4/ocssw/build/lib3/EOS -zlibhome /home/jaemoo/seadas-7.4/ocssw/build/lib3 -jpeghome /home/jaemoo/seadas-7.4/ocssw/build/lib3 -sziphome /home/jaemoo/seadas-7.4/ocssw/build/lib3
