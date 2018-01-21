# 
# set the toolkit home directory environment variable
# SDP installation done on Sun Jan 21 18:28:17 KST 2018 
# 
 
PGSHOME=/home/jaemoo/seadas-7.4/ocssw/build/lib3/src/sdptk/TOOLKIT	# set the toolkit home directory
pgs_formal=1		# formal directory structure (1=yes/0=no)
sgi_mode=32			# set SGI for standard mode
LINUX_BRAND=linux64		# set LINUX for -64 mode
AA_install=0    #set option for AA tool instalation
export AA_install
 
# file: pgs-env.ksh (Korn/Bourne shell version)
# set up the PGS directory environment variables
# this version is for Release A of the toolkit, DAAC and SCF versions


: ${sgi_mode:=32}	# by default, SGI mode is standard 32-bit

: ${use_flavor:=0}	# by default, do not use "flavors"


# save user path, then set path to a base subset of directories for unknown host

user_path=$PATH
PATH=/usr/local/bin:/bin:/usr/bin:/etc:/usr/etc:/usr/ucb:/usr/bin/X11
export PATH 

# get operating system type

case "`uname -m | awk '{print $1}'`" in
  i386) 	OSTYPE=sco386   	;; 	# SCO box
  CRAY) 	OSTYPE=UNICOS   	;; 	# CRAY
  *) 		OSTYPE="`uname`" 	;; 	# everybody else
esac

CYGPL="`uname | awk -F_ '{print $1}'`"  
   
    # Intel Macintosh is also i386 or i686 (?) machine

    if [ "`uname -m | awk '{print $1}'`" = "i386" ] ; then
	if [ "`uname`" = "Darwin" ] ; then 
	    OSTYPE=DarwinIntel
	fi
	if [ "$CYGPL" = "CYGWIN" ] ; then 
	    OSTYPE=Cygwin
	fi
    fi
    if [ "`uname -m | awk '{print $1}'`" = "i686" ] ; then
	if [ "`uname`" = "Darwin" ] ; then 
	    OSTYPE=DarwinIntel
	fi
	if [ "$CYGPL" = "CYGWIN" ] ; then 
	    OSTYPE=Cygwin
	fi
    fi
    if [ "`uname -m | awk '{print $1}'`" = "x86_64" ] || [ "`uname -m | awk '{print $1}'`" = "ia64" ] ; then
	if [ "$temp_ostype" = "Darwin" ] ; then 
	    	OSTYPE=DarwinIntel
		export OSTYPE
	    	if [ "$MAC_BRAND" = "" ] ; then 
			echo " Error: In 64-bit MAC platform the environment variable MAC_BRAND must be set to macintel32 or macintel64 before running this script."
	    	fi
	else	
		if [ "$LINUX_BRAND" = "" ] ; then 
	    		echo " Error: In 64-bit linux platform the environment variable LINUX_BRAND must be set to linux32 or linux64 before running this script."
		fi
	fi
    fi

# set variable BRAND,  used in achitecture-specific path names

case "$OSTYPE" in
  AIX) 	BRAND="ibm" ;;
  HP-UX) 	BRAND="hp" ;;
  IRIX) 	BRAND="sgi" ;;
  IRIX64) 
    case $sgi_mode in
      64 ) BRAND=sgi64 ;;
      n32) BRAND=sgi32 ;;
      32 ) BRAND=sgi ;;
      65 ) BRAND=irix65 ;;
      *  ) BRAND=sgi ;;  # just in case
    esac
    ;;
  OSF1) 	BRAND="dec" ;;
  sco386) 	BRAND="sco" ;;
  SunOS) 	# distinguish between SunOS 5.x and 4.x versions
    if [ "`uname -r | awk -F. '{print $1, $2}'`" = "5 10" ] ; then 
        BRAND="sun5.10"
   elif [ "`uname -r | awk -F. '{print $1, $2}'`" = "5 9" ] ; then 
        BRAND="sun5.9"
   elif [ "`uname -r | awk -F. '{print $1, $2}'`" = "5 8" ] ; then 
        BRAND="sun5.8"
   elif
        [ "`uname -r | awk -F. '{print $1}'`" = "4" ] ; then 
        BRAND="sun4"
   else
	BRAND="sun5"
   fi	
    ;;
  UNICOS) BRAND="cray" ;;
  Linux )
	BRAND=linux
	if [ "$LINUX_BRAND" = "linux64" ] ; then
	    BRAND=linux64
	elif [ "$LINUX_BRAND" = "linux32" ] ; then
	    BRAND=linux32
	elif [ "$LINUX_BRAND" = "linux" ] ; then
	    BRAND=linux
	fi
	;;
  Darwin )
        BRAND=macintosh
        ;;
  DarwinIntel )
	BRAND=macintel32
	if [ "$MAC_BRAND" = "macinel64" ] ; then
		BRAND=macintel64
    	elif [ "$MAC_BRAND" = "macintel32" ] ; then
		BRAND=macintel32
    	fi
        ;;
  Cygwin)
        BRAND=cygwin
        ;;
  *)	 echo "Operating system: $OSTYPE not supported" ;;
esac
export BRAND


# copy the machine-specific path to variable pgs_path

pgs_path=$PATH

# set PGS-related environment variables
# these may be referred to in makefiles and on compiler command lines

if [ "$PGSHOME" != "" ] ; then

# set up base set of PGS Toolkit directory variables.

    PGSBIN=$PGSHOME/bin/$BRAND		# executable files
    PGSDAT=$PGSHOME/database/$BRAND	# database files
    PGSINC=$PGSHOME/include		# include header files
    PGSLIB=$PGSHOME/lib/$BRAND 		# library files
    PGSMSG=$PGSHOME/message		# SMF message files
    PGSOBJ=$PGSHOME/obj/$BRAND		# object files
    PGSRUN=$PGSHOME/runtime		# runtime work files
    PGSSRC=$PGSHOME/src			# toolkit source files
    PGSTST=$PGSHOME/test		# test source files

    if [ $use_flavor = 1 ] ; then

        PGSBIN=${PGSBIN}_${flavor}		# executable files
        PGSLIB=${PGSLIB}_${flavor}  		# library files
        PGSOBJ=${PGSOBJ}_${flavor}		# object files

    fi

    export PGSHOME PGSBIN PGSDAT PGSINC PGSLIB 
    export PGSMSG  PGSOBJ PGSRUN PGSSRC PGSTST

# update path variables

    PATH=$PATH:$PGSBIN; export PATH		# add PGSBIN to path
    pgs_path=$pgs_path:$PGSBIN			# add PGSBIN to pgs  path
    user_path=$user_path:$PGSBIN		# add PGSBIN to user path

# set up variables needed by Process Control (PC) tools. 
# user may customize these in their private .cshrc file

    PGS_PC_INFO_FILE=$PGSRUN/$BRAND/PCF.relB0
    if [ $use_flavor = 1 ] ; then
	case $flavor in
	    *daac*)
		PGS_PC_INFO_FILE=$PGSRUN/$BRAND/PCF.relB0_daac
	    ;;
	    *scf*)
		PGS_PC_INFO_FILE=$PGSRUN/$BRAND/PCF.relB0_scf
	    ;;
	esac
    fi

    export PGS_PC_INFO_FILE 



else

    echo "You must first set the environment variable PGSHOME"

fi


# restore user path

PATH=$user_path
export PATH 

