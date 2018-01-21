
#-------------------------------------------------------------------------#
#                                                                         #
#  COPYRIGHT[copyright mark] 2000, Raytheon System Company, its vendors,  #
#  and suppliers.  ALL RIGHTS RESERVED.                                   #
#                                                                         #
#-------------------------------------------------------------------------#
#------------------------------------------------------------------------------
#BEGIN_FILE_PROLOG:
#
#FILENAME:
#	PGS_PC_Shell.sh
#
#DESCRIPTION:
#	This file contains all script functions neccessary
#	to run the PGS_PC_Shell.sh command.
#
#AUTHOR:
#	Ray Milburn / Applied Research Corp.
#
#HISTORY:
#	06-Dec-94 RM Initial Version
#	28-Dec-94 RM Fixing items found during code inspection held on
#			22-Dec-94.
#	17-Feb-95 RM Updating prolog.
#	28-Mar-95 RM Adding third argument to accept location of PCF.
#			This is a modification for TK5.
#	25-Apr-95 RM Added support for SMF message cache size parameter.
#			Added status message.  Revised -h option code.
#			Added -v option for verbose output.
#			All changes were copied from Mike Sucher.
#	13-Jul-95 RM Added checks to determine if PCF and PGE are regular
#			files and contain the permission necessary to run
#			this shell.  This was mandated by DR ECSed00997.
#	21-Jul-95 RM Added comments to be placed in User's Guide regarding
#			FORTRAN limitations in shared memory per DR
#			ECSed00949.
#	19-Oct-95 DH Added support for MSS Event Log Setup Error.
#	10-Jan-96 RM Added support for return values from the PGE.  This 
#			was accomplished by expanding the -v option and 
#			to install the -p option.
#	26-Mar-96 RM Fixed DR ECSed01826.  Removed square brackets from
#			SMF Cache Size command line variable listed in
#			usage function.
#	14-Aug-96 RM Fixed DR ECSed02851.  Added internal documentation 
#			in NOTES and RETURNS sections informing user
#			that return values may be different based on 
#			which platform and shell is being used.
#	20-Aug-96 RM DR ECSed00949 fix has been implemented.  The 
#			documentation in the NOTES section about using	
#			shared memory with FORTRAN is removed.
#
#END_FILE_PROLOG:
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
#BEGIN_PROLOG:
#
#TITLE:
#	Shell script to run the PGE.
#
#NAME:
#	PGS_PC_Shell.sh
#
#SYNOPSIS:
#	PGS_PC_Shell.sh [-h] <PGE file> <Init string> <PCF location> 
#			<SMF Cache Size> [-v] [-p]
#
#DESCRIPTION:
#	This shell script accepts four command line arguments as input.
#	The first argument is the PGE to run.  This may be a shell script
#	or an executable.  The second argument is the Init string which
#	contains 4 binary digits that define how the Tollkit will behave.  
#	Together, these instruct the shell about what to do in the case of 
#	using/not using shared memory or using/not using log files.  The
#	third argument is the location of the Process Control File (PCF).
#	The forth argument is the SMF cache size.  A fifth argument may
#	be used.  This argument will run this script in verbose mode.
#
#INPUTS:
#	Name		Description			Units	Min	Max
#
#	PGE file	The full path/file name of the
#			PGE to be run.
#
#	Init string	The string to be passed in with
#			the instructions about what to
#			do with shared memory and the 
#			log file.  See NOTES section
#			for complete description of each
#			field in the Init string flag.
#
#	PCF location	The full path/file name of the
#			Process Control File (PCF).
#       SMF Cache Size  size of SMF message cache in records
#
#       -v              Run in verbose mode. Output status
#                       messages displaying settings, current
#                       file being run.
#
#	-p		Make the return value of this script
#			be the return value of the PGE if
#			the PGE is run.  If the PGE does not
#			get run then revert to the normal 
#			method of return values for this
#			shell.
#
#	-h		Upon receiving the -h flag a 
#			short description of the usage
#			of PGS_PC_Shell.sh will be provided
#			to the user and the command
#			will exit.  
#
#OUTPUTS:
#	Name		Description			Units	Min	Max
# 
#	NONE
#
#RETURNS:
#	See NOTES section for explanation of values.
#
#	PGS_S_SUCCESS return is 0
#	PGS_SH_SMF_MSSLOGFILE return is 242 or -14
#	PGS_SH_PC_TRUNC return is 243 or -13
#	PGS_SH_PC_TOOLERROR return is 244 or -12
#	PGS_SH_PC_NODATA return is 245 or -11
#	PGS_SH_SYS_PARAM return is 246 or -10
#	PGS_SH_MEM_INIT return is 247 or -9
#	PGS_SH_PC_DELETETMP return is 248 or -8
#	PGS_SH_SMF_SENDRUNTIME return is 249 or -7
#	PGS_SH_SMF_SENDLOGFILE return is 250 or -6
#	PGS_SH_MEM_TERM return is 251 or -5
#	PGS_SH_SMF_LOGFILE return is 252 or -4
#	PGS_SH_PC_LOADDATA return is 253 or -3
#	PGS_SH_PC_ENV return is 254 or -2
#	PGS_SH_SMF_SHMMEM return is 255 or -1
#
#EXAMPLES:
#
#	PGS_PC_Shell.sh -h
#	PGS_PC_Shell.sh /usr/PGE/somePGE 1111 /usr/PGE/data/PCF.current 50 -v
#	PGS_PC_Shell.sh /usr/home/PGE/runFile 1010 /home/PCFDATA/pcf.data 200
#	PGS_PC_Shell.sh /usr/PGEhome/runThis 0000 /home/Data/MY.pcf 150 -p
#
#NOTES:
#	This shell script parses the input to ensure correctness and will
#	report any input problems to the user.
#
#	This shell script acts as the outer most shell for the PGE.
#
#	The Init string flag consists of four (4) fields.  Each field contains
#	a single digit.  The digits should  be a one (1) or a zero (0).  
#	Therefore the Init String would appears as "1010" or "1111", etc.  
#	For ease of use PGS_PC_Shell.sh will interpret any non-zero digit as 
#	a one.  Therefore, 8020 would be interpreted as 1010, and 5500 would 
#	be interpreted as 1100, etc.  The field descriptions are listed as 
#	follows:
#		
#		FIELD 1 - 1 (or any non-zero digit) = Use shared memory if available
#		          0 = Do not use shared memory
#
#		FIELD 2 - 1 (or any non-zero digit) = If shared memory fails continue
#                                                     using ASCII files
#		          0 = If shared memory fails stop now
#
#		FIELD 3 - 1 (or any non-zero digit) = Use Log Files
#		          0 = Do not use Log Files
#
#		FIELD 4 - 1 (or any non-zero digit) = If Log Files fail continue 
#                                                     anyway
#		          0 = If Log Files fail stop now
#
#       The return value from PGS_PC_Shell.sh is affected by the shell and
#       system being run on.  On some shells and systems a 255 will be 
#       returned as -1, a 254 will be -2, etc.  The user may look at the 
#	RETURNS section to match the value with the mnemonic.
#
#REQUIREMENTS:
#	PGSTK-1312
#
#DETAILS:
#	NONE
#
#GLOBALS:
#	NONE
#
#FILES:
#	NONE
#
#FUNCTIONS_CALLED:
#	Usage		- print a usage message
#       Error           - print an error message
#	Message		- prints status of Log and Shm flags
#	PGS_PC_InitCom	- C program which initializes the PGS Toolkit
#	PGS_PC_TermCom	- C program which terminates the PGS Toolkit
#
#END_PROLOG:
#------------------------------------------------------------------------------
# Define return values
: ${PGS_FALSE=0}
: ${PGS_TRUE=1}
: ${PGSd_PC_ATTRIBUTE_LOCATION=1}
: ${PGSd_PC_ATTRIBUTE_STRING=2}
: ${PGSd_IO_Gen_NoEndurance=0}
: ${PGSd_IO_Gen_Endurance=1}
: ${PGS_S_SUCCESS=0}
# error return values
: ${PGS_SH_SMF_MSSLOGFILE=242}
: ${PGS_SH_PC_TRUNC=243}
: ${PGS_SH_PC_TOOLERROR=244}
: ${PGS_SH_PC_NODATA=245}
: ${PGS_SH_SYS_PARAM=246}
: ${PGS_SH_MEM_INIT=247}
: ${PGS_SH_PC_DELETETMP=248}
: ${PGS_SH_SMF_SENDRUNTIME=249}
: ${PGS_SH_SMF_SENDLOGFILE=250}
: ${PGS_SH_MEM_TERM=251}
: ${PGS_SH_SMF_LOGFILE=252}
: ${PGS_SH_PC_LOADDATA=253}
: ${PGS_SH_PC_ENV=254}
: ${PGS_SH_SMF_SHMMEM=255}
# twos compliments of above values
: ${PGS_SH_SMF_MSSLOGFILE2=-14}
: ${PGS_SH_MEM_INIT2=-9}
: ${PGS_SH_SMF_LOGFILE2=-4}
: ${PGS_SH_SMF_SHMMEM2=-1}

# run PGE (1) or not (0)
: ${PGS_EXECUTE_PGE=1}

# Usage function
Usage()
{
    cat << !
 
usage:
 
    PGS_PC_Shell.sh <PGE file> <Init string> <PCF location> <SMF Cache Size> [-v] [-p] [-h]
 
    PGE file            - full path/file name of PGE script or executable.
 
    Init string         - 4-digit string containing flags for initialization.
 
                FIELD 1 - 1 (or any non-zero digit) = Use shared memory if available
                          0 = Do not use shared memory
 
                FIELD 2 - 1 (or any non-zero digit) = If shared memory fails continue
                                                      using ASCII files
                          0 = If shared memory fails stop now
 
                FIELD 3 - 1 (or any non-zero digit) = Use Log Files
                          0 = Do not use Log Files
 
                FIELD 4 - 1 (or any non-zero digit) = If Log Files fail continue 
                                                      anyway
                          0 = If Log Files fail stop now

    PCF location	- full path/file name of PCF.
 
    SMF Cache Size      - size of SMF message cache in records
 
    -v                  - verbose mode: output status messages

    -p			- return value of PGE return (if run)
 
    -h                  - usage help.
 
!
 
}


# Error function
Error()
{
        echo ""
        echo "Error:  $*"  >&2
}


# Message function
Message()
{
    echo "`basename $0`:"
    echo ""
    echo "    Running PGE file: $PGEFILE"
    echo ""
 
    if [ $SHMON -ne 0 ]
    then
        echo "    Shared memory is ENABLED"
    else
        echo "    Shared memory is DISABLED"
    fi
 
    if [ $SHMCONTINUE -ne 0 ]
    then
        echo "    Shared memory continue is ENABLED"
    else
        echo "    Shared memory continue is DISABLED"
    fi
 
    if [ $LOGON -ne 0 ]
    then
        echo "    Log files are ENABLED"
    else
        echo "    Log files are DISABLED"
    fi
 
    if [ $LOGCONTINUE -ne 0 ]
    then
        echo "    Log file continue is ENABLED"
    else
        echo "    Log file continue is DISABLED"
    fi
 
    echo ""
    echo "    SMF Message cache size set to $SMFCACHESIZE records"
    echo ""
}



### Main function

#  Determine if the user wants some help, if they do give it to them
#  and exit.

this_tool=`basename $0`

ARGS=$#
ARGL=$@
count=`expr 1`
HELPFLAG=-h
VERBOSEFLAG=-v
verbose=`expr 0`
PGERETURNFLAG=-p
forcepgeret=`expr 0`
pgerunflag=`expr 0`

while [ $count -le $ARGS ]
do
        ARGX=`echo $ARGL | awk '{print $'$count'}'`
 
        if [ "$ARGX" = "$HELPFLAG" ]
        then
                Usage
                exit
        fi
 
        if [ "$ARGX" = "$VERBOSEFLAG" ]
        then
                verbose=1
        fi
 
        if [ "$ARGX" = "$PGERETURNFLAG" ]
        then
                forcepgeret=1
        fi
 
        count=`expr ${count} + 1`
done


#  If there were not four command line arguments passed in then there 
#  is a problem.
if [ $ARGS -lt 4 ]
then
	Error "must have at least four command line arguments!"
	echo "" 
	Usage
	exit $PGS_SH_SYS_PARAM
fi

#  If the PGE file that was sent in does not exist or has a size of 
#  zero then we consider this a problem.  We are also checking to 
#  determine if the PGE is a regular file as opposed to a directory
#  and if it has execute permission.  The checks are all being performed
#  separately to give the user the proper message.
if [ -s $1 ]
then
	if [ -f $1 ]
	then
		if [ -x $1 ]
		then
			PGEFILE=$1
		else
			Error "PGE file must have execute permission!"
			echo ""
			Usage
			exit $PGS_SH_SYS_PARAM
		fi
	else
		Error "PGE file must be a regular file!"
		echo ""
		Usage
		exit $PGS_SH_SYS_PARAM
	fi
else
	Error "PGE file must exist and have size greater than zero!"
	echo ""
	Usage
	exit $PGS_SH_SYS_PARAM
fi

#  If the Init string is not equal to four characters then we have
#  a problem.  The reason we are check for five here is because wc
#  counts the delimiter.
NUMCHAR=`echo $2 |wc -c`
if [ $NUMCHAR -eq 5 ]
then
	INITSTRING=$2
else
        Error "Init string must contain four characters!"
        echo ""
        Usage
        exit $PGS_SH_SYS_PARAM
fi

#  If the PCF location that was sent in does not exist or has a size of 
#  zero then we consider this a problem.  We are also checking to determine
#  if the file is a regular file as opposed to a directory and if the file
#  has both read and write permissions.  The checks are being perfomed
#  separately to give the user the proper message.
if [ -s $3 ]
then
	if [ -f $3 ]
	then
		if [ -r $3 -a -w $3 ]
		then
			PGS_PC_INFO_FILE=$3
			export PGS_PC_INFO_FILE
		else
			Error "PCF file must be readable and writeable!"
			echo ""
			Usage
			exit $PGS_SH_SYS_PARAM
		fi
	else
		Error "PCF file must be a regular file!"
		echo ""
		Usage
		exit $PGS_SH_SYS_PARAM
	fi
else
	Error "PCF file must exist and have size greater than zero!"
	echo ""
	Usage
	exit $PGS_SH_SYS_PARAM
fi

#  Get the SMF Message cache size from the command line
#  Make sure it's a valid numeric argument
#  If no value is supplied, take the default
if [ "$4" = "" ]
then
        SMFCACHESIZE=50        # default
else
        SMFCACHESIZE=$4
        expr $SMFCACHESIZE + 0 > /dev/null 2>&1
        if [ $? -ne 0 ]
        then
            Error "SMF Cache Size string must be numeric!"
            echo ""
            Usage
            exit $PGS_SH_SYS_PARAM
        fi
fi


#  All input is determined ready to go, let's continue
#  with our processing.
#  Set our Process ID key here and export it.
PGSMEM_SHM_SYSKEY=$$
export PGSMEM_SHM_SYSKEY

#  Break out each digit of the Init String and load them
#  into separate variables for referencing later.
SHMON=`echo $INITSTRING | cut -c1`
SHMON=`expr $SHMON`
SHMCONTINUE=`echo $INITSTRING | cut -c2`
SHMCONTINUE=`expr $SHMCONTINUE`
LOGON=`echo $INITSTRING | cut -c3`
LOGON=`expr $LOGON`
LOGCONTINUE=`echo $INITSTRING | cut -c4`
LOGCONTINUE=`expr $LOGCONTINUE`

#  If the user wants shared memory set the environment variable
#  and the flag that gets passed in to PGS_PC_InitCom.
if [ $SHMON -ne 0 ]
then
	PGSMEM_USESHM=YES
	MEMFLAG=ShmOn
else
	PGSMEM_USESHM=NO
	MEMFLAG=ShmOff
fi
export PGSMEM_USESHM

#  If the user want to continue if shared memory fails let's 
#  make sure we have that variable set correctly.
if [ $SHMCONTINUE -ne 0 ]
then
	SHMCONTINUE=`expr 1`
fi

#  If the user want to continue if the log files fail let's 
#  make sure we have that variable set correctly.
if [ $LOGCONTINUE -ne 0 ]
then
	LOGCONTINUE=`expr 1`
fi

#  If the user wants to use the log files let's make sure we
#  set the variables that get passed in to PGS_PC_InitCom. 
if [ $LOGON -ne 0 ]
then
	LOGFLAG=LogOn
else
	LOGFLAG=LogOff
fi


#  Let the user know what values will be used in the run
#  Run PGS_PC_InitCom.  Capture the return value for use later.
#
if [ $verbose -ne 0 ]
then
    Message
    echo "$this_tool: running PGS_PC_InitCom"
    echo ""
fi
 
 
#  Run PGS_PC_InitCom.  Capture the return value for use later.
PGS_PC_InitCom $MEMFLAG $LOGFLAG $SMFCACHESIZE
RETVAL=`expr $?`

#  Display result of PGS_PC_InitCom.
if [ $verbose -ne 0 ]
then
    echo "$this_tool: return from PGS_PC_InitCom = $RETVAL"
    echo ""
fi

if [ $RETVAL -ne $PGS_S_SUCCESS ]
then
    PGS_EXECUTE_PGE=`expr 0`
fi

#  If the return value was a Shared Memory Init problem then
#  reset our Shared Memory flags.
if [ $RETVAL -eq $PGS_SH_MEM_INIT -o $RETVAL -eq $PGS_SH_MEM_INIT2 ]
then
	PGSMEM_USESHM=NO
	MEMFLAG=ShmOff
	export PGSMEM_USESHM
	if [ $SHMCONTINUE -eq 1 ]
	then
	    PGS_EXECUTE_PGE=`expr 1`
	fi
fi

#  If the return value was a SMF Log File problem then
#  reset our Log File flags.
if [ $RETVAL -eq $PGS_SH_SMF_LOGFILE -o $RETVAL -eq $PGS_SH_SMF_LOGFILE2 -o \
     $RETVAL -eq $PGS_SH_SMF_SHMMEM -o $RETVAL -eq $PGS_SH_SMF_SHMMEM2 -o \
     $RETVAL -eq $PGS_SH_SMF_MSSLOGFILE -o $RETVAL -eq $PGS_SH_SMF_MSSLOGFILE2 ]
then
	LOGFLAG=LogOff
	if [ $LOGCONTINUE -eq 1 ]
	then
	    PGS_EXECUTE_PGE=`expr 1`
	fi
fi

#  The following if statement checks for each circumstance in which we
#  want to run the PGE.  If any of these are successful then we do
#  want to run the PGE.  In any case we do want to run PGS_PC_TermCom.

if [ $PGS_EXECUTE_PGE -eq 1 ]
then
    if [ $verbose -ne 0 ]
    then
        echo "$this_tool: running $PGEFILE"
        echo ""
    fi

    $PGEFILE
    RETPGE=`expr $?`
    pgerunflag=1

#  Display return from the PGE.
    if [ $verbose -ne 0 ]
    then
        echo "$this_tool: return from PGE file:  $PGEFILE = $RETPGE"
        echo ""
    fi

fi

if [ $verbose -ne 0 ]
then
    echo "$this_tool: running PGS_PC_TermCom"
    echo ""
fi

#  Now we get to run PGS_PC_TermCom.
PGS_PC_TermCom $MEMFLAG $LOGFLAG
RETVAL2=`expr $?`

#  Display result of PGS_PC_TermCom.
if [ $verbose -ne 0 ]
then
    echo "$this_tool: return from PGS_PC_TermCom = $RETVAL2"
    echo ""
fi

#  Exit PGS_PC_Shell.sh
if [ $forcepgeret -ne 0 -a $pgerunflag -ne 0 ]
then
    exit $RETPGE
else
    if [ $RETVAL -ne 0 ]
    then
        exit $RETVAL
    else		# return Termination errors if all else OK.	
        exit $RETVAL2
    fi
fi
