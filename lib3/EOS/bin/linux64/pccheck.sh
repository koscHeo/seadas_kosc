#! /bin/sh
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
#	pccheck.sh
#
#DESCRIPTION:
#	This file contains all script functions neccessary
#	to run the pccheck.sh command.
#
#AUTHOR:
#	Ray Milburn / Applied Research Corp.
#
#HISTORY:
#	16-Sep-94 RM Initial Version
#	26-Sep-94 RM Corrected number alignment problem 
#			bug ECSed00123
#	07-Oct-94 RM Added compare option.
#	17-Feb-95 DH Modified Prolog
#	08-Mar-96 MS Updated usage message to properly output script name.
#
#END_FILE_PROLOG:
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
#BEGIN_PROLOG:
#
#TITLE:
#	Check Process Control Information file. (PCF)
#
#NAME:
#	pccheck.sh
#
#SYNOPSIS:
#	pccheck.sh [-h] <-i user-PCF> [-o numbered-PCF] [-c standard PCF] [-s]	
#
#C:
#	N/A
#
#FORTRAN:
#	N/A
#
#DESCRIPTION:
#	The purpose of this tool is to assist the developer in setting up
#       a Process Control File (PCF).  This utility will help to point out 
#	simple syntax and content errors which might lead to more serious 
#	runtime errors, if left uncorrected. This tool will not, however, 
#	detect errors in logic, nor will it correct PCF files.
#
#INPUTS:
#	Name		Description			Units	Min	Max
#
#	-i <PCF>	The -i flag will be followed
#			by the Process Control Information
#			File.  This flag is mandatory.
#
#	-o <outfile>	The -o flag will be followed by
#			a file name which will be output
#			by this command.  The name of 
#			output file must be a file which
#			does not already exist.  This 
#			flag is	optional.
#
#	-h		Upon receiving the -h flag a 
#			short description of the usage
#			of pccheck.sh will be provided
#			to the user and the command
#			will exit.  
#
#	-c		The -c option will cause a compare
#			to be run against a specfied 
#			template file.  The compare will
#			only compare the reserved Product
#			ID's.
#
#	-s		The -s flag will cause all 
#			output except for the output
#			from the -c flag to be suppressed.
#
#OUTPUTS:
#	Name		Description			Units	Min	Max
# 
#	NONE
#
#RETURNS:
#	0	- Normal completion.
#	1	- Error condition.
#
#EXAMPLES:
#
#	pccheck.sh -i $PGSHOME/runtime/pcf.fil -o out.fil
#	pccheck.sh -o out.fil -i $PGSHOME/runtime/pcf.fil
#	pccheck.sh -i $PGSHOME/runtime/pcf.fil -o out.fil -c $PGSRUN/PC/PCF.v3
#	pccheck.sh -i $PGSHOME/runtime/pcf.fil -c $PGSRUN/PC/PCF.v3 -s
#	pccheck.sh -i in.fil
#	pccheck.sh -h
#
#NOTES:
#	This shell script accepts an input file (PCF) and an optional output
#	file.  The output file will be an exact copy of the input file except
#	that line numbers are inserted into the file.  This output file is
#	provided as a convenience to the user when analyzing the generated 
#	report, which sometimes references line locations in the original PCF. 
#	This utility is also capable of comparing against a "standardized" 
#	PCF file to detect changes that have been made to the SDP Toolkit 
#	specific records (those with reserved logical identifiers in the 
#	10K - 11K range); the optional suppression flag prevents all output, 
#	other than the comparison results, from being reported.
#
#REQUIREMENTS:
#	PGSTK-1313
#
#DETAILS:
#	This shell script parses the input to ensure correctness and will
#	report any input problems to the user.  After ensuring the accuracy
#	of the input, the command then spawns the program which actually 
#	performs the check of the PCF.  Upon completion of the check of the 
#	PCF, if the user has requested an output file, the command will then
#	build an exact copy of the file containing line numbers.
#
#	The program that performs the check of the PCF, flags two types of
#	conditions, one is an error, the second is a warning.  There is no
#	directive stating what these conditions are, so a brief definition
#	is given here.  The pccheck.sh command considers a condition an 
#	error if during execution of any of the PC tools an error condition 
#	will be returned.  The pccheck.sh command considers a condition a
#	warning if during execution of any of the PC tools a successful 
#	condition will be returned, but it is possible that using what is
#	returned as output from any of the PC tools MAY result in an error
#	condition.  An example of an error condition would be a failure to 
#	include the version number in a standard input file data line.  An
#	example of a warning condition would be to have a blank space in a
#	file name.
#
#	The compare option is designed to be used to compare the PCF with the 
#	template PCF delivered with the PGS Toolkit.  The compare option will
#	compare the Product ID's reserved by the PGS Toolkit to find any
#	differences.  The differences displayed are for the user only.  This
#	tool makes no attempt to alter the input PCF.
#
#GLOBALS:
#	NONE
#
#FILES:
#	This script passes an input file to the check program and creates
#	an output file which is identical to the input file except that it
#	contains line numbers.  This script also creates and deletes a 
#	temporary file during the compare option.
#
#FUNCTIONS_CALLED:
#	Usage		- print a usage message
#	Error		- print an error message and return an error status
#	number		- make a copy of an input file and place line numbers
#			  in the copy file
#	pctcheck	- C program which checks the validity of a PCF
#
#END_PROLOG:
#------------------------------------------------------------------------------
# Define return values
: ${OK=0} ${FAIL=1}

# Usage function
Usage()
{
	echo "usage:  `basename $0` [-h] [-s] [-i input file] [-o output file] [-c compare file]"
	echo "-i	file name to be checked, this option must be used."
	echo "-o	output file name (optional)."
	echo "-c	compare against a template file (optional)."
	echo "-s	suppress all output except compare option."
	echo "-h	usage help."
}

# Error function
Error()
{
	echo "\nError:  $*"  >&2
	return $FAIL
}


number()
{
# adapted from the 'num' line-numbering script developed by Bob Stockler 12/01/90
#
# this version operates as a function within a larger script to support
# line-number generation for a file used as input to the master script.
# NOTE: This version assigns line numbers to all except BLANK lines.
#
# The following arguments are expected:
#       arg 1 - number of columns supported for line numbers
#       arg 2 - line number separator
#       arg 3 - original input file
#       arg 4 - line numbers prefixed to input file

	awk '
	BEGIN	{ format = "%" "'"$1"'" "d" "'"$2"'" "%s\n" }
		{ printf( format, ++i, $0 ) > "'"$4"'" } ' $3

}

### Main function

#  If there were no command line arguments passed in then there
#  is a problem.
if [ $# -eq 0 ]
then
	Usage
	echo "" 
	exit $FAIL
fi

#  Initialize our variables here.
INFLAG=NO
OUTFLAG=NO
SUPPFLAG=NO
ARGS=$#
count=`expr 1`

#  Loop the number of command line arguments times.
while [ $count -le $ARGS ]
do

#  Determine what the command line argument is.
	case $1 in

#  The argument is a -i, the user is going to give us our input
#  file name here.  Let's first make sure that there are more
#  arguments on the command line, if not then the user entered a
#  -i with nothing after it.
	-i)	if [ $count -ge $ARGS ]
		then
			Error "input file must be specified."
			Usage
			exit $FAIL

#  We know the user entered something after the -i, let's make
#  sure it is a file that actually exists before we move on.
		else
			if [ -f $2 ]
			then
				INFILE=$2
				INFLAG=YES
				shift
				count=`expr ${count} + 1`
			else
				Error "input file must exist!"
				Usage
				exit $FAIL
			fi
		fi
		;;

#  The argument is a -o, the user is going to give us our output
#  file name here.  Let's first make sure that there are more
#  arguments on the command line, if not then the user entered a
#  -o with nothing after it.
	-o)	if [ $count -ge $ARGS ]
		then
			Error "output flag specified but no file given."
			Usage
			exit $FAIL

#  We know the user entered something after the -o, let's make
#  sure it is a file that does not already exist before we move on.
		else
			if [ -f $2 ]
			then
				Error "output file given already exists!"
				Usage
				exit $FAIL
			else
				OUTFILE=$2
				shift
				count=`expr ${count} + 1`
				OUTFLAG=YES
			fi
		fi
		;;

	-c)	if [ $count -ge $ARGS ]
		then
			Error "template file for comparison must be specified."
			Usage
			exit $FAIL

#  We know the user entered something after the -i, let's make
#  sure it is a file that actually exists before we move on.
		else
			if [ -f $2 ]
			then
				COMPFILE=$2
				COMPFLAG=YES
				shift
				count=`expr ${count} + 1`
			else
				Error "template file for comparison must exist!"
				Usage
				exit $FAIL
			fi
		fi
		;;


	-s)	SUPPFLAG=YES
		;;

#  The user wants help, this option overrides all other options.
#  Let's give them some help and leave.
	-h)	Usage
		exit $OK
		;;
	esac

#  Move to the next command line argument and update our counter
#  so we know how many of these boys we have processed.
	shift
	count=`expr ${count} + 1`
done

#  Here is where we actually run the program the performs the 
#  checking - then we capture the return value sent back from
#  the program.
if [ "$SUPPFLAG" = NO ]
then
	if [ "$INFLAG" = YES ]
	then
		pctcheck $INFILE
		LINES=`expr $?`
	else
		Usage
		exit $FAIL
	fi
fi

if [ "$SUPPFLAG" = NO ]
then
	if [ "$OUTFLAG" = YES -a $LINES -gt 0 ]
	then
		SEP=':'
		COLS=`expr 0`
		while [ $LINES -gt 0 ]
		do
			LINES=`expr ${LINES} / 10`
			COLS=`expr ${COLS} + 1`
		done
 		number $COLS "$SEP" $INFILE $OUTFILE 
	fi
fi

if [ "$COMPFLAG" = YES ]
then
	diff $INFILE $COMPFILE | grep "^[\<\>] *10[0-9][0-9][0-9]|" > fpgspcrm08.tmp
	if [ -s fpgspcrm08.tmp ]
	then
		echo ""
		echo ""
		echo "The following lines were listed in the template file:"
		echo "$COMPFILE"
		echo "and have been altered or deleted from the input file."
		echo ""
		sed -e '/^</d' < fpgspcrm08.tmp
		echo ""
		echo ""
		echo "These are the lines in the input file:"
		echo "$INFILE"
		echo "that differ from the template file."
		echo ""
		sed -e '/^>/d' < fpgspcrm08.tmp
		echo ""
	else
		echo ""
		echo "No differences found in compare check."
		echo ""
	fi
	rm -f fpgspcrm08.tmp
fi

exit $OK
