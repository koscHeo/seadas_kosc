#!/bin/sh
#------------------------------------------------------------------------------
#BEGIN_FILE_PROLOG:
#
#FILENAME:
#	update_utcpole.sh
#
#DESCRIPTION:
# This script updates the utcpole.dat file by ftping to USNO
# and reformatting the information into the utcpole file:
# $PGSHOME/database/common/CSC/utcpole.dat
#
#AUTHOR:
#	Ed Larson /  Space Applications Corp
#	Peter D. Noerdlinger /  Space Applications Corp
#       David P. Heroux /  Space Applications Corp
#
#HISTORY:
#	16-Jun-95 EL Initial Version
#	06-Jul-95 EL Modified Prolog
#	05-Mar-95 PDN Modified script to access USNO data instead of IERS
#	20-Mar-95 PDN Modified script for better safety and messaging 
#	24-Feb-96 PDN added case statement for better messaging
#	07-Jul-97 PDN added section to trap default return from PGS_CSC_UT1_update;
#                     e.g. for case function is missing or corrupted.
#	08-Jul-97 PDN added section to trap case of already checked-out file
#	09-Jul-97 PDN added section to trap case of failed "mv" of new file to old name
#	20-Aug-97 PDN reworded message for return 11, which now includes a check on 
#                     string length before first record is written.
#	04-Sep-97 PDN replaced repeated strings with internal string variables.
#                     Added the user's ID to the mailing list when $USER exists.
#                     Added a check on length of file "finals.data"
#                     Added a check for out-of-order records (#14)
#	14-Sep-97 PDN Arranged to capture and return the return values from
#                     functions called in the script (formerly acted on these
#                     but failed to save them)
#       12-Nov-97 PDN Replaced first line with blank line for AUTOSYS compatibility
#                     Added exit stati in some omitted cases
#                     Updated some comments
#       13-Nov-97 PDN Improved error reporting 
#       13-Nov-97 DPH Added shell variables to capture certain returns
#                     added echo of return statuses
#       09-Dec-97 PDN Altered file length check to accommodate newer, shorter files
#       10-Dec-97 PDN Added check for exit value 15 - preparing for USNO going to 2
#                     lines a day about year 2000 or 201
#       18-Feb-98 PDN Replaced return value of 26 in case of success by UNIX std 0
#       03-Mar-98 PDN added -n flag on the ftp call and changed syntax so no 
#                     .netrc file is needed .  Deleted all references to netrc file
#       06-Sep-99 PDN added #!/bin/sh at the top; removed mails to pnoerdli
#       01-Jun-00 PDN fixed phone number list
#       18-Apr-02 XW  added users options to add user names on the email list and 
#		      to input ftp machine name and userid.
#
#	
#END_FILE_PROLOG:
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
#BEGIN_PROLOG:
#
#TITLE:
#	Update UT1 and polar motion file
#
#NAME:
#	update_utcpole.sh
#
#SYNOPSIS:
#	update_utcpole.sh	
#
#
#FORTRAN:
#	N/A
#
#DESCRIPTION:
#       This script updates the utcpole.dat file by ftping to USNO
#       and reformatting the information into the
#
#       utcpole file: $PGSHOME/database/common/CSC/utcpole.dat
#
#       To maintain a current utcpole.dat, this script must be run 
#       every week. It will, however, update the file to the proper
#       current state, even if you are more than a week behind.
#       If you have missed updating for a long time (generally several
#       years!), there could be a failure due to a gap between the old
#       data and the new USNO data.   Generally, the USNO data (see below)
#       start several years before "present time" - the time you access the
#       data files - so there should be no problem.  If, perchance,
#       the script has not been run for a very long time, or has failed
#       over a very long period, unnoticed, please contact the Toolkit
#       staff for up-to-date files.
#       
#       Update_utcpole.sh calls PGS_CSC_UT1_update, a C program that performs
#       most of the actual update work. 
#       
#	The update is done by collecting the latest information via ftp from 
#       the U.S. Naval Observatory in Washington, D.C.  Their file "finals.data" 
#       in "Series 7" contains information on UT1-UTC and the x and y pole
#       displacements, and the errors. The utcpole.dat header contains the date 
#       listed on the USNO server for the last "finals.data" used to update it. 
#
#       The function PGS_CSC_UT1_update reformats the new information and first 
#       overwrites, and then appends it to the utcpole.dat file.
#
#
#INPUTS:
#	Name		Description			Units	Min	Max
#
#	NONE 
#
#OUTPUTS:
#	Name		Description			Units	Min	Max
# 
#	NONE
#
#RETURNS:
#	
#	n/a
#
#EXAMPLES:
#
#       update_utcpole.sh	
#
#NOTES:
#
#      This update process depends on the constancy of the USNO
#      format and data server location/directory structure. If these change,
#      this script may not work.  The Naval Observatory announces any changes
#      well in advance, through its Series 14 bulletins.
#
#      The USNO anticipates changing "finals.data" to two lines a day, using
#      half-integral Modified Julian Dates in the alternate lines and repeating
#      the same, ordinary civil date in the first 7 characters, each line, with
#      the year rolling over to 00 at 2000 AD. The update function 
#      PGS_CSC_UT1_update has been modified to skip the alternate lines, as the
#      EOSDIS project does not need so much accuracy, and we wish to prevent
#      excessive growth of our data file. New return value 15 supports
#      detection of irregularities in the alternate-lines methodology.
#
#      If the updating process fails you must rerun this script.
#
#REQUIREMENTS:
#	PGSTK-1210
#
#DETAILS:
#	
# "finals.data" is issued twice a week. It contains anywhere from a year's
# worth to several years, depending on improvements that the IERS and USNO
# may have made to old solutions for Earth rotation.  The data labeled "I"
# for final generally extend to the day of posting, and the predicted "P"
# data run 1 year beyond that.  The data used for the utcpole.dat file are
# smoothed data for UT1-UTC and the x and y pole displacements.  The errors
# are also imported for the convenience of users who may want them.
# PGS_CSC_UT1_update reads the existing file utcpole.dat, and then it opens
# "finals.data" which it checks for readable format, for start date, and
# for length. The file is then rewound and these data are used to set up a
# copying/overwriting/appending process such that all old data in utcpole.dat
# are overwritten starting at the beginning of the USNO file; the rest of
# the USNO file is then copied over the old values and, when they are exhausted,
# to extend the data set.
#
#GLOBALS:
#	NONE
#
#FILES:
#	$PGSHOME/database/common/CSC/utcpole.dat
#
#FUNCTIONS_CALLED:
#	PGS_CSC_UT1_update
#
#END_PROLOG:
#------------------------------------------------------------------------------

#preliminary definitions and setup

#head string for e-mail
echo "Do you want to modify the preamble in the emails [yes/no]"
read response
case  "$response" in
    y* | Y* )
echo "Please enter the preamble:"
read user_preamble
headerstr="$user_preamble"
    ;;

    * )
headerstr="WARNING: Failure of shell script \"update_utcpole.sh\". 
This message was generated by the script \"update_utcpole.sh\" .
The file \"\$PGSHOME/database/common/CSC/utcpole.dat\"  
will gradually become out of date. 
Check file header for date of most recent update - older than 30 days 
could mean trouble; older than 60 days means inaccurate geolocation. "
    ;;

    esac

#tail string for end of e-mail

tailstr="Or contact ESDIS SDP Toolkit Staff. 
Tel: 1-800-ECS-DATA (1-800-327-3282)"

#set up mailing list
echo "Do you want to add users to the e-mail recipient list [yes/no]"
read add_user
case  "$add_user" in
    y* | Y* )
echo "Please enter the email:"
read user_email

if [ ${USER}0 != ""0 ] 
then userlist="$USER,landover_PGSTLKIT@raytheon.com, $user_email"

else 
userlist="landover_PGSTLKIT@raytheon.com, $user_email"
fi
    ;;
 
    * )                 
if [ ${USER}0 != ""0 ]
then userlist="$USER"

else
userlist="$user_email"
fi
    ;;
 
    esac

# change directory to the utcpole.dat directory:

cd $PGSHOME/database/common/CSC
# if original file is missing abort with a message
if [ ! -s utcpole.dat ] ; then

mail $userlist << QQQ
$headerstr
Old UT1 and polar motion file not found.
Find UT1-polar-motion file "$PGSHOME/database/common/CSC/utcpole.dat" 
and re-run script. 
$tailstr
QQQ

   exit 20
fi

# get file of recent and predicted data from Naval Observatory 
echo "Are you a DAAC-location [yes/no]"
read user_response
case  "$user_response" in
    y* | Y* )
echo "Please enter your DAAC-specific ftp machine:"
read ftp_machine
echo "Please enter your user name:"
read ftp_user
ftp -n <<GG>holder
#open p0fwi09
#user cmts2
open $ftp_machine
user $ftp_user
quote site maia.usno.navy.mil 
user anonymous EOSuser
prompt
cd ser7
dir
get finals.data
bye
GG
    ;;
    * )
echo "Please wait ...."
ftp -n <<RRR>holder
open maia.usno.navy.mil
user anonymous EOSuser
prompt
cd ser7
dir
get finals.data
bye
RRR
    ;;
    esac

#abort and return status if ftp failed

state=$?

if [ $state != 0 ];  then 

       mail $userlist << FTPB
$headerstr
WARNING: failure of ftp to maia.usno.navy.mil
Check with your system administrator
$tailstr
FTPB
exit $state
fi

#abort if file is missing or short (less than about a year's data)

if [ ! -s finals.data ] || [ `wc -c < finals.data` -lt  50000 ]
then
echo "Failure to get data file from maia.usno.navy.mil"
mail $userlist << ZYXW
$headerstr
The file "finals.data" from ftp is missing or too short. Rerun...
$tailstr
ZYXW
exit 22
 
fi
 
# parse the results of ftp for the date, to be used in the new header

grep "finals.data" holder | awk '{print $5 " " $6 " " $7 " " $8}'  > func_input

# overwrite all the data in the old file whose dates are duplicated in new data
# extend file length as needed to accommodate new data
# growth is a line a day. The "func_input" temp file has the header data

PGS_CSC_UT1_update < func_input

state=$?

echo "Status of PGS_CSC_UT1_update call was ($state)"

# test return value.  

    case $state in
      0)
 
#  return value OK but make one further check to be sure file 
#  reaches 30 days or more into the future

DAYTOTL=`date '+%m %Y' | awk '{ sss = 30*($1-1) + 365.2425*($2-1972) + 30;  print(sss) }'`
DAYTOTL=`echo ${DAYTOTL} | awk -F"." '{print $1}'`

   new_length=` cat utcpole.dat.NEW | wc -l `  ;
   grew=`expr $new_length - $DAYTOTL `  ;
   if [ $grew -ge 0 ] ; then

#      all is A-OK

       /bin/mv -f utcpole.dat.NEW utcpole.dat

   mv_state=$?
   echo "Status of MOVE command was ($mv_state)"

#BEGIN nested if - move file in if allowed - branch as error if mv fails

       if [ $mv_state -eq 0 ] ; then

       /bin/rm func_input holder
       /bin/rm finals.data
       exit 0

      else
#error branch within nested if - used if "mv" fails

state=21
mail $userlist << CBA
$headerstr
WARNING: Could not replace "utcpole.dat" by "/bin/mv" - old one 
retained. Check file permissions and rerun, ...
$tailstr
CBA

#      clean up
       /bin/rm func_input holder
       /bin/rm finals.data
       /bin/rm utcpole.dat.NEW

        fi
#EXITED nested if

   else

       mail $userlist << AAA
$headerstr
WARNING: Undersize "utcpole.dat" file generated-old one retained. Rerun ...
$tailstr
AAA

#      clean up
       /bin/rm func_input holder
       /bin/rm finals.data
       /bin/rm utcpole.dat.NEW
       exit 23
   fi

;;

      2)
       mail $userlist << BBB
$headerstr
WARNING: Failure of program "PGS_CSC_UT1_update".
File "utcpole.dat" found but could not be opened.
$tailstr
BBB

;;

      3)
       mail $userlist << CCC
$headerstr
WARNING: Failure of program "PGS_CSC_UT1_update".
Defective string returned from "dir" within ftp - no month name present.
This could also indicate a failure of ftp. Check with system administrator
$tailstr
CCC

;;


      4)
   mail $userlist << DDD
$headerstr
WARNING: Failure of program "PGS_CSC_UT1_update".
Could not create temporary file "utcpole.dat.NEW" used to hold
the new utcpole.dat file during processing. Check file permissions.
$tailstr
DDD


;;

      5)
   mail $userlist << EEE
$headerstr
WARNING: Failure of program "PGS_CSC_UT1_update".
Trying to create too many records for storage within updating program.
(Should not happen before 2054 AD.) Increase limit MAX_RECS and recompile 
function PGS_CSC_UT1_update.
$tailstr
EEE

;;

      6)
   mail $userlist << FFF
$headerstr
WARNING: Failure of program "PGS_CSC_UT1_update".
Could not close original file:  "utcpole.dat". Rerun
$tailstr
FFF

;;

      7)
   mail $userlist << GGG
$headerstr
Unable to open USNO data file "finals.data" (which presumably 
came in by ftp). Check directory structure, rerun, ..
$tailstr
GGG

;;

      8)
   mail $userlist << HHH
$headerstr
WARNING: Failure of program "PGS_CSC_UT1_update".
Fewer than 8 fields found on a data line of "finals.data". Rerun.
$tailstr
HHH

;;

      9)
   mail $userlist << III
$headerstr
WARNING: Failure of program "PGS_CSC_UT1_update".
Unexpected beginning of NEOS (USNO) data apparently before old utcpole 
data (which started in 1972).  Rerun, check original file "utcpole.dat"
$tailstr
III

;;

      10)
   mail $userlist << JJJ
$headerstr
WARNING: Failure of program "PGS_CSC_UT1_update".
Your old file is too old to update as such; you need a later 
utcpole.dat on hand to process. See the Toolkit ftp site.
$tailstr
JJJ

;;

      11)
   mail $userlist << KKK
$headerstr
WARNING: Failure of program "PGS_CSC_UT1_update".
String length or record length as written disagrees with required length 
for later reading; check LogStatus file for more.
$tailstr
KKK

;;

      12)
   mail  $userlist << LLL
$headerstr
WARNING: Failure of program "PGS_CSC_UT1_update".  Could not 
close data file "finals.data" from USNO after use. Rerun ..
$tailstr
LLL

;;

      13)
   mail $userlist << MMM
$headerstr
WARNING: Failure of program "PGS_CSC_UT1_update".
Could not close output file "utcpole.dat.NEW" . Rerun ...
$tailstr
MMM

;;

       14)
    mail $userlist << OAOA
$headerstr
WARNING: Failure of program "PGS_CSC_UT1_update".
Records in old utcpole file out of order. See User Log File
for details.  Obtain an uncorrupted file from Toolkit 
maintenance staff or web site.
$tailstr
OAOA

;;

       15)
   mail $userlist << OBOB
$headerstr
WARNING: Failure of program "PGS_CSC_UT1_update".
Two successive records in imported UNSO file "finals.data" were
skipped, indicating MJD values too close. See User Log File
for details.  Obtain an uncorrupted file from Toolkit 
maintenance staff or web site.
$tailstr
OBOB

;;

      *)
   mail  $userlist << NNN
$headerstr
WARNING: Failure of program "PGS_CSC_UT1_update".
Probable failure mode: executable "PGS_CSC_UT1_update" missing or
corrupted, or system was otherwise unable to execute program. Rerun ...
$tailstr
NNN

;;

    esac

#   clean up files 
   if [ -f func_input ] ; then /bin/rm func_input ; fi
   if [ -f holder ] ; then /bin/rm holder ; fi
   if [ -f finals.data ] ; then /bin/rm finals.data ; fi
   if [ -f utcpole.dat.NEW ] ; then /bin/rm utcpole.dat.NEW ; fi

   exit $state
