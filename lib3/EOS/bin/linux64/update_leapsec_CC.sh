#! /bin/sh
#-------------------------------------------------------------------------#
#                                                                         #
#  COPYRIGHT[copyright mark] 2000, Raytheon System Company, its vendors,  #
#  and suppliers.  ALL RIGHTS RESERVED.                                   #
#                                                                         #
#-------------------------------------------------------------------------#
#--------------------------------------------------------------------------
#BEGIN_FILE_PROLOG:
#
#FILENAME:
#	update_leapsec_CC.sh
#
#DESCRIPTION:
#This script updates the leapsec.dat file by ftping to the USNO
#and reformatting their file "ser7/tai-utc.dat" into the leapsec file:
#       $PGSHOME/database/common/TD/leapsec.dat
#
#AUTHOR:
#	Peter D. Noerdlinger / Space Applications Corp.
#       David P. Heroux / Space Applications Corp.
#
#HISTORY:
#	19-Mar-96 PDN Modified Prolog and failure mode response
#	16-Aug-96 PDN Rewrote to work with US Naval Observatory
#                     instead of IERS
#	23-Mar-97 PDN Reworked trapping of bad returns
#	08-Jul-97 PDN Added Clear Case capability; added failure mode for 
#                     cases where executable is missing, corrupt, or cannot 
#                     be executed (e.g. architecture mismatch)
#       04-Sep-97 PDN replaced repeated strings with internal string variables.
#                     Added the user's ID to the mailing list when $USER exists.
#       13-Nov-97 PDN Improved error reporting and changed 1st line per
#                     advice of DPH
#       13-Nov-97 DPH Added shell variables to capture certain returns
#                     added echo of return statuses
#       18-Feb-98 PDN Fixed return in case of success from 25 to 0 (UNIX Std)
#       05-Mar-98 PDN added -n flag to ftp and moved machine name to input
#                     of the ftp; added user and password EOSuser there.
#                     This eliminates need for .netrc file
#       05-Mar-98 PDN changed mail header to show file stale at 83 days
#       10-Sep-98 PDN Restored initial line to force Bourne shell (Toolkit Standard)
#       10-Sep-98 PDN Put in exit status 22 when "mv" fails. Fixed some comments.
#       10-Sep-98 PDN Fixed scratch file cleanup in the case of "exit 22"
#       06-Sep-99 PDN Deleted e-mail to self (pnoerdli)
#       01-Dec-99 PDN inserted the "-cvi" which is now present in 
#                     "cleartool lscheckout -cvi leapsec.dat" to prevent an
#                     aborted run due to checkout in a different view.
#       01-Jun-00 PDN Replaced phone number list with 1-800-ECS-DATA
#
#END_FILE_PROLOG:
#--------------------------------------------------------------------------
#BEGIN_PROLOG:
#
#TITLE:
#	Update Leap Seconds file
#
#NAME:
#	update_leapsec_CC.sh
#
#SYNOPSIS:
#	update_leapsec_CC.sh	
#
#C:
#	N/A
#
#FORTRAN:
#	N/A
#
#DESCRIPTION:
#       This script updates the leapsec.dat file by ftping to USNO
#       and reformatting the information into the
#       leapsec file: $PGSHOME/database/common/TD/leapsec.dat
#       
#       To maintain a current leapsec.dat, this script must be run at 
#       least every three months, but it is safer and just as easy to
#       run it weekly.  The three month requirement stems from the 
#       following considerations: 
#        
#       Leap Seconds are normally introduced about every 12 to 18 months
#       by decision of the  International Earth Rotation Service (IERS) 
#       in Paris, France, in order to keep UTC within 0.9 seconds of UT1,
#       They are normally added at the end of December and June, 
#       but could be added at the end of March or September instead. 
#       (The new leap second value TAI - UTC is effective after mid- 
#       night of that last day of the month, so the leap seconds are 
#       sometimes said to be effective, for example, Jan 1 or July). 
#        
#       While announcements of new leap seconds are nominally made 6 
#       months in advance, there can be a lag of a few days in posting. 
#       To allow for this, and any problems with ftp, and for the 
#       possibility of a March or September date, it is best to run
#       update_leapsec_CC.sh weekly.  If there is no new leap second 
#       announced, the update function called by this script will 
#       update the header to show that the USNO file was checked for 
#       a change, but will do no more.
#
#       The announcements are made via the posting of a new "bulletinc" 
#       by the IERS on its server 145.238.100.28 (or   hpiers.obspm.fr).
#       The change is echoed within a few days by the USNO on their 
#       server maia.usno.mavy.mil =  192.5.41.22, by posting a full file
#       of all leap seconds (historical and new) in the directory "ser7". 
#       Their file name is "tai-utc.dat".
#
#       Because the USNO file contains no incidental remarks, as does
#       "bulletinc", and seems more immune to format changes, and
#       because the IERS moved its server in July, 1996, it was
#       decided in August 1996 to change service to the USNO.
# 
#       The present script, after obtaining the required file invokes
#       PGS_TD_NewLeap, a C program that performs the actual update work. 
#
#       The function puts the current date in the header of the new 
#       leapsec.dat, with a remark that the file was either "Checked"
#       (no new leap second) or "Updated" (new leap second)
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
#	               The script will send e-mail to the SDP Toolkit mail
#                      address in case of failure
#
#EXAMPLES:
#
#       update_leapsec_CC.sh	
#
#NOTES:
#
#      To run this script you must own a Clear Case view and be in it.
#
#      This update process depends on the constancy of the format of the
#      USNO file. If that format changes, this script may not work.
#
#
#REQUIREMENTS:
#	PGSTK-1210
#
#DETAILS:
#	
#	The existing file leapsec.dat contains the historical and announced 
#leap seconds, all denoted actual. The IERS announcements are definitive and 
#result in an actual entry. They do not include what were called "predicted" 
#leap seconds in the SCF version data file. The latter predictions, designed
#to accommodate users needing to test software with far future dates,
#extend to the year 2010. Useful ONLY for simulations, they are eliminated
#in the DAAC Toolkit
#PGS_TD_NewLeap reads the existing leap seconds file, and the USNO file
#maia.usno.navy.mil/ser7/tai-utc.dat,and compares the last leap second in
#each. If the dates are the same, the the header of "leapsec.dat" is updated
#to show that the file was checked.  If the last USNO leap second is later
#than ours, the header is updated to show the change, and all the USNO leap 
#seconds are copied into leapsec.dat in Toolkit format
#
#GLOBALS:
#	NONE
#
#FILES:
#	$PGSHOME/database/common/TD/leapsec.dat
#
#FUNCTIONS_CALLED:
#	PGS_TD_NewLeap
#
#END_PROLOG:
#-------------------------------------------------------------------------------

#preliminary definitions and setup
 
#head string for e-mail
headerstr="WARNING: Failure of shell script \"update_leapsec_CC.sh\". 
This message was generated by the script \"update_leapsec_CC.sh\" .
Hand Check Leap seconds file \"$PGSHOME/database/common/TD/leapsec.dat\"
against USNO data. Use  \"ftp maia.usno.navy.mil\" and change to the ser7 
directory. See the file \"tai-utc.dat\" and edit the Toolkit file as needed
Check \$PGSRUN/LogStatus for more details.
The Leap Seconds File will gradually become out of date. Check file header 
for date of most recent update or successful checking operation - older 
than 83 days means that the file is too stale for use by the Toolkit."

#tail string for end of e-mail
 
tailstr="You may wish to contact ESDIS SDP Toolkit Staff. 
Tel: 1-800-ECS-DATA (1-800-327-3282)"
 
#set up mailing list
if [ ${USER}0 != ""0 ] 
 
then userlist="$USER,landover_PGSTLKIT@raytheon.com"
 
else 
userlist="landover_PGSTLKIT@raytheon.com"
fi

# change directory to the leapsec.dat directory:

cd $PGSHOME/database/common/TD

if  [ ! -s leapsec.dat ] ; then

mail $userlist << EEE
$headerstr
This message was generated by the script "update_leapsec_CC.sh"
Old leap seconds file "$PGSHOME/database/common/TD/leapsec.dat" not found
Find leap seconds file and re-run script, or
$tailstr
EEE
exit 21
fi

# if file is checked out abort with a message
 
isout=`cleartool lscheckout -cvi leapsec.dat`
if [  -n "$isout" ];  then
 
mail $userlist << QAQA
$headerstr
File  "$PGSHOME/database/common/TD/leapsec.dat" is already checked out.
Determine cause, ensure that a reliable recent version is checked in and
re-run script, or ...
$tailstr
QAQA
 
   exit 23
fi

# Check out the file

cleartool co -nc leapsec.dat

  if [ $? != 0 ] ; then

mail $userlist << RPM
$headerstr
Can't check out leapsec.dat - check with your Clear Case administrator or..
$tailstr
RPM
exit 24
  
fi

#
#ftp to USNO  (maia.usno.navy.mil 192.5.41.22 ) and get the current tai-utc.dat in ser7

ftp -n <<GG>holder
open maia.usno.navy.mil
user anonymous EOSuser
prompt
cd ser7
dir
get tai-utc.dat
bye
GG

#parse the results of ftp for the date, to be used in the new header

grep "tai-utc.dat" holder | awk '{print  " " $6 " " $7 " " $8}'  > func_input

#clean up a bit
/bin/rm holder

#run update C program PGS_TD_NewLeap:
#PGS_TD_NewLeap.c assumed compiled and in $PGSBIN (done by the INSTALL script)
# if the executable PGS_TD_NewLeap has been lost or corrupted, "cd" to the
# directory $PGSSRC/TD and at the keyboard enter "make utilities"

   PGS_TD_NewLeap < func_input

   state=$?

echo "Status of PGS_TD_NewLeap call was ($state)"

#test return value

  case $state in

       0)

#All is OK 

   /bin/mv -f leapsec.dat.new $PGSHOME/database/common/TD/leapsec.dat
 
#BEGIN if - check file in if allowed - branch as error if mv fails
 
   mv_state=$?
   echo "Status of MOVE command was ($mv_state)"
 
 if [ -f holder ]; then /bin/rm holder ; fi
 if [ -f tai-utc.dat ]; then /bin/rm tai-utc.dat ; fi
 if [ -f func_input ]; then /bin/rm func_input ; fi

       if [ $mv_state -eq 0 ] ; then

cleartool ci -nc $PGSHOME/database/common/TD/leapsec.dat

       exit 0
 
      else
 
   /bin/rm -f leapsec.dat.new

#error branch within if - used if "mv" fails
 
       mail $userlist << CBA
$headerstr
WARNING: Could not replace "leapsec.dat" by "/bin/mv" - old one 
retained. Check file system, or Clear Case, or ....
$tailstr
CBA
 
cleartool uncheckout leapsec.dat  << CBZZ
n
CBZZ
exit 22

fi
#EXITED if test on bin/mv
 

;;

       2)

mail $userlist << AAA
$headerstr
Failure of program "PGS_TD_NewLeap".  No month name found in 
parsing results of "dir" within "ftp". Rerun recommended.
$tailstr
AAA

cleartool uncheckout leapsec.dat  << PAAA
n
PAAA


;;

       3)

mail $userlist << BBB
$headerstr
Failure of program "PGS_TD_NewLeap": Can't open old "leapsec.dat".
$tailstr
BBB

cleartool uncheckout leapsec.dat  << PBBB
n
PBBB

;;

       4)

mail $userlist << CCC
$headerstr
Failure of program "PGS_TD_NewLeap": Can't open local scratch file.
$tailstr
CCC

cleartool uncheckout leapsec.dat  << PCCC
n
PCCC

;;
      5)

mail $userlist << TYY
$headerstr
Failure of program "PGS_TD_NewLeap": New header was too short.
$tailstr
TYY

cleartool uncheckout leapsec.dat  << PTYY
n
PTYY

;;
       6)


mail $userlist << PUPP
$headerstr
Failure of program "PGS_TD_NewLeap": Data line too short.
$tailstr
PUPP

cleartool uncheckout leapsec.dat  << PUPP
n
PUPP

;;

       7)


mail $userlist << DDD
$headerstr
Failure of program "PGS_TD_NewLeap": New USNO leap second before 
old final one.
$tailstr
DDD

cleartool uncheckout leapsec.dat  << PDDD
n
PDDD

;;

       9)

mail $userlist << EEE
$headerstr
Failure of program "PGS_TD_NewLeap": Can't open data file from USNO 
(result of ftp).
$tailstr
EEE

cleartool uncheckout leapsec.dat  << PEEE
n
PEEE

;;

       12)

mail $userlist << FFF
$headerstr
Failure of program "PGS_TD_NewLeap": Old leap seconds file too long.
Recompile with larger array size.  This event not expected before 
2150 AD
$tailstr
FFF

cleartool uncheckout leapsec.dat  << PFFF
n
PFFF

;;
       13)

mail $userlist << GGG
$headerstr
Failure of program "PGS_TD_NewLeap": Too many records in USNO data file.
This event not expected before 2150 AD. Recompile with larger array size - or
$tailstr
GGG

cleartool uncheckout leapsec.dat  << PGGG
n
PGGG

;;
       14)
 
mail $userlist << GHH
$headerstr
Failure of program "PGS_TD_NewLeap": Problem reading data line in USNO file.
$tailstr
GHH
 
cleartool uncheckout leapsec.dat  << PGHH
n
PGHH
 
;;

       15)

mail $userlist << HHH
$headerstr
Failure of program "PGS_TD_NewLeap": Can't close old leap seconds file.
$tailstr
HHH

cleartool uncheckout leapsec.dat  << PHHH
n
PHHH

;;
       18)

mail $userlist << III
$headerstr
Failure of program "PGS_TD_NewLeap": Can't close data file from 
USNO (result of ftp)
$tailstr
III

cleartool uncheckout leapsec.dat  << PIII
n
PIII

;;

      *)
mail $userlist << NNN
$headerstr
Failure of program "PGS_TD_NewLeap":
Probable failure mode: executable "PGS_TD_NewLeap" missing or
corrupted, or system was otherwise unable to execute program.
$tailstr
NNN
 
cleartool uncheckout leapsec.dat  << PNNN
n
PNNN
 
;;

esac

#clean up files

 if [ -f tai-utc.dat ]; then /bin/rm tai-utc.dat ; fi
 if [ -f func_input ]; then /bin/rm func_input ; fi
 if [ -f leapsec.dat.new ]; then /bin/rm leapsec.dat.new ; fi
 if [ -f holder ]; then /bin/rm holder ; fi

exit $state
