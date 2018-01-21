
/* determine if year is a leap year or not */
#define TRUE 1
#define FALSE 0

int isleap(int year)
{
if ( ((year % 400) == 0) || (((year % 4) == 0) && ((year % 100) != 0)) )
	return TRUE;
else
	return FALSE;
}
