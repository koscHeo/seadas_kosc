#include <setupflags.h>
#include <string.h>


void setupflags (char *flagdef, char *flaguse, uint32 *flagusemask, 
        uint32 *required, int *status)
{
        int32 BITS[32] = {  0x00000001, 0x00000002, 0x00000004, 0x00000008,
                           0x00000010, 0x00000020, 0x00000040, 0x00000080,
                           0x00000100, 0x00000200, 0x00000400, 0x00000800,
                           0x00001000, 0x00002000, 0x00004000, 0x00008000, 
                           0x00010000, 0x00020000, 0x00040000, 0x00080000, 
                           0x00100000, 0x00200000, 0x00400000, 0x00800000, 
                           0x01000000, 0x02000000, 0x04000000, 0x08000000, 
                           0x10000000, 0x20000000, 0x40000000, 0x80000000  };

	int bitNum;
        char *tmpFlags;
	char *ptr, *ptr2;

	*status = 0;
	*flagusemask = 0;
	*required = 0;

	bitNum = 0;
        tmpFlags = strdup(flagdef);
	ptr = strtok(tmpFlags, ",");
	while ( ptr ) {
	  if ( ptr ) {
	    if ((ptr2 = strstr(flaguse, ptr))) {
	      ptr2--;
	      if (*ptr2 == '~')
		*required = *required | BITS[bitNum];
	      else
		*flagusemask = *flagusemask | BITS[bitNum];
	    }
	  }
	  ptr = strtok(NULL,",");
	  bitNum++;
	  if (bitNum > 33) {
	    *status = -1;
	    break;
	  }
	}
        
        free(tmpFlags);
}

