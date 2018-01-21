//-*- Mode: C++; -*-  
//============================================================================
//  This software contains Caltech/JPL confidential information.
//  
//  Copyright (c) 2006 by the California Institute of Technology. ALL RIGHTS
//  RESERVED. United States Government Sponsorship Acknowledged. Any
//  commercial use must be negotiated with the Office of Technology
//  Transfer at the California Institute of Technology.
//============================================================================

#ifndef AqFwWsIoFileStager_h
#define AqFwWsIoFileStager_h

#include <stdlib.h>
#include <string.h>
#include <fstream>

namespace Aq {
  namespace Fw {
    namespace Ws {
      namespace Io {

	typedef int NumBytesAvail;

	/**
	 *  This class fetches a file, allocates a buffer for it and reads 
	 *  it into memory up front, and maintains a pointer to the next 
         *  still-unused data in its buffer.
	 */
	class FileStager {
	  char *buffer;
	  char *bufptr;
	  int32_t startbyte;
	  int32_t totalbytes;
	  int recordSize;

	public:
	  FileStager();
	  ~FileStager();

	  enum NameType { NAME_IS_FILENAME, NAME_IS_ENV_VAR };
	  typedef int NumBytesAvail;

	  /**
	   *  This function gets the file name, opens the file, finds out how 
           *  big it is, allocates a buffer for the entire content of the file 
           *  and reads the entire file.  If the file is not found or other 
           *  errors are encountered, this function returns prints appropriate 
           *  messages using one of the Error functions, and returns 0.
	   *
	   *  The name parameter is either a filename or the name of an 
           *  environment variable that is to be read to get a filename, 
           *  depending on the value of the ntype variable.
	   */
	  NumBytesAvail stageFile(const char* const name, NameType ntype);
	  NumBytesAvail stageFile(const char* const name, NameType ntype,
				  int32_t magic);


	  /**
	   *  This function returns the number of bytes left in the buffer - 
           *  that is, the number that the user still hasn't read.  (See 
           *  getNext).
	   *
	   */
	  NumBytesAvail bytesLeft() const;

	  /**
	   *  This function returns a pointer to the next byte that the caller 
           *  has not yet read.   The "bytesToRead" argument tells this object 
           *  how many bytes the caller intends to read, starting with the 
           *  pointer returned by this function.
	   *  On the next call to this function, the returned pointer will be
	   *  "bytesToRead" more than it was at the start of the current call.
           *  Also, the result returned by the next call to bytesLeft() will 
           *  be reduced by "bytesToRead".
	   *
	   *  If "bytesToRead" is larger than the current value of bytesLeft(),
           *  then that is an error, and this function returns 0.
	   *
	   */
	  char *getNext(int bytesToRead);
	  char *getNext(int bytesToRead, char *databuf);

	  char *findSync(int32_t magic, int32_t offset);

	  char *rewind();
	  char *rewind(int32_t magic);

	  int numRecs(int recSize);
          // int numRecs(int32_t magic);

	  NumBytesAvail createFile(const char* const name, 
				   NameType ntype,
				   int bytesToWrite, char *buffer);
	};
      }
    }
  }
}

#endif // AqFwWsIoFileStager_h
