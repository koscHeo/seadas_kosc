#include "FileStager.h"
#include <stdint.h>

namespace Aq {
  namespace Fw {
    namespace Ws {
      namespace Io {

	using namespace std;

	/**
	 *  These Error functions are for printing errors in the event that 
	 *  something goes wrong in opening or reading the file.   These are 
	 *  not the flight versions but have the same signature.
	 */
	void Error(const char* s) { printf("%s", s); }
	void Error(const char* s, int n1) { printf("%s %x\n", s, n1); }
        // void Error(const char* s, int n1, int n2) { printf("%s %x\n", s, n1, n2); }



	/**
	 *  This class fetches a file, allocates a buffer for it and reads 
	 *  it into memory up front, and maintains a pointer to the next 
         *  still-unused data in its buffer.
	 */
	FileStager::FileStager(){
	  bufptr = 0;
	  startbyte = 0;
	};

	FileStager::~FileStager(){
	  delete buffer;
	};


	// *** stageFile *** //
	NumBytesAvail FileStager::stageFile(const char* const name, 
					    NameType ntype) {
	  fstream filestr;

	  // Open file for input in binary mode
	  if (ntype == NAME_IS_FILENAME) {
	    filestr.open(name, fstream::in | fstream::binary);
	    if (!filestr.is_open()) {
	      Error("Error opening: ");
	      Error(name);
	      Error("\n");
	      return 0;
	    }
	  }

	  if (ntype == NAME_IS_ENV_VAR) {
	    char * envVar;
	    envVar = getenv (name);
	    filestr.open(envVar, fstream::in | fstream::binary);
	    if (!filestr.is_open()) {
	      Error("Error opening: ");
	      Error(envVar);
	      Error("\n");
	      return 0;
	    }
	  }

	  // Get filesize
	  filestr.seekg (0, ios::end);
	  totalbytes = filestr.tellg();
	  filestr.seekg (0, ios::beg);

	  // Allocate & clear buffer (4-byte aligned)
	  if ((totalbytes % 4) != 0) {
	    buffer = new char[4*((totalbytes/4) + 1)];
	    memset(buffer, 0, 4*((totalbytes/4) + 1));
	  } else {
	    buffer = new char[totalbytes];
	    memset(buffer, 0, totalbytes);
	  }

	  // Read entire file into buffer
	  filestr.read (buffer, totalbytes);

	  if (filestr.eof()) printf("bad\n");

	  // Close file
	  filestr.close();

	  bufptr = buffer;

	  // Set totalbyes to 4-byte-aligned value
	  if ((totalbytes % 4) != 0) totalbytes = 4*((totalbytes/4) + 1);

	  // Return total bytes read
	  return totalbytes;
	}


	// *** stageFile (with magic word search) *** //
	NumBytesAvail FileStager::stageFile(const char* const name, 
					    NameType ntype, int32_t magic) {

	  // Stage file & return number of bytes in stream buffer
	  NumBytesAvail n = this->stageFile(name, ntype);

	  // Find magic word
	  startbyte = 0;
	  bool found = false;
	  /*
	  for (int32_t j=0; j<totalbytes/4; j++) {
	    int32_t i;
	    memcpy(&i, &buffer[4*j], sizeof(int32_t));
	  */
	  // Step every byte rather than 4-bytes JMG 08/14/2011
	  for (int32_t j=0; j<totalbytes; j++) {
	    int32_t i;
	    memcpy(&i, &buffer[j], sizeof(int32_t));

	    // If magic word found then break else increase startbyte counter
	    if (i == magic) {
	      found = true;
	      break;
	    } else {
	      //	      startbyte += sizeof(int32_t);
	      startbyte += sizeof(uint8_t);
	    }
	  }
	    
	  // If magic word not found, exit
	  if (found == false) {
	      unsigned char *ptr;
	      ptr = (unsigned char *) &magic;
	      char strbuf[64];
	      // Write magic word in hex format
	      sprintf(strbuf, "Magic word: \"%x%x%x%x\" not found.\n",
		      (int) ptr[0], (int) ptr[1], (int) ptr[2], (int) ptr[3]);
	      Error(strbuf);
	      exit(2);
	  }

	  // Set buffer pointer to startbyte
	  bufptr = &buffer[startbyte];

	  // Return number of bytes left in buffer 
	  return (n - startbyte);
	}


	// *** getNext *** //
	// Seek
	char *FileStager::getNext(int bytesToRead) {

	  bufptr += bytesToRead;

	  if (bufptr >= &buffer[totalbytes]) {
	    return NULL;
	  } else {
	    // Return current streambuffer position
	    return bufptr;
	  }
	}


	// *** getNext *** //
	// Read
	char *FileStager::getNext(int bytesToRead, char *varbuf) {

	  // Record current stream buffer position
	  int pos = bufptr - buffer;

	  if ((totalbytes - pos) < bytesToRead) {
	    return(0);
	  }

	  // Read next bytesToRead bytes into variable
	  memcpy(varbuf, bufptr, bytesToRead);
	  bufptr += bytesToRead;

	  return bufptr;
	}


	// *** findSync *** //
	// Seek to next sync word
	char *FileStager::findSync( int32_t magic, int32_t offset) {

	  int pos = bufptr - buffer;

	  for (int32_t j=pos+1; j<totalbytes; j++) {
	    int32_t i;
	    memcpy(&i, &buffer[j], sizeof(int32_t));

	    // If magic word found then break else increase startbyte counter
	    if (i == magic) {
	      bufptr = buffer + j + offset;
	      break;
	    }
	  }

	  if (bufptr >= &buffer[totalbytes]) {
	    return NULL;
	  } else {
	    // Return current streambuffer position
	    return bufptr;
	  }
	}


	char *FileStager::rewind() {
	  bufptr = buffer;
	  return bufptr;
	}


	char *FileStager::rewind(int32_t magic) {

	  // Find magic word
	  startbyte = 0;
	  bool found = false;
	  /*
	  for (int32_t j=0; j<totalbytes/4; j++) {
	    int32_t i;
	    memcpy(&i, &buffer[4*j], sizeof(int32_t));
	  */
	  // Step every byte rather than 4-bytes JMG 08/14/2011
	  for (int32_t j=0; j<totalbytes; j++) {
	    int32_t i;
	    memcpy(&i, &buffer[j], sizeof(int32_t));

	    // If magic word found then break else increase startbyte counter
	    if (i == magic) {
	      found = true;
	      break;
	    } else {
	      //	      startbyte += sizeof(int32_t);
	      startbyte += sizeof(uint8_t);
	    }
	  }
	    
	  // If magic word not found, exit
	  if (found == false) {
	      unsigned char *ptr;
	      ptr = (unsigned char *) &magic;
	      char strbuf[64];
	      // Write magic word in hex format
	      sprintf(strbuf, "Magic word: \"%x%x%x%x\" not found.\n",
		      (int) ptr[0], (int) ptr[1], (int) ptr[2], (int) ptr[3]);
	      Error(strbuf);
	      exit(2);
	  }

	  // Set buffer pointer to startbyte
	  bufptr = &buffer[startbyte];

	  return bufptr;
	}


	// *** bytesLeft *** //
	NumBytesAvail FileStager::bytesLeft() const {

	  // Record current stream buffer position
	  int pos = bufptr - buffer;

	  return (totalbytes - pos);
	}


	// *** numRecs *** //
	// Determine # of records from record size
	int FileStager::numRecs(int recSize) {
	  recordSize = recSize;

	  int n = (totalbytes - startbyte) / recordSize;

	  if (n * recordSize != (totalbytes - startbyte))
	    printf("Warning: Incomplete records in data buffer\n");

	  return n;
	}



	// *** createFile *** //
	NumBytesAvail FileStager::createFile(const char* const name, 
					     NameType ntype,
					     int bytesToWrite, char *buffer) {
	  fstream filestr;

	  // Open file for input in binary mode
	  if (ntype == NAME_IS_FILENAME) {
	    filestr.open(name, fstream::out | fstream::binary);
	    if (!filestr.is_open()) {
	      Error("Error opening: ");
	      Error(name);
	      Error("\n");
	      return 0;
	    }
	  }

	  if (ntype == NAME_IS_ENV_VAR) {
	    char * envVar;
	    envVar = getenv (name);
	    filestr.open(envVar, fstream::out | fstream::binary);
	    if (!filestr.is_open()) {
	      Error("Error opening: ");
	      Error(envVar);
	      Error("\n");
	      return 0;
	    }
	  }

	  // Read entire file into buffer
	  filestr.write (buffer, bytesToWrite);

	  // Close file
	  filestr.close();

	  // Return total bytes written
	  return bytesToWrite;
	}

      }
    }
  }
}


