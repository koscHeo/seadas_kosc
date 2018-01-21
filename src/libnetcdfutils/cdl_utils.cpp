#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <libgen.h>

#include <sstream>
#include <iomanip> 
#include "cdl_utils.h"
#include "netcdf.h"
#include "timeutils.h"
#include "dfutils.h"

ncdfFile::ncdfFile() {
  ncid = -1;
}


ncdfFile::~ncdfFile() {
  
}

 /*----------------------------------------------------------------- */
 /* Create an Generic NETCDF4 level1 file                            */
 /* ---------------------------------------------------------------- */
int ncdfFile::cdlCreate( char* l1_filename, char* cdl_filename, 
                         int32_t numScans) {

   int status;
   idDS ds_id;
   status = nc_create( l1_filename, NC_NETCDF4, &ncid);
   check_err(status,__LINE__,__FILE__);

   ifstream cdl_data_structure;
   string line;
   string dataStructureFile( cdl_filename);

   expandEnvVar( &dataStructureFile);

   cdl_data_structure.open( dataStructureFile.c_str(), ifstream::in);
   if ( cdl_data_structure.fail() == true) {
     cout << "\"" << dataStructureFile.c_str() << "\" not found" << endl;
     exit(1);
   }

   // Find "dimensions" section of CDL file
   while(1) {
     getline( cdl_data_structure, line);

     size_t firstNonBlank = line.find_first_not_of(" ");
     if ( firstNonBlank != string::npos)
       if (line.compare( firstNonBlank, 2, "//") == 0) continue;

     size_t pos = line.find("dimensions:");
     if ( pos == 0) break;
   }

   // Define dimensions from "dimensions" section of CDL file
   ndims = 0;
   while(1) {
     getline( cdl_data_structure, line);

     size_t firstNonBlank = line.find_first_not_of(" ");
     if ( firstNonBlank != string::npos)
       if (line.compare( firstNonBlank, 2, "//") == 0) continue;

     size_t pos = line.find(" = ");
     if ( pos == string::npos) break;

     uint32_t dimSize;
     istringstream iss(line.substr(pos+2, string::npos));
     iss >> dimSize;

     iss.clear(); 
     iss.str( line);
     iss >> skipws >> line;

     cout << "Dimension Name: " << line.c_str() << " Dimension Size: "
          << dimSize << endl;

     status = nc_def_dim( ncid, line.c_str(), dimSize, &dimid[ndims++]);
     check_err(status,__LINE__,__FILE__);
   } // while loop

   ngrps = 0;
   // Loop through groups
   while(1) {
     getline( cdl_data_structure, line);

     size_t firstNonBlank = line.find_first_not_of(" ");
     if ( firstNonBlank != string::npos)
       if (line.compare( firstNonBlank, 2, "//") == 0) continue;

     // Check if end of CDL file
     // If so then close CDL file and return
     if (line.substr(0,1).compare("}") == 0) {
       cdl_data_structure.close();
       return 0;
     }

     // Check for beginning of new group
     size_t pos = line.find("group:");

     // If found then create new group and variables
     if ( pos == 0) {

       // Parse group name
       istringstream iss(line.substr(6, string::npos));
       iss >> skipws >> line;
       cout << "Group: " << line.c_str() << endl;

       // Create NCDF4 group
       status = nc_def_grp( ncid, line.c_str(), &this->gid[ngrps]);
       check_err(status,__LINE__,__FILE__);

       ngrps++;

       int numDims=0;
       int varDims[NC_MAX_DIMS];
       size_t dimSize[NC_MAX_DIMS];
       char dimName[NC_MAX_NAME+1];
       string sname;
       string lname;
       string standard_name;
       string units;
       string flag_values;
       string flag_meanings;
       double valid_min=0.0;
       double valid_max=0.0;
       double fill_value=0.0;

       int ntype=0;

       // Loop through datasets in group
       // Skip until "variables:" found
       while(1) {
         getline( cdl_data_structure, line);
         if ( line.find( "variables:") != string::npos) break;
       }

       while(1) {
         getline( cdl_data_structure, line);

         if (line.length() == 0) continue;
         if (line.substr(0,1).compare("\r") == 0) continue;
         if (line.substr(0,1).compare("\n") == 0) continue;

         size_t firstNonBlank = line.find_first_not_of(" ");
         if ( firstNonBlank != string::npos)
           if (line.compare( firstNonBlank, 2, "//") == 0) continue;

         size_t pos = line.find(":");

         // No ":" found, new dataset or empty line or end-of-group
         if ( pos == string::npos) {

           if ( numDims > 0) {
             // Create previous dataset
             createNCDF( gid[ngrps-1], 
                         sname.c_str(), lname.c_str(),
                         standard_name.c_str(), units.c_str(),
                         (void *) &fill_value, 
                         flag_values.c_str(), flag_meanings.c_str(),
                         valid_min, valid_max, ntype, numDims, varDims);

             flag_values.assign("");
             flag_meanings.assign("");
           }

           valid_min=0.0;
           valid_max=0.0;
           fill_value=0.0;

           if (line.substr(0,10).compare("} // group") == 0) break;

           // Parse variable type
           string varType;
           istringstream iss(line);
           iss >> skipws >> varType;

           // Get corresponding NC variable type
           if ( varType.compare("char") == 0) ntype = NC_CHAR;
           else if ( varType.compare("byte") == 0) ntype = NC_BYTE;
           else if ( varType.compare("short") == 0) ntype = NC_SHORT;
           else if ( varType.compare("int") == 0) ntype = NC_INT;
           else if ( varType.compare("long") == 0) ntype = NC_INT;
           else if ( varType.compare("float") == 0) ntype = NC_FLOAT;
           else if ( varType.compare("real") == 0) ntype = NC_FLOAT;
           else if ( varType.compare("double") == 0) ntype = NC_DOUBLE;
           else if ( varType.compare("ubyte") == 0) ntype = NC_UBYTE;
           else if ( varType.compare("ushort") == 0) ntype = NC_USHORT;
           else if ( varType.compare("uint") == 0) ntype = NC_UINT;
           else if ( varType.compare("int64") == 0) ntype = NC_INT64;
           else if ( varType.compare("uint64") == 0) ntype = NC_UINT64;

           // Parse short name (sname)
           pos = line.find("(");
           size_t posSname = line.substr(0, pos).rfind(" ");
           sname.assign(line.substr(posSname+1, pos-posSname-1));
           cout << "sname: " << sname.c_str() << endl;

           // Parse variable dimension info
           this->parseDims( line.substr(pos+1, string::npos), 
                            &numDims, varDims);
           for (int i=0; i<numDims; i++) { 
             nc_inq_dim( ncid, varDims[i], dimName, &dimSize[i]);
             cout << line.c_str() << " " << i << " " << dimName
                  << " " << dimSize[i] << endl;
           }

         } else {
           // Parse variable attributes
           size_t posEql = line.find("=");
           size_t pos1qte = line.find("\"");
           size_t pos2qte = line.substr(pos1qte+1, string::npos).find("\"");
           //    cout << line.substr(pos+1, posEql-pos-2).c_str() << endl;

           string attrName = line.substr(pos+1, posEql-pos-2);

           // Get long_name
           if ( attrName.compare("long_name") == 0) {
             lname.assign(line.substr(pos1qte+1, pos2qte));
             //             cout << "lname: " << lname.c_str() << endl;
           }

           // Get units
           else if ( attrName.compare("units") == 0) {
             units.assign(line.substr(pos1qte+1, pos2qte));
             //             cout << "units: " << units.c_str() << endl;
           }

           // Get _FillValue
           else if ( attrName.compare("_FillValue") == 0) {
             iss.clear(); 
             iss.str( line.substr(posEql+1, string::npos));
             iss >> fill_value;
             //             cout << "_FillValue: " << fill_value << endl;
           }

           // Get flag_values
           else if ( attrName.compare("flag_values") == 0) {
             flag_values.assign(line.substr(pos1qte+1, pos2qte));
           }

           // Get flag_meanings
           else if ( attrName.compare("flag_meanings") == 0) {
             flag_meanings.assign(line.substr(pos1qte+1, pos2qte));
           }

           // Get valid_min
           else if ( attrName.compare("valid_min") == 0) {
             iss.clear(); 
             iss.str( line.substr(posEql+1, string::npos));
             iss >> valid_min;
             //             cout << "valid_min: " << valid_min << endl;
           }

           // Get valid_max
           else if ( attrName.compare("valid_max") == 0) {
             iss.clear(); 
             iss.str( line.substr(posEql+1, string::npos));
             iss >> valid_max;
             //             cout << "valid_max: " << valid_max << endl;
           }

         } // if ( pos == string::npos)
       } // datasets in group loop
     } // New Group loop
   } // Main Group loop
   
   
   return 0;
}


int ncdfFile::parseDims( string dimString, int *numDims, int *varDims) {

  size_t dimSize, curPos=0;
  char dimName[NC_MAX_NAME+1];

  *numDims = 0;

  while(1) {
    size_t pos = dimString.find(",", curPos);
    if ( pos == string::npos) 
      pos = dimString.find(")");

    string varDimName;
    istringstream iss(dimString.substr(curPos, pos-curPos));
    iss >> skipws >> varDimName;

    for (int i=0; i<ndims; i++) {
      int status = nc_inq_dim( ncid, dimid[i], dimName, &dimSize);
      if ( varDimName.compare(dimName) == 0) {
        varDims[(*numDims)++] = dimid[i];
        break;
      }
    }
    if ( dimString.substr(pos, 1).compare(")") == 0) break;

    curPos = pos + 1;
  }

  return 0;
}

int ncdfFile::getGid( const char *grpName) {

  int status;
  int grpID;
  status = nc_inq_grp_ncid( ncid, grpName, &grpID);
  check_err(status,__LINE__,__FILE__);

  return grpID;
}

int ncdfFile::close() {
  
  int status = nc_close(ncid);

  return 0;
}

int createNCDF( int ncid, const char *sname, const char *lname, 
                const char *standard_name, const char *units,
                void *fill_value,
                const char *flag_values, const char *flag_meanings,
                double low, double high, int nt,
                int rank, int *dimids) {

  int32_t varid;
  int i;
  int status;
  size_t dimlength;
  size_t newchunk;
  size_t chunksize[3] = {50,50,50};
     
  /* Create the NCDF dataset */
  status = nc_def_var(ncid, sname, nt, rank, dimids, &varid);
  if( status != NC_NOERR) {
    printf("-E- %s %d: %s for %s\n", 
	   __FILE__, __LINE__, nc_strerror(status), sname);
    exit(1);
  } 

  // Set fill value
  double fill_value_dbl;
  memcpy( &fill_value_dbl, fill_value, sizeof(double));

  int8_t i8;
  uint8_t ui8;
  int16_t i16;
  int32_t i32;
  float f32;

  if ( (low < high) && (low != fill_value_dbl)) {
    if ( nt == NC_BYTE) {
      i8 = fill_value_dbl;
      status = nc_def_var_fill( ncid, varid, 0, (void *) &i8);
    } else if ( nt == NC_UBYTE) {
      ui8 = fill_value_dbl;
      status = nc_def_var_fill( ncid, varid, 0, (void *) &ui8);
    } else if ( nt == NC_SHORT) {
      i16 = fill_value_dbl;
      status = nc_def_var_fill( ncid, varid, 0, (void *) &i16);
    } else if ( nt == NC_INT) {
      i32 = fill_value_dbl;
      status = nc_def_var_fill( ncid, varid, 0, (void *) &i32);
    } else if ( nt == NC_FLOAT) {
      f32 = fill_value_dbl;
      status = nc_def_var_fill( ncid, varid, 0, (void *) &f32);
    } else {
      status = nc_def_var_fill( ncid, varid, 0, (void *) &fill_value_dbl);
    }
    check_err(status,__LINE__,__FILE__);
  }

#if 0 
  /* vary chunck size based on dimensions */ 
  if (rank > 1){
    for (i = 0; i<rank; i++){
      status = nc_inq_dimlen(ncid, dimids[i], &dimlength);
      newchunk = floor(dimlength/25);
      // newchunk = floor(dimlength/12.5);
      if (newchunk > chunksize[i]){
        if ( newchunk < 250){
	  chunksize[i] = newchunk;
        } else {
	  chunksize[i] = 250;
        }
      }
    }
  }
  /* decide whether it is worth compression - dims must be larger than chunks */
  int do_deflate = 1;
  for (i = 0; i<rank; i++){
      status = nc_inq_dimlen(ncid, dimids[i], &dimlength);
      if (dimlength < chunksize[i]) {
          do_deflate = 0;
          break;
      }
  }

  /* Set compression */
    if (ds_id.deflate > 0 && do_deflate ) {
        /* First set chunking */
        status = nc_def_var_chunking(ncid, varid, NC_CHUNKED, chunksize);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n", __FILE__, __LINE__,
                    nc_strerror(status), sname);
            exit(1);
        }

        /* Now we can set compression */
        status = nc_def_var_deflate(ncid, varid, NC_NOSHUFFLE, 1,
                ds_id.deflate);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n", __FILE__, __LINE__,
                    nc_strerror(status), sname);
            exit(1);
        }
    }
#endif

  /* Add a "long_name" attribute */
  status = nc_put_att_text(ncid, varid, "long_name", strlen(lname), lname);
  if( status != NC_NOERR) {
    printf("-E- %s %d: %s for %s\n", 
	   __FILE__, __LINE__, nc_strerror(status), "long_name");
    exit(1);
  } 

  /* Add a "flag_values" attribute if specified*/
  if ( strcmp( flag_values, "") != 0) {
    status = nc_put_att_text(ncid, varid, "flag_values", 
                             strlen(flag_values), flag_values);
    if( status != NC_NOERR) {
      printf("-E- %s %d: %s for %s\n", 
             __FILE__, __LINE__, nc_strerror(status), "flag_values");
      exit(1);
    } 
  }

  /* Add a "flag_meanings" attribute if specified*/
  if ( strcmp( flag_meanings, "") != 0) {
    status = nc_put_att_text(ncid, varid, "flag_meanings", 
                             strlen(flag_meanings), flag_meanings);
    if( status != NC_NOERR) {
      printf("-E- %s %d: %s for %s\n", 
             __FILE__, __LINE__, nc_strerror(status), "flag_meanings");
      exit(1);
    } 
  }

  /* Add a "valid_range" attribute if one is specified */
  if (low < high) {
    switch(nt) {              /* Use the appropriate number type */
    case NC_BYTE:
      {
	uint8_t vr[2];
	vr[0] = (uint8_t)low;
	vr[1] = (uint8_t)high;
	status = nc_put_att_uchar(ncid, varid,"valid_min",NC_BYTE,1,&vr[0]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_min");
	  exit(1);
	} 
	status = nc_put_att_uchar(ncid, varid,"valid_max",NC_BYTE,1,&vr[1]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_max");
	  exit(1);
	} 
      }
      break;
    case NC_UBYTE:
      {
	uint8_t vr[2];
	vr[0] = (uint8_t)low;
	vr[1] = (uint8_t)high;
	status = nc_put_att_uchar(ncid, varid,"valid_min",NC_UBYTE,1,&vr[0]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_min");
	  exit(1);
	} 
	status = nc_put_att_uchar(ncid, varid,"valid_max",NC_UBYTE,1,&vr[1]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_max");
	  exit(1);
	} 
      }
      break;
    case NC_SHORT:
      {
	int16 vr[2];
	vr[0] = (int16_t)low;
	vr[1] = (int16_t)high;
	status = nc_put_att_short(ncid, varid,"valid_range",NC_SHORT,1,&vr[0]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_min");
	  exit(1);
	} 
	status = nc_put_att_short(ncid, varid,"valid_max",NC_SHORT,1,&vr[1]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_max");
	  exit(1);
	} 
      }
      break;
    case NC_INT:
      {
	int32_t vr[2];
	vr[0] = (int32_t)low;
	vr[1] = (int32_t)high;
	status = nc_put_att_int(ncid, varid,"valid_min",NC_INT,1,&vr[0]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_min");
	  exit(1);
	} 
	status = nc_put_att_int(ncid, varid,"valid_max",NC_INT,1,&vr[1]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_max");
	  exit(1);
	} 
      }
      break;
    case NC_FLOAT:
      {
	float vr[2];
	vr[0] = (float)low;
	vr[1] = (float)high;
	status = nc_put_att_float(ncid, varid,"valid_min",NC_FLOAT,1,&vr[0]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_min");
	  exit(1);
	} 
	status = nc_put_att_float(ncid, varid,"valid_max",NC_FLOAT,1,&vr[1]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_max");
	  exit(1);
	} 
      }
      break;
    case NC_DOUBLE:
      {
	double vr[2];
	vr[0] = low;
	vr[1] = high;
	status = nc_put_att_double(ncid, varid,"valid_min",NC_DOUBLE,1,&vr[0]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_min");
	  exit(1);
	} 
	status = nc_put_att_double(ncid, varid,"valid_max",NC_DOUBLE,1,&vr[1]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_max");
	  exit(1);
	} 
      }
      break;
    default:
      fprintf(stderr,"-E- %s line %d: ",__FILE__,__LINE__);
      fprintf(stderr,"Got unsupported number type (%d) ",nt);
      fprintf(stderr,"while trying to create NCDF variable, \"%s\", ",sname);
      return(PROGRAMMER_BOOBOO);
    }
  }           
    
  /* Add a "units" attribute if one is specified */
  if(units != NULL && *units != 0) {
    status = nc_put_att_text(ncid, varid, "units", strlen(units), units);
    if( status != NC_NOERR) {
      printf("-E- %s %d: %s for %s\n", 
	     __FILE__, __LINE__, nc_strerror(status), "units");
      exit(1);
    } 
  }

  /* Add a "standard_name" attribute if one is specified */
  if(standard_name != NULL && *standard_name != 0) {
    status = nc_put_att_text(ncid, varid, "standard_name", 
			     strlen(standard_name), standard_name);
    if( status != NC_NOERR) {
      printf("-E- %s %d: %s for %s\n", 
	     __FILE__, __LINE__, nc_strerror(status), "standard_name");
      exit(1);
    } 
  }
  
  return 0;
}

