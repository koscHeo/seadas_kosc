/** @file sstval_extract.h
	@brief Process a small section of a Level-2 NetCDF file.
*/

#ifndef ___SSTVAL_EXTRACT_H_
#define ___SSTVAL_EXTRACT_H_

#include "val_extract.h"

#include <stdbool.h>

#ifndef SSTVALEXTRACT_API_VERSION
/** @brief Current API version.  If already defined, use API defined first. */
#define SSTVALEXTRACT_VERSION 1000002
/** @brief Current API version as a string. */
#define SSTVALEXTRACT_VERSION_STR "1.0.2"
#endif

/** @brief Returned on successful processing. */
#define SSTVALEXTRACT_ERR_NONE EXIT_SUCCESS
/** @brief Returned when the desired point is not in the file boundaries. */
#define SSTVALEXTRACT_ERR_POINT_NOT_FOUND 109
/** @brief Returned for unexpected errors like malloc failures or, possibly, permissions
		problems and the like. */
#define SSTVALEXTRACT_ERR_UNKNOWN 110
/** @brief Returned when the NetCDF file can't be opened (due to errors or corruption). */
#define SSTVALEXTRACT_ERR_NCFILE_ERR 111
/** @brief Returned when the NetCDF file isn't in the format expected (not an L2, etc). */
#define SSTVALEXTRACT_ERR_NCFILE_INVALID 112
/** @brief Returned when the something goes wrong processing l2_flags. */
#define SSTVALEXTRACT_ERR_FLAG 113
/** @brief Returned when the something goes wrong processing a product or when finding a
		product specified on the command line. */
#define SSTVALEXTRACT_ERR_VARIABLE 114
/** @brief Returned when given bad or no arguments. */
#define SSTVALEXTRACT_ERR_INPUT 115
/** @brief Not a real error, but returned when the L2QC step says it's a poor quality file. */
#define SSTVALEXTRACT_ERR_L2QC 116

/** @brief argpar structure used for making programs that inherit options from this library. */
extern argpar sstval_extract_argpar;

/** @brief Passed to parser when beginning the processing. */
#define SSTVALEXTRACT_KEY_INIT    1
/** @brief Passed to parser when the data has been processed and is passed a pointer to a
		sstval_extract_output structure. */
#define SSTVALEXTRACT_KEY_OUTPUT  2
/** @brief Passed to parser when the processing is finished and no errors were encountered. */
#define SSTVALEXTRACT_KEY_SUCCESS 3
/** @brief Passed to parser when the processing is finished and an error was encountered. */
#define SSTVALEXTRACT_KEY_ERROR   4
/** @brief Passed to parser when ending the processing. */
#define SSTVALEXTRACT_KEY_FINI    5

/** @typedef typedef struct sstval_extract_output
	@brief For the lazy. */
/** @struct sstval_extract_output
	@brief Passed into the sstval_extract_parser when the data's been processed.
*/
typedef struct sstval_extract_output {
	/** @brief Region within the NetCDF taht was processed. */
	nc_region region;
	/** @brief Options passed in from the user. */
	const char *request_id, *request_file, *buoy_id, *buoy_time;
	bool print_header;
	/** @brief Single character sensor ID, grabbed from the filename. */
	char sensor;
	/** @brief Granule identifier, the date grabbed from the filename, starting with the
		two-digit year. */
	char *granule;
	/** @brief The version identifier pulled from the calibration_data attribute in the NetCDF
		file.  Specifically, it is the version at the end of MYD02_Reflective_LUTs file. */
	char *reflective_lut_version;
	/** @brief Date of the center of the region. */
	unsigned year, day_of_year;
	/** @brief Time in fractions of an hour of the center of the region. */
	double hour;
	/** @brief Mirror side (zero-based). */
	unsigned mirror;
	/** @brief Detector number (zero-based). */
	unsigned detector_number;
	/** @brief Scan pixel number (zero-based). */
	unsigned glint_flag;
	/** @brief The Level-2 processing flags of the center pixel of the region. */
	unsigned pixel_number;
	/** @brief The HIGLINT flag of the center pixel of the region. */
	unsigned l2_flags;
	/** @brief Product-specific flags, Sea Surface Temperature, of the center pixel of the region. */
	unsigned flags_sst;
	/** @brief Product-specific flags, 4um Sea Surface Temperature, of the center pixel of the region. */
	unsigned flags_sst4;
	/** @brief Location of the center pixel of the region. */
	double lat, lon;
	/** @brief Solar zenith angle of the center pixel of the region. */
	double solz;
	/** @brief Sensor zenith angle of the center pixel of the region. */
	double senz;
	/** @brief Angle of incidence of what now?. */
	double angle_of_incidence;
	/** @brief Solar azimuth angle of the center pixel of the region. */
	double sola;
	/** @brief Sensor azimuth angle of the center pixel of the region. */
	double sena;
	/** @brief Relative azimuth angle of the center pixel of the region. */
	double relaz;
	/** @brief Quality levels, Sea Surface Temperature, of the center pixel of the region. */
	double quality_sst;
	/** @brief Quality levels, 4um Sea Surface Temperature, of the center pixel of the region. */
	double quality_sst4;
	/** @brief Statistics for like-named variables found in the NetCDF file. */
	nc_var_stats sst, sst4, bt_3750, bt_3959, bt_4050, bt_6715, bt_7325, bt_8550, bt_11000, bt_12000, rhot_678, rho_cirrus, sstref;
	/** @brief Statistics for pixel-by-pixel subtraction of bt_3959 and bt_4050. */
	nc_var_stats bt3959minus4050;

	/** @brief  Full variables, with data, for their like-named NetCDF variables. */
	nc_var *bt_3959_var, *bt_4050_var, *pixel_number_var;
} sstval_extract_output;

/** @brief Pointer to a callback function to call for each argument parsed.

	@param[in] key One of the SSTVALEXTRACT_KEY_ macros.
	@param[in] sstval_input Either NULL or a sstval_extract_output, depending on they key.
	@param[in] user_input Input pointer from the sstval_extract_arguments structure.

	@return 0 on success, 1 on error.
*/
typedef int (*sstval_extract_parser)(int key, void *sstval_input, void *user_input);

/** @typedef typedef struct sstval_extract_arguments
	@brief For the lazy. */
/** @struct sstval_extract_arguments
	@brief Passed into the library function to control the processing.
*/
typedef struct sstval_extract_arguments {
	/** @brief Print the header before data.  If true, an input file for val-extract is not required. */
	bool header;
	/** @brief Request ID that goes into the header.  Not required if header is false. */
	const char *req_id;
	/** @brief Request file that goes into the header.  Not required if header is false. */
	const char *req_file;
	/** @brief Buoy ID to print with the header. */
	const char *buoy_id;
	/** @brief Buoy time to print with the header, in the format of YYYYDDDHHMMSS. */
	const char *buoy_time;
	/** @brief Arguments for val-extract to control where the data comes from. */
	val_extract_arguments val_extract_arguments;

	/** @brief Required, this is the function that is given the sstval_extract_output once the data is
		processed. */
	sstval_extract_parser sstval_extract_parser;
	/** @brief If given, this is passed to the sstval_extract_parser for arbitrary use by the user. */
	void *user_input;
} sstval_extract_arguments;

/** @brief Process a small section of a Level-2 NetCDF file.

	For actual usage arguments, run argc=0 and review STDOUT.

	@param[in] argc Argument count.
	@param[in] argv Parameter-style argument list (key=value) followed by regular key arguments.
	@param[in] envp Environment variable list, NULL on end.

	@return 0 on success, a negative integer if a problem occurs, or a positive
		integer if the input file fails the L2QC check.
*/
int sstval_extract(sstval_extract_arguments *arguments);

/** @brief Clean up stuff malloc'd by the argpar callback.  Should be called at the end of processing
		if sstval_extract_argpar is used. */
int sstval_extract_clean(sstval_extract_arguments *arguments);

#endif /* ___SSTVAL_EXTRACT_H_ */
