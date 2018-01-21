// implemented version numbers
#define VALEXTRACT_IMPLEMENTATION 2005009
#define VALEXTRACT_IMPLEMENTATION_STR "2.5.9"
#define VALEXTRACT_API_VERSION 2005008
#define VALEXTRACT_API_VERSION_STR "2.5.8"

#include "argpar.h"
#include "pqueue.h"
#include "val_extract.h"

#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <libgen.h>
#include <limits.h>
#include <math.h>
#include <netcdf.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define str(s) #s
#define expanded_str(s) str(s)

static double strtod_strict(const char*); // should make a strict parser library

static const char *FLAG_MEANINGS_SEP = " ";

#define VALEXTRACT_OPTICS_DEFAULT 20

//const char *ignore_geospatial[] = { "l2_flags", NULL };
const char *ignore_geospatial[] = { NULL };

#define START_LINE_KEY  1
#define END_LINE_KEY    2
#define START_PIXEL_KEY 3
#define END_PIXEL_KEY   4
#define BOX_SIZE_KEY    5

#define SLAT_KEY 6
#define ELAT_KEY 7
#define SLON_KEY 8
#define ELON_KEY 9
#define L2QC_KEY 10
#define RANGES_KEY 11
#define IGNORE_KEY 12
#define OPTICS_KEY 13
#define COUNT_KEY 14

static argpar_option options[] = {
	{ "ifile", 'f', "FILE", 0, "input NetCDF file" },
	{ "optics_threshold", OPTICS_KEY, "N", OPTION_DBL | OPTION_HIDDEN, "OPTICS cluster threshold" },
	{ 0,0,0,0, "For validation box extracts, using pixels/lines:" },
	{ "sline", START_LINE_KEY, 0, OPTION_INT, "target line, zero-based (positive integer)" },
	{ "spixl", START_PIXEL_KEY, 0, OPTION_INT, "target pixel, zero-based (positive integer)" },
	{ "boxsize", BOX_SIZE_KEY, "N", OPTION_INT, "size of box (N by N pixels) (positive integer)" },
	{ 0,0,0,0, "For region extracts, using pixels/lines:" },
	{ "sline", START_LINE_KEY, 0, OPTION_INT, "starting line, zero-based (positive integer)" },
	{ "eline", END_LINE_KEY, 0, OPTION_INT, "ending line, zero-based (positive integer)" },
	{ "spixl", START_PIXEL_KEY, 0, OPTION_INT, "starting pixel, zero-based (positive integer)" },
	{ "epixl", END_PIXEL_KEY, 0, OPTION_INT, "ending pixel, zero-based (positive integer)" },
	{ 0,0,0,0, "For validation box extracts, using lat/lon:" },
	{ "slat", SLAT_KEY, 0, OPTION_DBL, "latitude target lies on" },
	{ "slon", SLON_KEY, 0, OPTION_DBL, "longitude target lies on" },
	{ "boxsize", BOX_SIZE_KEY, "N", OPTION_INT, "size of box (N by N pixels) (positive integer)" },
	{ 0,0,0,0, "For region extracts, using lat/lon:" },
	{ "slat", SLAT_KEY, 0, OPTION_DBL, "starting latitude (south-most boundary)" },
	{ "elat", ELAT_KEY, 0, OPTION_DBL, "ending latitude (north-most boundary)" },
	{ "slon", SLON_KEY, 0, OPTION_DBL, "starting longitude (west-most boundary)" },
	{ "elon", ELON_KEY, 0, OPTION_DBL, "ending longitude (east-most boundary)" },
	{ 0,0,0, OPTION_DOC, "If no location information is given, the entire file is processed." },
	{ 0,0,0,0, "Filtering options:" },
	{ "l2qc", L2QC_KEY, "FLAG=decimal ...", 0, "Max percentages (in decimal form) an extract can contain "
		"of a flag before failing the Level-2 Quality Checker (L2QC)."},
	{ "valid_ranges", RANGES_KEY, "var_prefix=min:max ...", 0, "Valid ranges for variables.  Any values "
		"outside these ranges will not be used for statistics."},
	{ "ignore_flags", IGNORE_KEY, "FLAG ...", 0, "Ignore values for ALL variables in pixels containing "
		"the flag(s) listed."},
	{ "count_flags", COUNT_KEY, "FLAG ...", 0, "Count the number of pixels containing ANY of the flag(s) "
			"listed."},
	{ 0,0,0, OPTION_DOC | OPTION_DOC_NO_BREAK, "Examples:" },
	{ 0,0,0, OPTION_DOC | OPTION_DOC_NO_BREAK, " l2qc=CLDICE=0.5 HIGLINT=0.25" },
	{ 0,0,0, OPTION_DOC | OPTION_DOC_NO_BREAK, " valid_ranges=Rrs_=0.001:0.02 chlor_a=0:10" },
	{ 0,0,0, OPTION_DOC, " ignore_flags=CLDICE HIGLINT" },
	{ 0,0,0, OPTION_DOC, "The flag names for the L2QC thresholds and to ignore them will be matched exactly "
		"and this will error if that flag is not found.  valid_ranges uses case-sensitive prefix matching.  "
		"All numbers are inclusive, e.g., 50% of pixels having the cloud flag passes QC in the "
		"above example and 10 is a valid chlor_a value." },
	{ 0,0,0, OPTION_DOC, "Failing the L2QC will successfully process the file, but will not process any "
		"variables within the file." },
	{ 0,0,0,0, "Caveats:" },
	{ 0,0,0, OPTION_DOC, "The \"center pixel\" for lat-/lon-based regions is only guaranteed to be inside the "
		"actual data, but isn't necessarily the actual center due to weirdly shaped granules and endless edge "
		"cases involving partial polygonal collision.  Line- and pixel-based should be more accurate, but "
		"YMMV.  The detected center is guaranteed to be consistent and repeatable." },
	{ 0 } // tell argpar to stop checking options
};

#ifdef DEBUG
#define dprintf(...) do { printf(__VA_ARGS__); } while(0)
#else
#define dprintf(...) do {} while (0)
#endif

#define eprintf(...) do { fprintf(stderr, __VA_ARGS__); } while(0)

static int parse_options(int key, char *argv, struct argpar_state *state) {
	int ret = 0;

	val_extract_arguments *arguments = state->input;
	// key > 0 means it has to be an actual option
	if (key >= START_LINE_KEY && key <= ELON_KEY && ((key < SLAT_KEY && state->argv_as_int_err) || (key >= SLAT_KEY && state->argv_as_dbl_err))){
		eprintf("Invalid cast for %s\n", state->arg_name);
		ret = ARGPAR_ERR_ABORT;
	} else {
		switch (key){
		case 'h':
			return ARGPAR_ERR_USAGE;
		case 'f':
			arguments->ifile = argv;
			break;
		case BOX_SIZE_KEY:
			arguments->box_size = state->argv_as_int;
			break;
		case START_LINE_KEY:
			arguments->start_line = state->argv_as_int;
			break;
		case END_LINE_KEY:
			arguments->end_line = state->argv_as_int;
			break;
		case START_PIXEL_KEY:
			arguments->start_pixel = state->argv_as_int;
			break;
		case END_PIXEL_KEY:
			arguments->end_pixel = state->argv_as_int;
			break;
		case SLAT_KEY:
			arguments->start_lat = state->argv_as_dbl;
			break;
		case ELAT_KEY:
			arguments->end_lat = state->argv_as_dbl;
			break;
		case SLON_KEY:
			arguments->start_lon = state->argv_as_dbl;
			break;
		case ELON_KEY:
			arguments->end_lon = state->argv_as_dbl;
			break;
		case L2QC_KEY:
			if (strlen(argv)){
				char tmp[strlen(argv)+1];
				strncpy(tmp, argv, strlen(argv));
				tmp[strlen(argv)] = '\0';
				char *tmp_p = argv;
				int space_count = 0;
				if (!isspace(*tmp_p)){
					++space_count;
				}
				while (*tmp_p){
					if (isspace(*tmp_p)){
						++space_count;
					}
					++tmp_p;
				}

				arguments->l2qc_thresholds = calloc(space_count, sizeof(const char*));
				if (!arguments->l2qc_thresholds){
					eprintf("Error allocating memory in '%s', near line %d\n", __FILE__, __LINE__);
					return ARGPAR_ERR_ABORT;
				}
				arguments->l2qc_flags = calloc(space_count, sizeof(const char*));
				if (!arguments->l2qc_flags){
					eprintf("Error allocating memory in '%s', near line %d\n", __FILE__, __LINE__);
					return ARGPAR_ERR_ABORT;
				}

				int i = 0;
				char *token = strtok(tmp, " ");
				while (token){
					if (strlen(token)){
						char *tmp_v = token;
						char *tmp_k = strsep(&tmp_v, "=");
						if (!tmp_v || !tmp_k){
							eprintf("Invalid L2QC spec: %s.\n", token);
							return ARGPAR_ERR_ABORT;
						}
						arguments->l2qc_flags[i] = malloc((strlen(tmp_k)+1) * sizeof(char));
						if (!arguments->l2qc_flags[i]){
							eprintf("Error allocating memory in '%s', near line %d\n", __FILE__, __LINE__);
							return ARGPAR_ERR_ABORT;
						}
						strcpy(arguments->l2qc_flags[i], tmp_k);
						arguments->l2qc_threshold_count++;
						arguments->l2qc_thresholds[i] = strtod_strict(tmp_v);
						if (errno){
							eprintf("Invalid L2QC spec: %s=%s.\n", tmp_k, tmp_v);
							return ARGPAR_ERR_ABORT;
						}
						if (arguments->l2qc_thresholds[i] < 0){
							eprintf("L2QC thresholds cannot be negative: %s=%s.\n", tmp_k, tmp_v);
							return ARGPAR_ERR_ABORT;
						}
						++i;
					}
					token = strtok(NULL, " ");
				}
			}
			break;
		case RANGES_KEY:
			if (strlen(argv)){
				char tmp[strlen(argv)+1];
				strncpy(tmp, argv, strlen(argv));
				tmp[strlen(argv)] = '\0';
				char *tmp_p = argv;
				int space_count = 0;
				if (!isspace(*tmp_p)){
					++space_count;
				}
				while (*tmp_p){
					if (isspace(*tmp_p)){
						++space_count;
					}
					++tmp_p;
				}

				arguments->valid_ranges = calloc(space_count, sizeof(val_extract_valid_range));
				if (!arguments->valid_ranges){
					eprintf("Error allocating memory in '%s', near line %d\n", __FILE__, __LINE__);
					return ARGPAR_ERR_ABORT;
				}

				int i = 0;
				char *token = strtok(tmp, " ");
				while (token){
					if (strlen(token)){
						char *tmp_v = token;
						char *tmp_k = strsep(&tmp_v, "=");
						if (!tmp_v || !tmp_k){
							eprintf("Invalid L2QC spec: %s.\n", token);
							return ARGPAR_ERR_ABORT;
						}
						arguments->valid_ranges[i].prefix = malloc((strlen(tmp_k)+1) * sizeof(char));
						if (!arguments->valid_ranges[i].prefix){
							eprintf("Error allocating memory in '%s', near line %d\n", __FILE__, __LINE__);
							return ARGPAR_ERR_ABORT;
						}
						strcpy(arguments->valid_ranges[i].prefix, tmp_k);
						arguments->valid_range_count++;

						char *tmp_max = tmp_v;
						char *tmp_min = strsep(&tmp_max, ":");

						arguments->valid_ranges[i].min = strtod_strict(tmp_min);
						int error = errno;
						arguments->valid_ranges[i].max = strtod_strict(tmp_max);
						error |= errno;
						if (error){
							eprintf("Invalid range spec: %s=%s:%s.\n", tmp_k, tmp_min, tmp_max);
							return ARGPAR_ERR_ABORT;
						}
						++i;
					}
					token = strtok(NULL, " ");
				}
			}
			break;
		case IGNORE_KEY:
			if (strlen(argv)){
				char tmp[strlen(argv)+1];
				strncpy(tmp, argv, strlen(argv));
				tmp[strlen(argv)] = '\0';
				char *tmp_p = argv;
				int space_count = 0;
				if (!isspace(*tmp_p)){
					++space_count;
				}
				while (*tmp_p){
					if (isspace(*tmp_p)){
						++space_count;
					}
					++tmp_p;
				}

				arguments->ignore_flags = calloc(space_count, sizeof(const char*));
				if (!arguments->ignore_flags){
					eprintf("Error allocating memory in '%s', near line %d\n", __FILE__, __LINE__);
					return ARGPAR_ERR_ABORT;
				}

				int i = 0;
				char *token = strtok(tmp, " ");
				while (token){
					if (strlen(token)){
						arguments->ignore_flags[i] = malloc((strlen(token)+1) * sizeof(char));
						if (!arguments->ignore_flags[i]){
							eprintf("Error allocating memory in '%s', near line %d\n", __FILE__, __LINE__);
							return ARGPAR_ERR_ABORT;
						}
						strcpy(arguments->ignore_flags[i], token);
						arguments->ignore_flag_count++;
						++i;
					}
					token = strtok(NULL, " ");
				}
			}
			break;
		case COUNT_KEY:
			if (strlen(argv)){
				char tmp[strlen(argv)+1];
				strncpy(tmp, argv, strlen(argv));
				tmp[strlen(argv)] = '\0';
				char *tmp_p = argv;
				int space_count = 0;
				if (!isspace(*tmp_p)){
					++space_count;
				}
				while (*tmp_p){
					if (isspace(*tmp_p)){
						++space_count;
					}
					++tmp_p;
				}

				arguments->count_flags = calloc(space_count, sizeof(const char*));
				if (!arguments->count_flags){
					eprintf("Error allocating memory in '%s', near line %d\n", __FILE__, __LINE__);
					return ARGPAR_ERR_ABORT;
				}

				int i = 0;
				char *token = strtok(tmp, " ");
				while (token){
					if (strlen(token)){
						arguments->count_flags[i] = malloc((strlen(token)+1) * sizeof(char));
						if (!arguments->count_flags[i]){
							eprintf("Error allocating memory in '%s', near line %d\n", __FILE__, __LINE__);
							return ARGPAR_ERR_ABORT;
						}
						strcpy(arguments->count_flags[i], token);
						arguments->count_flag_count++;
						++i;
					}
					token = strtok(NULL, " ");
				}
			}
			break;
		case ARGPAR_KEY_INIT:
			arguments->ifile = NULL;
			arguments->box_size = 0;
			arguments->optics_threshold = VALEXTRACT_OPTICS_DEFAULT;
			arguments->start_line = VALEXTRACT_UNSET;
			arguments->end_line = VALEXTRACT_UNSET;
			arguments->start_pixel = VALEXTRACT_UNSET;
			arguments->end_pixel = VALEXTRACT_UNSET;
			arguments->start_lat = VALEXTRACT_UNSET;
			arguments->end_lat = VALEXTRACT_UNSET;
			arguments->start_lon = VALEXTRACT_UNSET;
			arguments->end_lon = VALEXTRACT_UNSET;
			arguments->count_flag_count = 0;
			arguments->ignore_flag_count = 0;
			arguments->l2qc_threshold_count = 0;
			arguments->valid_range_count = 0;
			break;
		case ARGPAR_KEY_SUCCESS:
			if (arguments->ifile == NULL){
				eprintf("input file required\n");
				ret = ARGPAR_ERR_USAGE;
			} else if (arguments->box_size > 0){
				if (arguments->start_line == VALEXTRACT_UNSET && arguments->start_pixel == VALEXTRACT_UNSET && arguments->start_lat == VALEXTRACT_UNSET && arguments->start_lon == VALEXTRACT_UNSET){
					eprintf("sline and spixl or slat and slon are required if boxsize is provided\n");
					ret = ARGPAR_ERR_USAGE;
				} else if (arguments->start_lat != VALEXTRACT_UNSET && arguments->start_lon != VALEXTRACT_UNSET){
					arguments->lat_and_lon = true;
					if (arguments->start_lat < -90 || arguments->start_lat > 90 || arguments->start_lon < -180 || arguments->start_lon > 180){
						eprintf("slat or slon is out of bounds\n");
						ret = ARGPAR_ERR_USAGE;
					}
				} else if (arguments->start_line < 0 || arguments->start_pixel < 0){
					eprintf("sline and spixl are both required and must be positive integers\n");
					ret = ARGPAR_ERR_USAGE;
				} else {
					arguments->line_and_pixel = true;
				}
			} else if (arguments->box_size < 0){
				eprintf("boxsize must be a positive integer\n");
			} else {
				if (arguments->start_lat != VALEXTRACT_UNSET || arguments->end_lat != VALEXTRACT_UNSET || arguments->start_lon != VALEXTRACT_UNSET || arguments->end_lon != VALEXTRACT_UNSET){
					if (arguments->start_lat < -90 || arguments->start_lat > 90 || arguments->start_lon < -180 || arguments->start_lon > 180){
						eprintf("slat or slon is out of bounds or was unspecified\n");
						ret = ARGPAR_ERR_USAGE;
					} else if (arguments->end_lat < -90 || arguments->end_lat > 90 || arguments->end_lon < -180 || arguments->end_lon > 180){
						eprintf("elat or elon is out of bounds or was unspecified\n");
						ret = ARGPAR_ERR_USAGE;
					} else if (arguments->start_lat > arguments->end_lat){
						eprintf("slat must be less than elat\n");
						ret = ARGPAR_ERR_USAGE;
					} else if (arguments->start_lon > arguments->end_lon){
						eprintf("slon must be less than elon\n");
						ret = ARGPAR_ERR_USAGE;
					}
					arguments->lat_and_lon = true;
				} else if (arguments->start_line != VALEXTRACT_UNSET || arguments->end_line != VALEXTRACT_UNSET || arguments->start_pixel != VALEXTRACT_UNSET || arguments->end_pixel != VALEXTRACT_UNSET){
					if (arguments->start_line < 0 || arguments->end_line < 0 || arguments->start_pixel < 0 || arguments->end_pixel < 0){
						eprintf("sline, eline, spixl, and epixl are all required and must be positive integers\n");
						ret = ARGPAR_ERR_USAGE;
					} else if (arguments->start_line > arguments->end_line){
						eprintf("sline must be less than eline\n");
						ret = ARGPAR_ERR_USAGE;
					} else if (arguments->start_pixel > arguments->end_pixel){
						eprintf("spixl must be less than epixl\n");
						ret = ARGPAR_ERR_USAGE;
					}
					arguments->line_and_pixel = true;
				} else {
					arguments->start_lat = -90;
					arguments->end_lat = 90;
					arguments->start_lon = -180;
					arguments->end_lon = 180;
					arguments->lat_and_lon = true;
				}
			}
			break;
		}
	}
	return ret;
}

argpar val_extract_argpar = { options, parse_options };

static void print_time();
static bool starts_with(const char *str, const char *prefix);
static void calc_dim_multipliers(int ndims, const int *dims, int *dim_multipliers);
static int flatten_index(int ndims, const int *dim_multipliers, int *index);
static int convert_lat_lons(val_extract_arguments *arguments, nc_region *region);
static int process_variable(val_extract_arguments *arguments, nc_var *var);
static double process_scan_line_att(val_extract_arguments *arguments, nc_region *region, const char *name);
static double process_location_var(val_extract_arguments *arguments, nc_region *region, const char *name);
static int process_file_attrs(val_extract_arguments *arguments, nc_region *region);
static int get_var_attributes(val_extract_arguments *arguments, nc_var *var);
static void free_var(val_extract_arguments *arguments, nc_var *var);
static nc_var_stats generate_stats(const double *data, const int offset, const int length);
void get_array_ranges(const double *sorted_data, const int length, const double start, const double end, int *start_i, int *end_i);

#define clean_and_return(ret) do { \
	if (ret == VALEXTRACT_ERR_NONE){ \
		arguments->val_extract_parser(VALEXTRACT_KEY_SUCCESS, &region, arguments->user_input); \
	} else { \
		arguments->val_extract_parser(VALEXTRACT_KEY_ERROR, &region, arguments->user_input); \
	} \
	nc_close(ncid); \
	arguments->val_extract_parser(VALEXTRACT_KEY_FINI, &region, arguments->user_input); \
	if (all_vars){ \
		for (product_count--; product_count >= 0; product_count--){ \
			free_var(arguments, all_vars[product_count]); \
			free(all_vars[product_count]); \
		} \
	} \
	free(all_vars); \
	if (region.boxes){ \
		free(region.boxes); \
	} \
	if (region.flag_counts){ \
		free(region.flag_counts); \
	} \
	if (region.pixel_flags){ \
		free(region.pixel_flags); \
	} \
	if (file.flag_masks){ \
		free(file.flag_masks); \
	} \
	if (file.flag_string){ \
		free(file.flag_string); \
	} \
	if (file.platform){ \
		free(file.platform); \
	} \
	if (file.instrument){ \
		free(file.instrument); \
	} \
	if (file.flag_meanings){ \
		free(file.flag_meanings); \
	} \
	return ret; \
} while (0)

int val_extract(val_extract_arguments *arguments) {
	if (!arguments || !arguments->val_extract_parser){
		return VALEXTRACT_ERR_INPUT;
	}

	arguments->val_extract_parser(VALEXTRACT_KEY_INIT, NULL, arguments->user_input);

	dprintf("ifile: %s\n", arguments->ifile);

	int ncid, nc_ret, i;
	if ((nc_ret = nc_open(arguments->ifile, NC_NOWRITE, &ncid)) != NC_NOERR){
		eprintf("Error opening file, error %d\n", nc_ret);
		return VALEXTRACT_ERR_NCFILE_ERR;
	}

	int ndims, nvars, ngatts, unlimdimid;
	if ((nc_ret = nc_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid)) != NC_NOERR){
		eprintf("Error getting file information, error %d\n", nc_ret);
		return VALEXTRACT_ERR_NCFILE_ERR;
	}
//	dprintf("nc_inq: %d | %d | %d | %d\n", ndims, nvars, ngatts, unlimdimid);

	char dim_name[64];
	int dim_lengths[ndims];
	int line_count = -1, pixel_count = -1;
	int line_dimid = -1, pixel_dimid = -1, pixel_control_dimid = -1;
	for (i=0;i<ndims;i++){
		size_t dim_length;
		if ((nc_ret = nc_inq_dim(ncid, i, dim_name, &dim_length)) != NC_NOERR){
			eprintf("Error getting dimensions, error %d\n", nc_ret);
			return VALEXTRACT_ERR_NCFILE_ERR;
		}
		dprintf("nc_inq_dim: %d | %zd | %s\n", i, dim_length, dim_name);
		dim_lengths[i] = dim_length;
		if ((strcmp(dim_name, "number_of_lines") == 0) || (strcmp(dim_name, "Number_of_Scan_Lines") == 0))   {
			line_count = dim_length;
			line_dimid = i;
		} else if ((strcmp(dim_name, "pixels_per_line") == 0) || (strcmp(dim_name, "Pixels_per_Scan_Line") == 0)) {
			pixel_count = dim_length;
			pixel_dimid = i;
		} else if ((strcmp(dim_name, "pixel_control_points") == 0)) {
			pixel_control_dimid = i;
		}
	}

	dprintf("line dim: %d, pixel dim: %d\n", line_count, pixel_count);
	if (line_count < 0 || pixel_count < 0){
		eprintf("Could not identify line/pixel dims, exiting\n");
		return VALEXTRACT_ERR_NCFILE_INVALID;
	}

	int ngrps;
	if ((nc_ret = nc_inq_grps(ncid, &ngrps, NULL)) != NC_NOERR){
		eprintf("Error finding group count, error %d\n", nc_ret);
		return VALEXTRACT_ERR_NCFILE_ERR;
	}
	int gids[ngrps];
	if ((nc_ret = nc_inq_grps(ncid, &ngrps, gids)) != NC_NOERR){
		eprintf("Error finding gids, error %d\n", nc_ret);
		return VALEXTRACT_ERR_NCFILE_ERR;
	}

	nc_file file = {
		.file_path = arguments->ifile, .ncid = ncid,
		.ndims = ndims, .dim_lengths = dim_lengths,
		.ngrps = ngrps, .gids = gids,
		.line_dimid = line_dimid, .pixel_dimid = pixel_dimid, .pixel_control_dimid = pixel_control_dimid,
		.control_point_size = dim_lengths[pixel_control_dimid] / dim_lengths[pixel_dimid]
	};
	nc_region region = {
		.file = &file,
		.box_count = 0
	};
	int product_count = 0;
	nc_var **all_vars = NULL;

	if (arguments->lat_and_lon){
		int convert_err = convert_lat_lons(arguments, &region);
		if (convert_err){
			if (convert_err != VALEXTRACT_ERR_POINT_NOT_FOUND){
				convert_err = VALEXTRACT_ERR_NCFILE_INVALID;
			}
			clean_and_return(convert_err);
		}
	} else {
		int start_line, end_line, start_pixel, end_pixel;
		if (arguments->box_size <= 0) {
			nc_point center_pixel = {
				.line = (arguments->start_line + arguments->end_line) >> 1,
				.pixel = (arguments->start_pixel + arguments->end_pixel) >> 1
			};
			region.center = center_pixel;

			start_line = arguments->start_line;
			end_line = arguments->end_line;
			start_pixel = arguments->start_pixel;
			end_pixel = arguments->end_pixel;
		} else {
			nc_point center_pixel = {
				.line = arguments->start_line,
				.pixel = arguments->start_pixel
			};
			region.center = center_pixel;

			int geo_box_margin = arguments->box_size >> 1; // if odd, -1, then divide by 2; else divide by 2
			start_line = arguments->start_line - geo_box_margin;
			end_line = arguments->start_line + geo_box_margin;
			start_pixel = arguments->start_pixel - (geo_box_margin * file.control_point_size);
			end_pixel = arguments->start_pixel + (geo_box_margin * file.control_point_size);
		}

		if (start_line < 0){
			start_line = 0;
		}
		if (start_pixel < 0){
			start_pixel = 0;
		}
		if (end_line >= line_count){
			end_line = line_count-1;
		}
		if (end_pixel >= pixel_count){
			end_pixel = pixel_count-1;
		}

		region.box_count = 1;
		region.boxes = malloc(sizeof(nc_box));
		if (region.boxes == NULL){
			eprintf("Error allocating memory in '%s', near line %d\n", __FILE__, __LINE__);
			clean_and_return(VALEXTRACT_ERR_UNKNOWN);
		}
		region.boxes[0].start[0] = start_line;
		region.boxes[0].start[1] = start_pixel;
		region.boxes[0].count[0] = end_line - start_line + 1;
		region.boxes[0].count[1] = end_pixel - start_pixel + 1;
	}

	if (region.box_count == 0){
		eprintf("Point not found in input file\n");
		clean_and_return(VALEXTRACT_ERR_NCFILE_INVALID);
	}

	if ((nc_ret = process_file_attrs(arguments, &region)) != 0){
		if (nc_ret == VALEXTRACT_CMD_STOP){
			nc_ret = 0;
		}
		clean_and_return(nc_ret);
	}

	if (arguments->product_count){
		product_count = arguments->product_count;
	} else {
		int gid_i;
		for (gid_i=0;gid_i<ngrps;gid_i++){
			int g_nvars;
			if ((nc_ret = nc_inq_varids(gids[gid_i], &g_nvars, NULL)) != NC_NOERR){
				eprintf("Error finding variable count, error %d\n", nc_ret);
				clean_and_return(VALEXTRACT_ERR_VARIABLE);
			}
			product_count += g_nvars;
		}
	}
	all_vars = malloc(product_count * sizeof(nc_var));
	product_count = 0;

	if (arguments->product_count){
		for (i=0;i<arguments->product_count;i++){
			int varid, gid, gid_i;
			varid = gid = -1;
			for (gid_i=0;gid_i<ngrps;gid_i++){
				if ((nc_ret = nc_inq_varid(gids[gid_i], arguments->products[i], &varid)) == NC_NOERR){
					gid = gids[gid_i];
					break;
				}
			}
			if (gid != -1 && varid != -1){
				dprintf("%s, gid = %d, varid = %d, nc_ret = %d\n", arguments->products[i], gid, varid, nc_ret);
				nc_var *var = malloc(sizeof(nc_var));
				var->file = &file;
				var->region = &region;
				var->name = arguments->products[i];
				var->gid = gid;
				var->varid = varid;
				if (get_var_attributes(arguments, var)){
					clean_and_return(VALEXTRACT_ERR_VARIABLE);
				} else {
					nc_ret = process_variable(arguments, var);
					if (nc_ret){
						if (nc_ret == VALEXTRACT_CMD_STOP){
							nc_ret = 0;
						}
						clean_and_return(VALEXTRACT_ERR_VARIABLE);
					}
					all_vars[product_count++] = var;
				}
			} else {
				eprintf("Variable %s not found\n", arguments->products[i]);
				clean_and_return(VALEXTRACT_ERR_VARIABLE);
			}
		}
	} else {
		int gid_i;
		for (gid_i=0;gid_i<ngrps;gid_i++){
			int g_nvars;
			if ((nc_ret = nc_inq_varids(gids[gid_i], &g_nvars, NULL)) != NC_NOERR){
				eprintf("Error finding variable count, error %d\n", nc_ret);
				clean_and_return(VALEXTRACT_ERR_VARIABLE);
			}
			int varids[g_nvars];
			if ((nc_ret = nc_inq_varids(gids[gid_i], &g_nvars, varids)) != NC_NOERR){
				eprintf("Error finding variables, error %d\n", nc_ret);
				clean_and_return(VALEXTRACT_ERR_VARIABLE);
			}
			dprintf("gid %d, variable count: %d\n", gids[gid_i], g_nvars);
			for (i=0;i<g_nvars;i++){
				nc_var *var = malloc(sizeof(nc_var));
				var->file = &file;
				var->region = &region;
				var->gid = gids[gid_i];
				var->varid = varids[i];
				if (get_var_attributes(arguments, var)){
					clean_and_return(VALEXTRACT_ERR_VARIABLE);
				} else {
					nc_ret = process_variable(arguments, var);
					if (nc_ret){
						if (nc_ret == VALEXTRACT_CMD_STOP){
							nc_ret = 0;
						}
						clean_and_return(VALEXTRACT_ERR_VARIABLE);
					}
					all_vars[product_count++] = var;
				}
			}
		}
	}

	dprintf("Exiting... ");
	print_time();
	clean_and_return(VALEXTRACT_ERR_NONE);
}

static uint32_t find_closest_index(double target_lat, double target_lon, uint32_t length, double *lats, double *lons){
	uint32_t i = 0, closest_i = -1;
	double best_diff = DBL_MAX;
	for (;i<length;i++){
		double diff = sqrt(pow(lats[i] - target_lat, 2) + pow(lons[i] - target_lon, 2));
		if (diff < best_diff){
			closest_i = i;
			best_diff = diff;
		}
	}
//	dprintf("Closest diff: %f, index %d\n", best_diff, closest_i);
	return closest_i;
}

#if 0
const long double pi = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899L;
static double rad2deg(const double rad){
      return (180.0 * rad / (pi));
}
static double deg2rad(const double deg){
      return (pi * deg / 180.0);
}
static long double rad2degl(const long double rad){
      return (180.0 * rad / (pi));
}
static long double deg2radl(const long double deg){
      return (pi * deg / 180.0);
}

static uint32_t find_closest_index_polar(double target_lat, double target_lon, uint32_t length, double *lats, double *lons){
	double target_x = cos(deg2rad(target_lat)) * cos(deg2rad(target_lon));
	double target_y = cos(deg2rad(target_lat)) * sin(deg2rad(target_lon));
	double target_z = sin(deg2rad(target_lat));

	uint32_t i = 0, closest_i = -1;
	double best_diff = DBL_MAX;
	for (;i<length;i++){
		double this_x = cos(deg2rad(lats[i])) * cos(deg2rad(lons[i]));
		double this_y = cos(deg2rad(lats[i])) * sin(deg2rad(lons[i]));
		double this_z = sin(deg2rad(lats[i]));

		double diff = this_x * target_x + this_y * target_y + this_z * target_z;
		if (diff > best_diff){
			closest_i = i;
			best_diff = diff;
		}
	}
	return closest_i;
}
#endif

static nc_box *add_to_boxes(nc_box *boxes, unsigned *box_shift, int *box_i, int *box_max, int line_i, int start_pixel, int end_pixel){
	if (*box_i == *box_max){
		nc_box *new_boxes_ptr = realloc(boxes, sizeof(nc_box) * (1 << (++*box_shift)));
		if (new_boxes_ptr == NULL){
			eprintf("Error allocating memory in '%s', near line %d\n", __FILE__, __LINE__);
			return NULL;
		}
		boxes = new_boxes_ptr;
		*box_max = (1 << *box_shift) - 1;
	}
	nc_box new_box = {.start = {line_i, start_pixel}, .count = {1, end_pixel - start_pixel + 1}};
	boxes[(*box_i)++] = new_box;
	return boxes;
}

// http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
static bool polygon_contains(double *poly_x, double *poly_y, const int poly_n, double x, double y){
	int i, j, c = 0;
	for (i = 0, j = poly_n-1; i < poly_n; j = i++) {
		if (((poly_x[i]>x) != (poly_x[j]>x)) && (y < (poly_y[j]-poly_y[i]) * (x-poly_x[i]) / (poly_x[j]-poly_x[i]) + poly_y[i])){
			c = !c;
		}
	}
	return c;
}
#if 0
// https://en.wikipedia.org/wiki/Centroid#Centroid_of_polygon
static void polygon_center(double *poly_x, double *poly_y, const int poly_n, double *x, double *y){
	if (poly_n == 0){
		*x = *y = 0;
	} else if (poly_n == 1){
		*x = poly_x[0];
		*y = poly_y[0];
	} else if (poly_n == 2){
		*x = (poly_x[0] + poly_x[1]) / 2;
		*y = (poly_y[0] + poly_y[1]) / 2;
	} else {
		int i;
		double area = 0;
		*x = *y = 0;
		for (i=0;i<poly_n;i++){
			double diff = (poly_x[i] * poly_y[i+1]) - (poly_x[i+1] * poly_y[i]);
			*x += (poly_x[i] + poly_x[i+1]) * diff;
			*y += (poly_y[i] + poly_y[i+1]) * diff;
			area += diff;
		}
		if (area){
			area /= 2;
			*x = *x / (6 * area);
			*y = *y / (6 * area);
		}
	}
}
#endif

#if 0
// DBSCAN is a relatively simple algorithm to detect clusters of points.  It's no longer used here due to requiring a density
//   factor that would either need a look up table or a bunch of extra analysis to find it first.  The oddities documented
//   for OPTICS below regarding neighbors apply to this algorithm, as well.
struct dbscan {
	double x, y;
	bool noise, visited;
	int cluster;
	int neighbors[8], neighbor_count;
};


static void dbscan_expand(struct dbscan *db, struct dbscan *p, const int cluster, const double epsilon, const int min_points){
	p->cluster = cluster;
	int neighbor_i;
	for (neighbor_i=0;neighbor_i<p->neighbor_count;neighbor_i++){
		struct dbscan *neighbor = &db[p->neighbors[neighbor_i]];
		if (!neighbor->visited){
			neighbor->visited = true;
			if (neighbor->neighbor_count >= min_points){
				dbscan_expand(db, neighbor, cluster, epsilon, min_points);
			}
		}
		if (!neighbor->cluster){
			neighbor->cluster = cluster;
		}
	}
}

static void dbscan(struct dbscan *db, const int db_count, const double epsilon, const int min_points){
	int db_i, cluster = 0;
	for (db_i=0;db_i<db_count;db_i++){
		db[db_i].noise = false;
		db[db_i].visited = false;
		db[db_i].cluster = 0;
	}
	for (db_i=0;db_i<db_count;db_i++){
		if (!db[db_i].visited){
			db[db_i].visited = true;
			if (db[db_i].neighbor_count < min_points){
				db[db_i].noise = true;
			} else {
				dbscan_expand(db, &db[db_i], ++cluster, epsilon, min_points);
			}
		}
	}
}
#endif

//https://en.wikipedia.org/wiki/OPTICS_algorithm
// OPTICS is an algorithm to identify clusters of points in 2D space.  The advantage of it is that it's
//   insensitive of cluster density, which is why it was chosen (pixel resolution = density).  The reason
//   it's needed at all is because of gaps in data, such as when the SeaWiFS sensor switches tilt.  When
//   that happens, a naive pixel grab could group pixels thousands of kilometers apart.
//
// Some weirdness that should be documented:
//   Epsilon (eps, in Wikipedia's pseudo-code) isn't even used.  I don't even use a getNeighbors function.
//   I already know the neighbors because it's just a grid, so each point just has all of its gridly
//   neighbors in its neighbor list.  This is slightly incorrect because I'm using this to detect gaps and the
//   appropriate epsilon would simply discard them as neighbors.  This entire bit is used specifically because
//   I can't (easily) know the appropriate epsilon; one might be able to use sensor lookup tables and geospatial
//   geometry to calculate the correct grid size and, therefore, epsilon.  To keep this library as generic as
//   possible, I use OPTICS's density-independent concept and a maximum ratio to determine if a pixel is too
//   far away from the core and should be considered a new cluster.  If a reachability-distance is X times
//   larger than the current cluster's moving average, it's a new cluster.   For these reasons, and any
//   assumptions I've made for this specific application, the struct for the data point is a kind of weird
//   format (for convenience), so this algorithm isn't in a nice, generic library like it should be (nor do I
//   plan to move it to its own library).
struct optic {
	double x, y;
	bool processed;
	double reachability, core;
	int neighbors[8], neighbor_count;
	unsigned i,j; // not part of the algorithm
};

static double optics_core_distance(const struct optic *db, const struct optic *p, const int min_points){
	if (p->neighbor_count < min_points || p->i == UINT_MAX){
		return DBL_MIN;
	}
	double distances[p->neighbor_count];

	int i;
	for (i=0;i<p->neighbor_count;i++){
		const struct optic *neighbor = &db[p->neighbors[i]];
		distances[i] = sqrt(pow(p->x - neighbor->x, 2) + pow(p->y - neighbor->y, 2));
	}
	qsort(distances, p->neighbor_count, sizeof(double), compare_double);
	if (min_points <= p->neighbor_count){
		return distances[min_points - 1];
	} else {
		return distances[p->neighbor_count - 1];
	}
}
static void optics_update(struct optic *db, struct optic *p, const int min_points, pqueue *seeds){
	int i;
	for (i=0;i<p->neighbor_count;i++){
		struct optic *neighbor = &db[p->neighbors[i]];
		if (!neighbor->processed){
			double dist_to_neighbor = sqrt(pow(p->x - neighbor->x, 2) + pow(p->y - neighbor->y, 2));
			double reachability = (dist_to_neighbor > p->core ? dist_to_neighbor : p->core);
			if (neighbor->reachability == DBL_MIN){
				pqueue_push(seeds, neighbor->reachability = reachability, neighbor);
			} else if (reachability < neighbor->reachability) {
				pqueue_repush(seeds, neighbor->reachability = reachability, neighbor);
			}
		}
	}
}
static void optics(struct optic *db, const int db_count, const int min_points, struct optic **out_db){
	int db_i;
	for (db_i=0;db_i<db_count;db_i++){
		db[db_i].processed = false;
		db[db_i].reachability = db[db_i].core = DBL_MIN;
	}
	pqueue *seeds = pqueue_create();
	int out_db_i = 0;
	for (db_i=0;db_i<db_count;db_i++){
		if (!db[db_i].processed){
			db[db_i].processed = true;
			out_db[out_db_i++] = &db[db_i];
			if ((db[db_i].core = optics_core_distance(db, &db[db_i], min_points)) != DBL_MIN){
				optics_update(db, &db[db_i], min_points, seeds);
				struct optic *seed = NULL;
				while ((seed = pqueue_pull(seeds)) != NULL){
					seed->processed = true;
					out_db[out_db_i++] = seed;
					if ((db[db_i].core = optics_core_distance(db, seed, min_points)) != DBL_MIN){
						optics_update(db, seed, min_points, seeds);
					}
				}
			}
		}
	}
	pqueue_destroy(seeds);
}

// This assumes the latitude and longitude variables do not have any scaling or offset and that the
// fill value is -999.  The fill value is currently listed incorrectly in the data files, so it's a
// moot point to do it correctly here.  I only check the latitude array for -999; if one is -999,
// then they both are.
static int convert_lat_lons(val_extract_arguments *arguments, nc_region *region){
	nc_file *file = region->file;
	int latlon_length = file->dim_lengths[file->line_dimid] * file->dim_lengths[file->pixel_dimid];
	double *lats = malloc(sizeof(double) * latlon_length);
	if (lats == NULL){
		eprintf("Error allocating memory in '%s', near line %d\n", __FILE__, __LINE__);
		return -1;
	}
	double *lons = malloc(sizeof(double) * latlon_length);
	if (lons == NULL){
		eprintf("Error allocating memory in '%s', near line %d\n", __FILE__, __LINE__);
		free(lats);
		return -1;
	}
	bool found_lat, found_lon;
	found_lat = found_lon = false;
	int nc_ret, varid, gid_i;
	varid = -1;
	for (gid_i=0;gid_i<file->ngrps && !(found_lat && found_lon);gid_i++){
		if (!found_lat && (nc_ret = nc_inq_varid(file->gids[gid_i], "latitude", &varid)) == NC_NOERR){
			nc_get_var_double(file->gids[gid_i], varid, lats);
			found_lat = true;
		}
		if (!found_lon && (nc_ret = nc_inq_varid(file->gids[gid_i], "longitude", &varid)) == NC_NOERR){
			nc_get_var_double(file->gids[gid_i], varid, lons);
			found_lon = true;
		}
	}
	if (!found_lat || !found_lon){
		eprintf("Failed finding latitude and/or longitude variable");
		free(lats);
		free(lons);
		return -1;
	}

	if (arguments->box_size <= 0){
		unsigned box_shift = 3;
		nc_box *boxes = malloc(sizeof(nc_box) << box_shift);
		if (boxes == NULL){
			eprintf("Error allocating memory in '%s', near line %d\n", __FILE__, __LINE__);
			free(lats);
			free(lons);
			return -1;
		}
		int box_max = (1 << box_shift) - 1;
		int box_i = 0;

		// NW NE SE SW NW
		double poly_lats[] = {
			arguments->end_lat, arguments->end_lat,
			arguments->start_lat, arguments->start_lat, arguments->end_lat
		};
		double poly_lons[] = {
			arguments->start_lon, arguments->end_lon,
			arguments->end_lon, arguments->start_lon, arguments->start_lon,
		};

		// make one-line wide strips of pixels (I call 'em boxes)
		// It might be worth pursuing finding the least amount of rectangles, which may significantly
		//     reduce the calls to nc_get_vara_<type>.  But, I'm lazy and uneducated.
		int pixel_i, line_i, latlon_i = 0;
		int start_pixel, end_pixel;
		for (line_i=0;line_i<file->dim_lengths[file->line_dimid];line_i++){
			start_pixel = end_pixel = -1;
			for (pixel_i=0;pixel_i<file->dim_lengths[file->pixel_dimid];pixel_i++, latlon_i++){
				if (lats[latlon_i] != -999 && polygon_contains(poly_lats, poly_lons, 5, lats[latlon_i], lons[latlon_i])){
					if (start_pixel == -1){
						start_pixel = end_pixel = pixel_i;
					} else {
						end_pixel = pixel_i;
					}
				} else if (start_pixel != -1){
					nc_box *new_box_ptr = add_to_boxes(boxes, &box_shift, &box_i, &box_max, line_i, start_pixel, end_pixel);
					if (new_box_ptr == NULL){
						eprintf("Error allocating memory in '%s', near line %d\n", __FILE__, __LINE__);
						free(boxes);
						free(lats);
						free(lons);
						return -1;
					}
					boxes = new_box_ptr;
					start_pixel = end_pixel = -1;
				}
			}
			if (start_pixel != -1){
				nc_box *new_box_ptr = add_to_boxes(boxes, &box_shift, &box_i, &box_max, line_i, start_pixel, end_pixel);
				if (new_box_ptr == NULL){
					eprintf("Error allocating memory in '%s', near line %d\n", __FILE__, __LINE__);
					free(boxes);
					free(lats);
					free(lons);
					return -1;
				}
				boxes = new_box_ptr;
			}
		}
		region->box_count = box_i;
		region->boxes = boxes;

		if (box_i == 0){
			eprintf("Region not found in input file\n");
			free(region->boxes);
			free(lats);
			free(lons);
			return -1;
		}
		nc_box center_box = region->boxes[box_i / 2];
		nc_point center_pixel = {
			.line = center_box.start[0] + (center_box.count[0] / 2),
			.pixel = center_box.start[1] + (center_box.count[1] / 2)
		};
		region->center = center_pixel;
	} else { // boxsize > 0
		int i;
		const int dims[2] = {file->dim_lengths[file->line_dimid], file->dim_lengths[file->pixel_dimid]};
		int dim_multipliers[2];
		calc_dim_multipliers(2, dims, dim_multipliers);

		int closest_index, result[2];
		unflatten_index(closest_index = find_closest_index(arguments->start_lat, arguments->start_lon, latlon_length, lats, lons), 2, dims, result);
		arguments->start_line = result[0];
		arguments->start_pixel = result[1];
		nc_point center_pixel = { .line = result[0], .pixel = result[1] };
		region->center = center_pixel;
		double distance_to_closest = sqrt(pow(lats[closest_index] - arguments->start_lat, 2) + pow(lons[closest_index] - arguments->start_lon, 2));

		dprintf("Closest pixel: %d (%d, %d), %f degrees away\n", closest_index, result[0], result[1], distance_to_closest);
		unsigned half_box = arguments->box_size >> 1;
		int line_start = result[0]-half_box;
		int line_end = result[0]+half_box;
		int pixel_start = result[1]-half_box;
		int pixel_end = result[1]+half_box;
		if (line_start < 0){
			line_start = 0;
		}
		if (line_end >= dims[0]){
			line_end = dims[0]-1;
		}
		if (pixel_start < 0){
			pixel_start = 0;
		}
		if (pixel_end >= dims[1]){
			pixel_end = dims[1]-1;
		}
		const int pixel_count = pixel_end-pixel_start+1;
		const int line_count = line_end-line_start+1;
		double **grid_lats = malloc(sizeof(double*) * line_count);
		double **grid_lons = malloc(sizeof(double*) * line_count);
		int **grid_cluster = malloc(sizeof(int*) * line_count);
		for (i=0;i<line_count;i++){
			grid_lats[i] = malloc(sizeof(double) * pixel_count);
			grid_lons[i] = malloc(sizeof(double) * pixel_count);
			grid_cluster[i] = malloc(sizeof(int) * pixel_count);
		}

		double dateline_offset = copysign(360, lons[closest_index]);
		int j;
		for (i=0;i<line_count;i++){
			int coords[2] = {line_start+i, pixel_start};
			int index = flatten_index(2, dim_multipliers, coords);
			for (j=0;j<pixel_count;j++){
				grid_lats[i][j] = lats[index + j];
				grid_lons[i][j] = lons[index + j];
				grid_cluster[i][j] = -1;
				if (grid_lats[i][j] != -999 && copysign(1, grid_lons[i][j]) != copysign(1, lons[closest_index])){
					if (fabs(grid_lons[i][j] - lons[closest_index]) > fabs(grid_lons[i][j] + dateline_offset - lons[closest_index])){
						grid_lons[i][j] += dateline_offset;
					}
				}
			}
		}

		int closest_line = -1, closest_pixel = -1;
//		struct optic db[line_count * pixel_count], *out_db[line_count * pixel_count];
		struct optic *db = malloc(sizeof(struct optic) * line_count * pixel_count);
		struct optic **out_db = malloc(sizeof(struct optic*) * line_count * pixel_count);
		int db_i = 0;
		for (i=0;i<line_count;i++){
			for (j=0;j<pixel_count;j++, db_i++){
				if (grid_lats[i][j] != -999){
					db[db_i].i = i;
					db[db_i].j = j;
					db[db_i].x = grid_lats[i][j];
					db[db_i].y = grid_lons[i][j];
					if (grid_lats[i][j] == lats[closest_index] && grid_lons[i][j] == lons[closest_index]){
						closest_line = i;
						closest_pixel = j;
					}
					int neighbor_i = 0;
					int n_i, n_j;
					for (n_i=i > 0 ? i-1 : i;n_i <= (i+1) && n_i<line_count;n_i++){
						for (n_j=j > 0 ? j-1 : j;n_j <= (j+1) && n_j<pixel_count;n_j++){
							if (grid_lats[n_i][n_j] != -999 && (i != n_i || j != n_j)){
								db[db_i].neighbors[neighbor_i++] = (n_i * pixel_count) + n_j;
							}
						}
					}
					db[db_i].neighbor_count = neighbor_i;
				} else {
					db[db_i].i = UINT_MAX;
				}
			}
		}

		if (closest_line == -1 || closest_pixel == -1){
			eprintf("Error finding closest pixel, somehow\n");
			for (i=0;i<line_count;i++){
				free(grid_lats[i]);
				free(grid_lons[i]);
				free(grid_cluster[i]);
			}
			free(grid_lats);
			free(grid_lons);
			free(grid_cluster);
			free(lats);
			free(lons);
			free(db);
			free(out_db);
			return -1;
		}

		optics(db, db_i, 3, out_db);

		double current_total = out_db[0]->core;
		int current_count = 0;

		db_i = 0;
		int cluster = 0;
		for (i=0;i<line_count * pixel_count;i++, db_i++){
			if (out_db[db_i]->i != UINT_MAX){
				if (current_count && current_total && (out_db[db_i]->reachability / (current_total/current_count)) > arguments->optics_threshold){
					dprintf("Gap detected, cluster %d\n", cluster);
					dprintf("\tcount %d, total %f, threshold %f > %f\n", current_count, current_total, (out_db[db_i]->reachability / (current_total/current_count)), arguments->optics_threshold);
					dprintf("\t(%f / (%f/%d)) > %f\n", out_db[db_i]->reachability, current_total, current_count, arguments->optics_threshold);
					cluster++;
					// Using the actual reachability would desensitize the detection because it would start off with too large of a total.
					// Using the previous average assumes that clusters have similar density.  Assuming anything isn't ideal, but this
					// should hold up pretty well.
					current_total = current_total / current_count;
					current_count = 1;
				} else {
					current_total += out_db[db_i]->reachability;
					current_count++;
				}
				grid_cluster[out_db[db_i]->i][out_db[db_i]->j] = cluster;
			}
		}
		if (cluster){
			dprintf("Gap detected, last cluster (%d)\n", cluster);
			dprintf("\tcount %d, total %f, threshold %f <= %f\n", current_count, current_total, (current_total/current_count), arguments->optics_threshold);
			dprintf("\t(%f/%d) <= %f\n", current_total, current_count, arguments->optics_threshold);
		}
		free(db);
		free(out_db);

		cluster = grid_cluster[closest_line][closest_pixel];

		double max_neighbor_distance = -1;
		for (i=0;i<line_count;i++){
			for (j=0;j<pixel_count;j++){
				if (grid_lats[i][j] != -999 && grid_cluster[i][j] == cluster){
					int n_i, n_j;
					for (n_i=i > 0 ? i-1 : i;n_i <= (i+1) && n_i<line_count;n_i++){
						for (n_j=j > 0 ? j-1 : j;n_j <= (j+1) && n_j<pixel_count;n_j++){
							if (grid_lats[n_i][n_j] != -999 && (i != n_i || j != n_j) && grid_cluster[n_i][n_j] == cluster){
								double neighbor_distance = sqrt(pow(grid_lats[i][j] - grid_lats[n_i][n_j], 2) + pow(grid_lons[i][j] - grid_lons[n_i][n_j], 2));
								if (neighbor_distance > max_neighbor_distance){
									max_neighbor_distance = neighbor_distance;
								}
							}
						}
					}
				}
			}
		}
		if (max_neighbor_distance != -1 && distance_to_closest > max_neighbor_distance){
			eprintf("Point not found in file, %f degrees (distance to closest pixel) > %f (farthest neighbor)\n", distance_to_closest, max_neighbor_distance);
			for (i=0;i<line_count;i++){
				free(grid_lats[i]);
				free(grid_lons[i]);
				free(grid_cluster[i]);
			}
			free(grid_lats);
			free(grid_lons);
			free(grid_cluster);
			free(lats);
			free(lons);
			return VALEXTRACT_ERR_POINT_NOT_FOUND;
		} else {
			dprintf("Point found in file, %f degrees (distance to closest pixel) <= %f (farthest neighbor)\n", distance_to_closest, max_neighbor_distance);
		}

		unsigned box_shift = 3;
		nc_box *boxes = malloc(sizeof(nc_box) << box_shift);
		if (boxes == NULL){
			eprintf("Error allocating memory in '%s', near line %d\n", __FILE__, __LINE__);
			for (i=0;i<line_count;i++){
				free(grid_lats[i]);
				free(grid_lons[i]);
				free(grid_cluster[i]);
			}
			free(grid_lats);
			free(grid_lons);
			free(grid_cluster);
			free(lats);
			free(lons);
			return -1;
		}
		int box_max = (1 << box_shift) - 1;
		int box_i = 0;

		for (i=0;i<line_count;i++){
			for (j=0;j<pixel_count;j++){
				if (grid_cluster[i][j] == cluster){
					int this_start_pixel = j;
					for (;j<pixel_count && grid_cluster[i][j] == cluster;j++);
					nc_box *new_box_ptr = add_to_boxes(boxes, &box_shift, &box_i, &box_max, i+line_start, this_start_pixel+pixel_start, j+pixel_start-1);
					if (new_box_ptr == NULL){
						eprintf("Error allocating memory in '%s', near line %d\n", __FILE__, __LINE__);
						for (i=0;i<line_count;i++){
							free(grid_lats[i]);
							free(grid_lons[i]);
							free(grid_cluster[i]);
						}
						free(grid_lats);
						free(grid_lons);
						free(grid_cluster);
						free(boxes);
						free(lats);
						free(lons);
						return -1;
					}
					boxes = new_box_ptr;
				}
			}
		}

		region->box_count = box_i;
		region->boxes = boxes;

		for (i=0;i<line_count;i++){
			free(grid_lats[i]);
			free(grid_lons[i]);
			free(grid_cluster[i]);
		}
		free(grid_lats);
		free(grid_lons);
		free(grid_cluster);
	}

	free(lats);
	free(lons);
	return 0;
}

static void free_var(val_extract_arguments *arguments, nc_var *var){
	if (var){
		if (arguments->product_count == 0){
			free((void*)var->name);
		}
		if (var->units){
			free((void*)var->units);
		}
		if (var->data){
			free(var->data);
		}
		free(var->dim_ids);
	}
}

static int get_var_attributes(val_extract_arguments *arguments, nc_var *var){
	int nc_ret;
	var->data = NULL;

	if (arguments->product_count == 0){
		if ((var->name = malloc(sizeof(char) * 64)) == NULL){
			eprintf("Error allocating memory in '%s', near line %d\n", __FILE__, __LINE__);
			return -1;
		}
		if ((nc_ret = nc_inq_varname(var->gid, var->varid, (char*)var->name)) != NC_NOERR){
			eprintf("Error getting variable name, error %d\n", nc_ret);
			free((void*)var->name);
			return -1;
		}
	}

	if ((nc_ret = nc_inq_vartype(var->gid, var->varid, &var->data_type)) != NC_NOERR){
		eprintf("Error getting data type for %s, error %d\n", var->name, nc_ret);
		if (arguments->product_count == 0){
			free((void*)var->name);
		}
		return -1;
	}

	// this switch should go into the nc_read library that doesn't exist
	switch(var->data_type){
	case NC_BYTE:
	case NC_UBYTE:
	case NC_CHAR:
		var->data_sizeof = 1;
		break;
	case NC_SHORT:
	case NC_USHORT:
		var->data_sizeof = 2;
		break;
	case NC_INT: // case NC_LONG:
	case NC_UINT:
	case NC_FLOAT:
		var->data_sizeof = 4;
		break;
	case NC_DOUBLE:
	case NC_INT64:
	case NC_UINT64:
		var->data_sizeof = 8;
		break;
	case NC_STRING:
		var->data_sizeof = sizeof(char*);
		break;
	default:
		free(var->dim_ids);
		if (arguments->product_count == 0){
			free((void*)var->name);
		}
		return -1;
	}

	if ((nc_ret = nc_inq_varndims(var->gid, var->varid, &var->ndims)) != NC_NOERR){
		eprintf("Error getting dimension count for %s, error %d\n", var->name, nc_ret);
		if (arguments->product_count == 0){
			free((void*)var->name);
		}
		return -1;
	}
	if ((var->dim_ids = malloc(sizeof(int) * var->ndims)) == NULL){
		eprintf("Error allocating memory in '%s', near line %d\n", __FILE__, __LINE__);
		if (arguments->product_count == 0){
			free((void*)var->name);
		}
		return -1;
	}
	if ((nc_ret = nc_inq_vardimid(var->gid, var->varid, var->dim_ids)) != NC_NOERR){
		eprintf("Error getting dimensions for %s, error %d\n", var->name, nc_ret);
		free(var->dim_ids);
		if (arguments->product_count == 0){
			free((void*)var->name);
		}
		return -1;
	}

	var->is_geospatial = (var->ndims == 2 && var->dim_ids[0] == var->file->line_dimid && (var->dim_ids[1] == var->file->pixel_dimid || var->dim_ids[1] == var->file->pixel_control_dimid));
	var->uses_control_points = (var->ndims == 2 && var->dim_ids[1] == var->file->pixel_control_dimid);
	int fake_geospatial_i = 0;
	while (var->is_geospatial && ignore_geospatial[fake_geospatial_i]){
		if (!strcmp(var->name, ignore_geospatial[fake_geospatial_i++])){
			var->is_geospatial = 0;
		}
	}

	if ((nc_ret = nc_get_att_double(var->gid, var->varid, "scale_factor", &var->scale)) != NC_NOERR){
		if ((nc_ret = nc_get_att_double(var->gid, var->varid, "slope", &var->scale)) != NC_NOERR){
			var->has_scale = false;
		} else {
			var->has_scale = true;
		}
	} else {
		var->has_scale = true;
	}
	if ((nc_ret = nc_get_att_double(var->gid, var->varid, "add_offset", &var->offset)) != NC_NOERR){
		if ((nc_ret = nc_get_att_double(var->gid, var->varid, "intercept", &var->offset)) != NC_NOERR){
			var->has_offset = false;
		} else {
			var->has_offset = true;
		}
	} else {
		var->has_offset = true;
	}
	if ((nc_ret = nc_get_att_double(var->gid, var->varid, "_FillValue", &var->fill)) != NC_NOERR){
		if ((nc_ret = nc_get_att_double(var->gid, var->varid, "bad_value_scaled", &var->fill)) != NC_NOERR){
			var->has_fill = false;
		} else {
			var->has_fill = true;
		}
	} else {
		var->has_fill = true;
	}

	return 0;
}

static int process_variable(val_extract_arguments *arguments, nc_var *var){
	int nc_ret;
	nc_region *region = var->region;
	nc_file *file = var->file;
	int pixel_count = 0;
	//TODO: fix region stuff for geospatial vars that use control points
	if (var->is_geospatial || !arguments->geospatial_only){
		if (var->is_geospatial && arguments->geospatial_to_double){
			var->data_type = NC_DOUBLE;
			var->data_sizeof = sizeof(double);
		}
		void *data;
		double *sorted_data;
		double center_value = 0;
		if (var->is_geospatial){
//			pixel_count = nc_get_region_size(region);
			pixel_count = region->pixel_count;
			data = calloc(pixel_count, var->data_sizeof);
			sorted_data = calloc(pixel_count, sizeof(double));

			if (arguments->geospatial_to_double){
				if ((nc_ret = nc_get_varr_double(var->gid, var->varid, region, data)) != NC_NOERR){
					eprintf("Error getting data area for %s, error %d\n", var->name, nc_ret);
					free(sorted_data);
					return -1;
				}
			} else {
				if ((nc_ret = nc_get_varr(var->gid, var->varid, region, data)) != NC_NOERR){
					eprintf("Error getting data area for %s, error %d\n", var->name, nc_ret);
					free(sorted_data);
					return -1;
				}
			}

			if (var->data_type == NC_DOUBLE){
				memcpy(sorted_data, data, pixel_count << 3);
			} else {
				if ((nc_ret = nc_get_varr_double(var->gid, var->varid, region, sorted_data)) != NC_NOERR){
					eprintf("Error getting data area for %s, error %d\n", var->name, nc_ret);
					free(sorted_data);
					return -1;
				}
			}
//			dprintf("Center pixel coords: %zd, %zd\n", region->center.line, region->center.pixel);
			size_t center_point[2] = {region->center.line, region->center.pixel};
			if ((nc_ret = nc_get_var1_double(var->gid, var->varid, center_point, &center_value)) != NC_NOERR){
				eprintf("Error getting data center for %s, error %d\n", var->name, nc_ret);
				free(sorted_data);
				return -1;
			}
			if (!(var->has_fill && center_value == var->fill)){
				if (var->has_scale || var->has_offset){
					if (var->has_scale){
						center_value *= var->scale;
					}
					if (var->has_offset){
						center_value += var->offset;
					}
				}
			}
		} else {
			if (var->ndims == 1){
				if (var->dim_ids[0] == file->line_dimid){
					pixel_count = nc_get_region_dim_size(region, 0);
					data = calloc(pixel_count, var->data_sizeof);
					sorted_data = calloc(pixel_count, sizeof(double));
					if ((nc_ret = nc_get_varrd(var->gid, var->varid, region, 0, data)) != NC_NOERR){
						eprintf("Error getting data for %s, error %d\n", var->name, nc_ret);
						free(sorted_data);
						return -1;
					}
					if ((nc_ret = nc_get_varrd_double(var->gid, var->varid, region, 0, sorted_data)) != NC_NOERR){
						eprintf("Error getting data for %s, error %d\n", var->name, nc_ret);
						free(sorted_data);
						return -1;
					}
				} else if (var->dim_ids[0] == file->pixel_dimid || var->dim_ids[0] == file->pixel_control_dimid){
					pixel_count = nc_get_region_dim_size(region, 1);
					data = calloc(pixel_count, var->data_sizeof);
					sorted_data = calloc(pixel_count, sizeof(double));
					if ((nc_ret = nc_get_varrd(var->gid, var->varid, region, 1, data)) != NC_NOERR){
						eprintf("Error getting data for %s, error %d\n", var->name, nc_ret);
						free(sorted_data);
						return -1;
					}
					if ((nc_ret = nc_get_varrd_double(var->gid, var->varid, region, 1, sorted_data)) != NC_NOERR){
						eprintf("Error getting data for %s, error %d\n", var->name, nc_ret);
						free(sorted_data);
						return -1;
					}
				}
			}
			if (!pixel_count){
				uint8_t i;
				pixel_count = 1;
				for (i=0;i<var->ndims;i++){
					pixel_count *= file->dim_lengths[var->dim_ids[i]];
				}
				data = calloc(pixel_count, var->data_sizeof);
				sorted_data = calloc(pixel_count, sizeof(double));
				if ((nc_ret = nc_get_var(var->gid, var->varid, data)) != NC_NOERR){
					eprintf("Error getting data for %s, error %d\n", var->name, nc_ret);
					free(sorted_data);
					return -1;
				}
				if ((nc_ret = nc_get_var_double(var->gid, var->varid, sorted_data)) != NC_NOERR){
					eprintf("Error getting data for %s, error %d\n", var->name, nc_ret);
					free(sorted_data);
					return -1;
				}
			}
		}
		var->data = data;

		size_t attr_length;
		nc_type att_xtype;
		if ((nc_ret = nc_inq_att(var->gid, var->varid, "units", &att_xtype, &attr_length)) != NC_NOERR){
//			eprintf("Error inquiring units for %s: %d\n", var->name, nc_ret);
//			return VALEXTRACT_ERR_FLAG;
			var->units = NULL;
		} else {
			char *units_string = malloc((attr_length+1) * sizeof(char));
			var->units = units_string;
			if ((nc_ret = nc_get_att_text(var->gid, var->varid, "units", units_string)) != NC_NOERR){
				eprintf("Error getting units: %d\n", nc_ret);
				return VALEXTRACT_ERR_FLAG;
			}
			units_string[attr_length] = '\0';
		}

		val_extract_valid_range *valid_range = NULL;
		int i;
		for (i=0;i<arguments->valid_range_count && !valid_range;i++){
			if (starts_with(var->name, arguments->valid_ranges[i].prefix)){
				valid_range = &arguments->valid_ranges[i];
			}
		}

		int good_pixels = 0;
		for (i=0;i<pixel_count;i++){
			if (!(var->has_fill && sorted_data[i] == var->fill)){
				double v = sorted_data[i];
				if (var->is_geospatial && (region->pixel_flags[i] & arguments->ignore_flags_mask)){
					continue;
				}
				if (var->has_scale || var->has_offset){
					if (var->has_scale){
						v *= var->scale;
					}
					if (var->has_offset){
						v += var->offset;
					}

					switch (var->data_type){
					case NC_BYTE:
						((int8_t*)data)[i] = (int8_t)v;
						break;
					case NC_CHAR:
					case NC_UBYTE:
						((uint8_t*)data)[i] = (uint8_t)v;
						break;
					case NC_SHORT:
						((int16_t*)data)[i] = (int16_t)v;
						break;
					case NC_USHORT:
						((uint16_t*)data)[i] = (uint16_t)v;
						break;
					case NC_INT: // case NC_LONG:
						((int32_t*)data)[i] = (int32_t)v;
						break;
					case NC_UINT:
						((uint32_t*)data)[i] = (uint32_t)v;
						break;
					case NC_INT64:
						((int64_t*)data)[i] = (int64_t)v;
						break;
					case NC_UINT64:
						((uint64_t*)data)[i] = (uint64_t)v;
						break;
					case NC_FLOAT:
						((float*)data)[i] = (float)v;
						break;
					case NC_DOUBLE:
						((double*)data)[i] = (double)v;
						break;
					default:
						free(sorted_data);
						return -1;
					}
				}

				if (!valid_range || (v >= valid_range->min && v <= valid_range->max)){
					sorted_data[good_pixels++] = v;
				}
			}
		}

		if (good_pixels){
			if (!var->is_geospatial){
				center_value = sorted_data[good_pixels >> 1];
			}
			qsort(sorted_data, good_pixels, sizeof(double), compare_double);

			var->stats = generate_stats(sorted_data, 0, good_pixels);
			var->stats.center_value = center_value;

			int filter_start_i = -1, filter_end_i = -1;
			const double filter_start = var->stats.median - (1.5 * var->stats.stddev);
			const double filter_end = var->stats.median + (1.5 * var->stats.stddev);
			get_array_ranges(sorted_data, good_pixels, filter_start, filter_end, &filter_start_i, &filter_end_i);
			int filtered_length = filter_end_i - filter_start_i + 1;
			var->filtered_stats = generate_stats(sorted_data, filter_start_i, filtered_length);

			int iqr_start_i = -1, iqr_end_i = -1;
			const double iqr_start = sorted_data[(int)(((good_pixels + 1) / 4))];
			const double iqr_end = sorted_data[(int)(((good_pixels + 1) * 3 / 4) - 1)];
			get_array_ranges(sorted_data, good_pixels, iqr_start, iqr_end, &iqr_start_i, &iqr_end_i);
			int iqr_length = iqr_end_i - iqr_start_i + 1;
			if (iqr_length >= 8){
				var->iqr_stats = generate_stats(sorted_data, iqr_start_i, iqr_length);
				var->iqr_stats.count = iqr_length;
			} else {
				nc_var_stats iqr_stats = {
					.initialized = true, .count = 0, .center_value = 0,
					.min = 0, .max = 0, .mean = 0, .median = 0, .stddev = 0, .rms = 0
				};
				var->iqr_stats = iqr_stats;
			}
		} else {
			nc_var_stats stats = {
				.initialized = true, .count = 0, .center_value = 0,
				.min = 0, .max = 0, .mean = 0, .median = 0, .stddev = 0, .rms = 0
			};
			var->stats = var->filtered_stats = var->iqr_stats = stats;
			var->stats.center_value = center_value;
		}

		free(sorted_data);
	}

	return arguments->val_extract_parser(VALEXTRACT_KEY_VAR, var, arguments->user_input);
}

void get_array_ranges(const double *sorted_data, const int length, const double start, const double end, int *start_i, int *end_i){
	*start_i = *end_i = (length >> 1);
	int i;
	for (i=0;i<length;i++){
		if (sorted_data[i] >= start){
			*start_i = i;
			break;
		}
	}
	for (i=length-1;i>=0;i--){
		if (sorted_data[i] <= end){
			*end_i = i;
			break;
		}
	}
}


nc_var_stats generate_stats(const double *sorted_data, const int offset, const int length){
	int i;
	double total = 0;
	for (i=offset;i<offset+length;i++){
		total += sorted_data[i];
	}
	double min, max, mean, median;
	if (length & 1){
		median = sorted_data[offset + (length >> 1)];
	} else {
		median = (sorted_data[offset + (length >> 1) - 1] + sorted_data[offset + (length >> 1)]) / 2;
	}
	min = sorted_data[offset];
	max = sorted_data[offset + length - 1];
	mean = total / length;
	double stddev = 0, rms = 0;
	for (i=offset;i<offset+length;i++){
		stddev += pow(mean - sorted_data[i], 2);
		rms += pow(sorted_data[i], 2);
	}
	if (length > 1){
//		stddev = sqrt(stddev / length);     // population standard deviation
		stddev = sqrt(stddev / (length-1)); // sample standard deviation
		rms = sqrt(rms / length);
	} else {
		stddev = rms = 0;
	}

	nc_var_stats stats = {
		.initialized = true, .count = length, .center_value = median,
		.min = min, .max = max, .mean = mean, .median = median, .stddev = stddev, .rms = rms
	};
	return stats;
}

static int process_file_attrs(val_extract_arguments *arguments, nc_region *region){
	nc_file *file = region->file;
	int nc_ret, pixel_count = nc_get_region_size(region);
	int flags_varid, flags_gid, gid_i;
	flags_varid = flags_gid = -1;
	for (gid_i=0;gid_i<file->ngrps;gid_i++){
		if ((nc_ret = nc_inq_varid(file->gids[gid_i], "l2_flags", &flags_varid)) == NC_NOERR){
			flags_gid = file->gids[gid_i];
			break;
		}
	}
	if (gid_i == -1 || flags_varid == -1){
		eprintf("l2_flags not found, exiting\n");
		return VALEXTRACT_ERR_FLAG;
	}

	int i, flag_i;
	nc_type att_xtype;
	size_t u_flag_count;
    int flag_count;

    uint32_t *flag_masks;
	char **flag_meanings;
	int *flag_counts;

	if ((nc_ret = nc_inq_att(flags_gid, flags_varid, "flag_masks", &att_xtype, &u_flag_count)) != NC_NOERR){
		flag_masks = calloc(32, sizeof(uint32_t));
		flag_meanings = malloc(32 * sizeof(char*));
		flag_counts = calloc(32, sizeof(int));
		file->flag_masks = flag_masks;
		file->flag_meanings = flag_meanings;
		region->flag_counts = flag_counts;

		int flag_string_length = 0;
		char att_name[10];
		int i=0;
		for (i=0,u_flag_count=0;i<32;i++){
			sprintf(att_name, "f%02zd_name", u_flag_count+1);
			size_t attr_length;
			if ((nc_ret = nc_inq_att(flags_gid, flags_varid, att_name, &att_xtype, &attr_length)) != NC_NOERR){
				break;
			} else {
				u_flag_count++;
				flag_string_length += attr_length;
				flag_meanings[i] = malloc((attr_length+1) * sizeof(char));
				if ((nc_ret = nc_get_att_text(flags_gid, flags_varid, att_name, flag_meanings[i])) != NC_NOERR){
					eprintf("Error getting %s: %d\n", att_name, nc_ret);
					return VALEXTRACT_ERR_FLAG;
				}
			}
		}
		if (i){
			if (u_flag_count == 33){
				--u_flag_count;
			}
		    flag_count = u_flag_count;
			file->flag_string = malloc((flag_string_length+1) * sizeof(char));
			int flag_string_pos = 0;
			for (i=0;i<flag_count;i++){
				flag_masks[i] = 1 << i;
				if (i){
					sprintf(&file->flag_string[flag_string_pos++], " ");
				}
				sprintf(&file->flag_string[flag_string_pos], "%s", flag_meanings[i]);
				flag_string_pos += strlen(flag_meanings[i]);
			}
		} else {
			eprintf("Error inquiring flag_masks: %d\n", nc_ret);
			return VALEXTRACT_ERR_FLAG;
		}
	} else {
	    flag_count = u_flag_count;
		flag_masks = malloc(u_flag_count * sizeof(uint32_t));
		flag_meanings = malloc(u_flag_count * sizeof(char*));
		flag_counts = calloc(u_flag_count, sizeof(int));
		file->flag_masks = flag_masks;
		file->flag_meanings = flag_meanings;
		region->flag_counts = flag_counts;
		if ((nc_ret = nc_get_att_int(flags_gid, flags_varid, "flag_masks", (int*)flag_masks)) != NC_NOERR){
			eprintf("Error getting flag_masks: %d\n", nc_ret);
			return VALEXTRACT_ERR_FLAG;
		}

		size_t attr_length;
		if ((nc_ret = nc_inq_att(flags_gid, flags_varid, "flag_meanings", &att_xtype, &attr_length)) != NC_NOERR){
			eprintf("Error inquiring flag_meanings: %d\n", nc_ret);
			return VALEXTRACT_ERR_FLAG;
		}
		char *flag_string = malloc((attr_length+1) * sizeof(char));
		file->flag_string = flag_string;
		if ((nc_ret = nc_get_att_text(flags_gid, flags_varid, "flag_meanings", flag_string)) != NC_NOERR){
			eprintf("Error getting flag_meanings: %d\n", nc_ret);
			return VALEXTRACT_ERR_FLAG;
		}
		flag_string[attr_length] = '\0';
		dprintf("%s\n", flag_string);

		i = 0;
	    char *token = strtok(flag_string, FLAG_MEANINGS_SEP);
	    while (token) {
	    	flag_meanings[i++] = token;
	        token = strtok(NULL, FLAG_MEANINGS_SEP);
	    }
	}


	size_t attr_length;
	if ((nc_ret = nc_inq_att(file->ncid, NC_GLOBAL, "platform", &att_xtype, &attr_length)) != NC_NOERR){
		eprintf("Error inquiring platform: %d\n", nc_ret);
		return VALEXTRACT_ERR_FLAG;
	}
	char *platform_string = malloc((attr_length+1) * sizeof(char));
	file->platform = platform_string;
	if ((nc_ret = nc_get_att_text(file->ncid, NC_GLOBAL, "platform", platform_string)) != NC_NOERR){
		eprintf("Error getting platform: %d\n", nc_ret);
		return VALEXTRACT_ERR_FLAG;
	}
	platform_string[attr_length] = '\0';
	dprintf("Platform: %s\n", platform_string);
	
	if ((nc_ret = nc_inq_att(file->ncid, NC_GLOBAL, "instrument", &att_xtype, &attr_length)) != NC_NOERR){
		eprintf("Error inquiring instrument: %d\n", nc_ret);
		return VALEXTRACT_ERR_FLAG;
	}
	char *instrument_string = malloc((attr_length+1) * sizeof(char));
	file->instrument = instrument_string;
	if ((nc_ret = nc_get_att_text(file->ncid, NC_GLOBAL, "instrument", instrument_string)) != NC_NOERR){
		eprintf("Error getting instrument: %d\n", nc_ret);
		return VALEXTRACT_ERR_FLAG;
	}
	instrument_string[attr_length] = '\0';
	dprintf("Instrument: %s\n", instrument_string);



    file->flag_count = flag_count;

    uint32_t *pixel_flags = malloc(sizeof(uint32_t) * pixel_count);
	region->pixel_flags = pixel_flags;
	if ((nc_ret = nc_get_varr_uint(flags_gid, flags_varid, region, pixel_flags)) != NC_NOERR){
		eprintf("Error getting flags data: %d\n", nc_ret);
		return VALEXTRACT_ERR_FLAG;
	}

	dprintf("Ignore count: %d\n", arguments->ignore_flag_count);
	for (i=0;i<arguments->ignore_flag_count;i++){
		dprintf("Ignore: %s\n", arguments->ignore_flags[i]);
		for (flag_i=0;flag_i<flag_count;flag_i++){
			if (strcmp(flag_meanings[flag_i], arguments->ignore_flags[i]) == 0){
				arguments->ignore_flags_mask |= flag_masks[flag_i];
				break;
			}
		}
		if (flag_i == flag_count){
			eprintf("Ignore flag %s not found, exiting\n", arguments->l2qc_flags[i]);
			return VALEXTRACT_ERR_FLAG;
		}
	}

	for (i=0;i<pixel_count;i++){
		for (flag_i=0;flag_i<flag_count;flag_i++){
			if (pixel_flags[i] & flag_masks[flag_i]){
				flag_counts[flag_i]++;
			}
		}
	}

	int unflagged_pixel_count = pixel_count;
	uint32_t flag_mask = arguments->ignore_flags_mask;
	dprintf("Flag mask: %u\n", arguments->ignore_flags_mask);
	for (i=0;i<pixel_count;i++){
		if (pixel_flags[i] & flag_mask){
			--unflagged_pixel_count;
		}
	}

	dprintf("Flags to count: %d\n", arguments->count_flag_count);
	for (i=0;i<arguments->count_flag_count;i++){
		dprintf("Count: %s\n", arguments->count_flags[i]);
		for (flag_i=0;flag_i<flag_count;flag_i++){
			if (strcmp(flag_meanings[flag_i], arguments->count_flags[i]) == 0){
				arguments->count_flags_mask |= flag_masks[flag_i];
				break;
			}
		}
		if (flag_i == flag_count){
			eprintf("Ignore flag %s not found, exiting\n", arguments->l2qc_flags[i]);
			return VALEXTRACT_ERR_FLAG;
		}
	}
	int flagged_pixel_count = 0;
	uint32_t count_flag_mask = arguments->count_flags_mask;
	dprintf("Flag count mask: %u\n", arguments->count_flags_mask);
	for (i=0;i<pixel_count;i++){
		if (pixel_flags[i] & count_flag_mask){
			++flagged_pixel_count;
		}
	}

	dprintf("L2QC count: %d\n", arguments->l2qc_threshold_count);
	int flag_ret = 0;
	for (i=0;i<arguments->l2qc_threshold_count;i++){
		for (flag_i=0;flag_i<flag_count;flag_i++){
			if (strcmp(flag_meanings[flag_i], arguments->l2qc_flags[i]) == 0){
				if (arguments->l2qc_thresholds[i] >= 0){
					double flag_percent = (double)flag_counts[flag_i]/pixel_count;
					if (flag_percent > arguments->l2qc_thresholds[i]){
						flag_ret |= flag_masks[flag_i];
						dprintf("%s failed (%.3f < %f)\n", flag_meanings[flag_i], arguments->l2qc_thresholds[i], flag_percent);
					} else {
						dprintf("%s clear (%.3f > %f)\n", flag_meanings[flag_i], arguments->l2qc_thresholds[i], flag_percent);
					}
				}
				break;
			}
		}
		if (flag_i == flag_count){
			eprintf("L2QC flag %s not found, exiting\n", arguments->l2qc_flags[i]);
			return VALEXTRACT_ERR_FLAG;
		}
	}

	errno = 0;
	int year = process_scan_line_att(arguments, region, "year");
	int day = process_scan_line_att(arguments, region, "day");
	int msec = process_scan_line_att(arguments, region, "msec");
	double lat = process_location_var(arguments, region, "latitude");
	double lon = process_location_var(arguments, region, "longitude");

	if (errno){
		eprintf("Error processing scan line attributes\n");
		return VALEXTRACT_ERR_NCFILE_INVALID;
	}

	struct tm t = { .tm_year = year - 1900, .tm_mday = day, .tm_sec = msec / 1000 };
	region->utime = mktime(&t);
	if (t.tm_isdst){
		t.tm_sec -= 3600;
		region->utime = mktime(&t);
	}
	sprintf(region->ascii_time, "%04d-%02d-%02d %02d:%02d:%02d.%03d", t.tm_year + 1900, t.tm_mon+1, t.tm_mday, t.tm_hour, t.tm_min, t.tm_sec, msec % 1000);
	region->lat = lat;
	region->lon = lon;
	region->pixel_count = pixel_count;
	region->unflagged_pixel_count = unflagged_pixel_count;
	region->flagged_pixel_count = flagged_pixel_count;

	int callback_ret = arguments->val_extract_parser(VALEXTRACT_KEY_FILE, region, arguments->user_input);

	if (callback_ret){
		if (callback_ret != VALEXTRACT_CMD_STOP){
			eprintf("val_extract_parser error\n");
			return VALEXTRACT_ERR_UNKNOWN;
		} else {
			return callback_ret;
		}
	}
	if (flag_ret){
		eprintf("File failed L2QC check\n");
		return VALEXTRACT_ERR_L2QC;
	}
	return 0;
}

static double process_scan_line_att(val_extract_arguments *arguments, nc_region *region, const char *name){
	int nc_ret, varid, gid, gid_i;
	varid = gid = -1;
	for (gid_i=0;gid_i<region->file->ngrps;gid_i++){
		if ((nc_ret = nc_inq_varid(region->file->gids[gid_i], name, &varid)) == NC_NOERR){
			gid = region->file->gids[gid_i];
			break;
		}
	}
	if (gid_i == -1 || varid == -1){
		eprintf("scan line attribute %s not found\n", name);
		errno |= 1;
		return 0;
	}

	double line_value = 0;
	if ((nc_ret = nc_get_var1_double(gid, varid, &region->center.line, &line_value)) != NC_NOERR){
		eprintf("Error getting %s data, %d\n", name, nc_ret);
		errno |= 1;
		return 0;
	}

	return line_value;
}

static double process_location_var(val_extract_arguments *arguments, nc_region *region, const char *name){
	int nc_ret, varid, gid, gid_i;
	varid = gid = -1;
	for (gid_i=0;gid_i<region->file->ngrps;gid_i++){
		if ((nc_ret = nc_inq_varid(region->file->gids[gid_i], name, &varid)) == NC_NOERR){
			gid = region->file->gids[gid_i];
			break;
		}
	}
	if (gid_i == -1 || varid == -1){
		eprintf("location variable %s not found\n", name);
		errno |= 1;
		return 0;
	}

	double center_value = 0;
	if ((nc_ret = nc_get_var1_double(gid, varid, &region->center.line, &center_value)) != NC_NOERR){
		eprintf("Error getting %s data, %d\n", name, nc_ret);
		errno |= 1;
		return 0;
	}

	return center_value;
}

int val_extract_clean(val_extract_arguments *arguments){
	if (arguments->products){
		free(arguments->products);
	}
	if (arguments->ignore_flags){
		int i;
		for (i=0;i<arguments->ignore_flag_count;i++){
			free(arguments->ignore_flags[i]);
		}
		free(arguments->ignore_flags);
	}
	if (arguments->count_flags){
		int i;
		for (i=0;i<arguments->count_flag_count;i++){
			free(arguments->count_flags[i]);
		}
		free(arguments->count_flags);
	}
	if (arguments->l2qc_flags){
		int i;
		for (i=0;i<arguments->l2qc_threshold_count;i++){
			free(arguments->l2qc_flags[i]);
		}
		free(arguments->l2qc_flags);
		free(arguments->l2qc_thresholds);
	}
	if (arguments->valid_ranges){
		int i;
		for (i=0;i<arguments->valid_range_count;i++){
			free(arguments->valid_ranges[i].prefix);
		}
		free(arguments->valid_ranges);
	}
	return 0;
}


static void print_time(){
	time_t epoch_time = time(NULL);
	struct tm *current_time = localtime(&epoch_time);
	char time_buffer[64];
	if (strftime(time_buffer, sizeof(time_buffer), "%A %c", current_time)){
		dprintf("%s\n", time_buffer);
	}
}

// http://stackoverflow.com/a/4771038/1236508 (with reordered and renamed args)
static bool starts_with(const char *str, const char *prefix){
    size_t lenstr = strlen(str), lenpre = strlen(prefix);
    return lenstr < lenpre ? false : strncmp(prefix, str, lenpre) == 0;
}

void unflatten_index(int index, int ndims, const int *dims, int *result){
	int j;
	int multiple[ndims];
	multiple[ndims-1] = 1;
	for (j = ndims-2; j >= 0; j--){
		multiple[j] = multiple[j+1]*dims[j+1];
	}
	for (j=0;j<ndims;j++){
		result[j] = (index / multiple[j]) % dims[j];
	}
}

static void calc_dim_multipliers(int ndims, const int *dims, int *dim_multipliers){
	int j;
	dim_multipliers[ndims-1] = 1;
	for (j = ndims-2; j >= 0; j--){
		dim_multipliers[j] = dim_multipliers[j+1]*dims[j+1];
	}
}

static int flatten_index(int ndims, const int *dim_multipliers, int *index){
	int ret = 0, j;
	for (j = ndims-1; j >= 0; j--){
		ret += (index[j] * dim_multipliers[j]);
	}
	return ret;
}

const char* val_extract_version(){
	return VALEXTRACT_IMPLEMENTATION_STR;
}
const char* val_extract_api_version(){
	return VALEXTRACT_API_VERSION_STR;
}

// The following should go into a math library

#define define_compare(type, suffix) \
int compare ## suffix(const void *a, const void *b){ \
	if (*(const type*)a < *(const type*)b){ \
		return -1; \
	} else if (*(const type*)a > *(const type*)b){ \
		return 1; \
	} else { \
		return 0; \
	} \
}

define_compare(int8_t, _int8)
define_compare(uint8_t, _uint8)
define_compare(int16_t, _int16)
define_compare(uint16_t, _uint16)
define_compare(int32_t, _int32)
define_compare(uint32_t, _uint32)
define_compare(int64_t, _int64)
define_compare(uint64_t, _uint64)
define_compare(float, _float)
define_compare(double, _double)

static double strtod_strict(const char *str){
	char *endptr;
	errno = 0;
	double ret = strtod(str, &endptr);
	if (!errno && *endptr != '\0'){
		errno = 1;
	}
	return ret;
}

// The following should go into an nc_read library

size_t nc_get_region_size(nc_region *region){
	int i;
	size_t ret = 0;
	for (i = 0; i < region->box_count; i++){
		ret += region->boxes[i].count[0] * region->boxes[i].count[1];
	}
	return ret;
}
size_t nc_get_region_dim_size(nc_region *region, int dim_index){
	int i;
	size_t ret = 0;
	for (i = 0; i < region->box_count; i++){
		ret += region->boxes[i].count[dim_index];
	}
	return ret;
}

#define define_nc_get_varr(type, suffix) \
int nc_get_varr ## suffix(int ncid, int varid, nc_region *region, type *p){ \
	int i, nc_ret = 0; \
	for (i=0;i<region->box_count;i++){ \
		if ((nc_ret = nc_get_vara ## suffix(ncid, varid, region->boxes[i].start, region->boxes[i].count, p)) != NC_NOERR){ \
			break; \
		} \
		p += region->boxes[i].count[0] * region->boxes[i].count[1]; \
	} \
	return nc_ret; \
} \
int nc_get_varrd ## suffix(int ncid, int varid, nc_region *region, int dim_index, type *p){ \
	int i, nc_ret = 0; \
	for (i=0;i<region->box_count;i++){ \
		if ((nc_ret = nc_get_vara ## suffix(ncid, varid, region->boxes[i].start + dim_index, region->boxes[i].count + dim_index, p)) != NC_NOERR){ \
			break; \
		} \
		p += region->boxes[i].count[dim_index]; \
	} \
	return nc_ret; \
}

define_nc_get_varr(char, _text)
define_nc_get_varr(unsigned char, _uchar)
define_nc_get_varr(signed char, _schar)
define_nc_get_varr(short, _short)
define_nc_get_varr(int, _int)
define_nc_get_varr(long, _long)
define_nc_get_varr(float, _float)
define_nc_get_varr(double, _double)
define_nc_get_varr(unsigned short, _ushort)
define_nc_get_varr(unsigned int, _uint)
define_nc_get_varr(long long, _longlong)
define_nc_get_varr(unsigned long long, _ulonglong)
define_nc_get_varr(char *, _string)
define_nc_get_varr(void, )
