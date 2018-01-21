// implemented version numbers
#define SSTVALEXTRACT_IMPLEMENTATION 1000002
#define SSTVALEXTRACT_IMPLEMENTATION_STR "1.0.2"
#define SSTVALEXTRACT_API_VERSION 1000002
#define SSTVALEXTRACT_API_VERSION_STR "1.0.2"

#include "argpar.h"
#include "sstval_extract.h"
#include "val_extract.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static argpar_option options[] = {
	{ "header", 'd', "1", OPTION_INT, "print header before data" },
	{ "buoyid", 'b', "ID", 0, "buoy ID" },
	{ "buoytime", 't', "YYDDDHHMMSS", 0, "buoy time" },
	{ 0, 0, 0, OPTION_DOC, "If header=1, ifile is not required.  The header will be printed and the program will exit." },
	{ 0 } // tell argpar to stop checking options
};

#ifdef DEBUG
#define dprintf(...) do { printf(__VA_ARGS__); } while(0)
#else
#define dprintf(...) do {} while (0)
#endif

#define eprintf(...) do { fprintf(stderr, __VA_ARGS__); } while(0)

static const char *blank_string_for_ifile = "";

static int parse_options(int key, char *argv, struct argpar_state *state) {
	int ret = 0;

	sstval_extract_arguments *arguments = state->input;
	switch (key){
	case 'd':
		arguments->header = state->argv_as_int;
		break;
	case 'b':
		arguments->buoy_id = argv;
		break;
	case 't':
		arguments->buoy_time = argv;
		break;
	case ARGPAR_KEY_ARG:
		if (arguments->req_file){
			eprintf("Too many arguments.\n");
			return ARGPAR_ERR_USAGE;
		} else if (arguments->req_id){
			arguments->req_file = argv;
		} else {
			arguments->req_id = argv;
		}
		break;
	case ARGPAR_KEY_END:
		if (!arguments->val_extract_arguments.ifile){
			arguments->val_extract_arguments.ifile = blank_string_for_ifile;
		}
		break;
	case ARGPAR_KEY_SUCCESS:
		if (arguments->val_extract_arguments.ifile){
			if (!arguments->buoy_id || !arguments->buoy_time){
				ret = SSTVALEXTRACT_ERR_INPUT;
			}
		} else {
			if (!arguments->header){
				ret = SSTVALEXTRACT_ERR_INPUT;
			}
		}
		break;
	case ARGPAR_KEY_INIT:
		state->child_inputs[0] = &arguments->val_extract_arguments;
		break;
	}
	return ret;
}

static const argpar_child argpar_children[] = {{&val_extract_argpar}, {0}};
argpar sstval_extract_argpar = { options, parse_options, NULL, NULL, argpar_children };

typedef struct val_extract_parser_input {
	char *ofile;
	unsigned printed_products;
	FILE *base_ofile_h;
} val_extract_main_input;

static const double ref_bt3959minus4050[] = {
		4.44, 4.22, 4.02, 3.97, 3.90, 3.73, 3.55, 3.44, 3.35, 3.29, 3.18, 3.09, 3.03, 2.99,
		2.93, 2.88, 2.81, 2.77, 2.74, 2.70, 2.67, 2.64, 2.58, 2.56, 2.54, 2.51, 2.46, 2.46,
		2.43, 2.41, 2.38, 2.35, 2.35, 2.33, 2.32, 2.29, 2.27, 2.25, 2.24, 2.26, 2.24, 2.21,
		2.20, 2.21, 2.19, 2.19, 2.18, 2.17, 2.17, 2.16, 2.17, 2.15, 2.16, 2.16, 2.15, 2.16,
		2.16, 2.16, 2.16, 2.16, 2.17, 2.15, 2.14, 2.14, 2.17, 2.21, 2.20, 2.19, 2.19, 2.19,
		2.18, 2.21, 2.26, 2.30, 2.31, 2.30, 2.34, 2.32, 2.35, 2.40, 2.42, 2.45, 2.42, 2.50,
		2.55, 2.58, 2.65, 2.71, 2.77, 2.82, 2.89, 2.96, 3.02, 3.10, 3.16, 3.22, 3.31, 3.38,
		3.45, 3.55, 3.66, 3.80, 3.94, 4.08, 4.30};

static int calculate_bt3959minus4050(double bt3959, double bt4050, int pixel_number){
	// linear interpolation of ref_bt3959minus4050 using pixel number
	double aijpx = ((pixel_number - 1.0)/13.0)+1.0;
	int ijpx = (int)aijpx;
	aijpx -= ijpx;
	double bijpx = 1.0 - aijpx;
	double tref = ref_bt3959minus4050[ijpx];
	if (ijpx < 105){ // 105 == size of ref_bt3959minus4050 - 1?
		tref = (ref_bt3959minus4050[ijpx] * bijpx) + (ref_bt3959minus4050[ijpx+1] * aijpx);
	}
	return bt3959 - bt4050 - tref;
}

static int read_val_extract(int key, void *nc_input, void *user_input){
	switch(key){
	case VALEXTRACT_KEY_INIT:
		if (user_input){ // useless check for dumb compiler problem
			sstval_extract_output *o = (sstval_extract_output*)user_input;
			o->lat = o->lon = o->solz = o->senz = o->angle_of_incidence = o->sena = o->sola = o->relaz = o->quality_sst = o->quality_sst4 = FLT_MAX;
			o->year = o->day_of_year = o->hour = o->mirror = o->detector_number = o->pixel_number = o->glint_flag = o->l2_flags = UINT_MAX;
		}
		break;
	case VALEXTRACT_KEY_SUCCESS:
		if (user_input){ // useless check for dumb compiler problem
			nc_region *region = nc_input;
			nc_file *file = region->file;
			sstval_extract_output *output = (sstval_extract_output*)user_input;

			output->angle_of_incidence = 10.5 + output->pixel_number / 1354.0 * (65.5 - 10.5);
			output->relaz = output->sena - 180 - output->sola;
			unsigned glint_mask = 0;
			int flag_i;
			for (flag_i = 0; flag_i < file->flag_count && !glint_mask; flag_i++){
				if (strcmp(file->flag_meanings[flag_i], "HIGLINT") == 0){
					glint_mask = file->flag_masks[flag_i];
				}
			}
			output->glint_flag = output->l2_flags & glint_mask;

			char *basename = (char*)file->file_path + strlen(file->file_path);
			while (basename >= file->file_path && *basename != '/'){
				basename--;
			}
			output->sensor = *(basename + 1);
			output->granule = malloc(12 * sizeof(char));
			output->granule[11] = '\0';
			strncpy(output->granule, basename + 4, 11);
			int grp_i;
			size_t cal_length = 0;
			for (grp_i = 0; grp_i < file->ngrps && !cal_length; grp_i++){
				if (nc_inq_attlen(file->gids[grp_i], NC_GLOBAL, "calibration_data", &cal_length) == NC_NOERR){
					char tmp[cal_length+1];
					nc_get_att_text(file->gids[grp_i], NC_GLOBAL, "calibration_data", tmp);
					tmp[cal_length] = '\0';
					char *version_start = strstr(tmp, "MYD02_Reflective_LUTs.") + strlen("MYD02_Reflective_LUTs.");
					char *version_end = strstr(version_start, ".hdf");
					int version_length = (version_end - version_start) / sizeof(char);
					output->reflective_lut_version = malloc((version_length + 1) * sizeof(char));
					output->reflective_lut_version[version_length] = '\0';
					strncpy(output->reflective_lut_version, version_start, version_length);
				}
			}

			if (output->pixel_number != UINT_MAX && output->bt_3959_var && output->bt_4050_var){
				nc_var *bt_3959 = output->bt_3959_var;
				nc_var *bt_4050 = output->bt_4050_var;
				nc_var *pixel_numbers = output->bt_4050_var;
				int i;
				double total = 0;
				double *values = malloc(region->pixel_count * sizeof(double));
				double *at = values;
				for (i=0;i<region->pixel_count;i++){
					if (!(bt_3959->has_fill && ((double*)bt_3959->data)[i] == bt_3959->fill) && !(bt_4050->has_fill && ((double*)bt_4050->data)[i] == bt_4050->fill)){
						double v22 = ((double*)bt_3959->data)[i];
						if (bt_3959->has_scale || bt_3959->has_offset){
							if (bt_3959->has_scale){
								v22 *= bt_3959->scale;
							}
							if (bt_3959->has_offset){
								v22 += bt_3959->offset;
							}
						}
						double v23 = ((double*)bt_4050->data)[i];
						if (bt_4050->has_scale || bt_4050->has_offset){
							if (bt_4050->has_scale){
								v23 *= bt_4050->scale;
							}
							if (bt_4050->has_offset){
								v23 += bt_4050->offset;
							}
						}

						*at++ = calculate_bt3959minus4050(v22, v23, ((double*)pixel_numbers->data)[i]);
					}
				}
				int unflagged_pixel_count = at - values;

				qsort(values, unflagged_pixel_count, sizeof(double), compare_double);

				double min, max, mean, median;
				median = values[unflagged_pixel_count >> 1];
				min = values[0];
				max = values[unflagged_pixel_count - 1];
				mean = total / unflagged_pixel_count;
				double stddev = 0, rms = 0;
				for (i=0;i<unflagged_pixel_count;i++){
					stddev += pow(mean - values[i], 2);
					rms += pow(values[i], 2);
				}
//				stddev = sqrt(stddev / unflagged_pixel_count); // population standard deviation
				stddev = sqrt(stddev / (unflagged_pixel_count-1)); // sample standard deviation
				rms = sqrt(rms / unflagged_pixel_count);

				free(values);

				output->bt3959minus4050.initialized = 1;
				output->bt3959minus4050.count = unflagged_pixel_count;
				output->bt3959minus4050.min = min;
				output->bt3959minus4050.max = max;
				output->bt3959minus4050.mean = mean;
				output->bt3959minus4050.median = median;
				output->bt3959minus4050.stddev = stddev;
				output->bt3959minus4050.rms = rms;
				output->bt3959minus4050.center_value = calculate_bt3959minus4050(bt_3959->stats.center_value, bt_4050->stats.center_value, pixel_numbers->stats.center_value);
			}
		}
		break;
	case VALEXTRACT_KEY_ERROR:
		break;
	case VALEXTRACT_KEY_FINI:
		break;
	case VALEXTRACT_KEY_FILE:
		if (nc_input){ // useless check for dumb compiler problem
			sstval_extract_output *output = (sstval_extract_output*)user_input;
			nc_region *region = nc_input;
//			nc_file *file = region->file;

			output->region = *region;
			output->lat = region->lat;
			output->lon = region->lon;
		}
		break;
	case VALEXTRACT_KEY_VAR:
		if (nc_input){ // useless check for dumb compiler problem
			sstval_extract_output *output = (sstval_extract_output*)user_input;
			nc_var *var = nc_input;
//			nc_file *file = var->file;
			if (!var->is_geospatial){
				if (strcmp(var->name, "year") == 0){
					output->year = (unsigned)var->stats.mean;
				} else if (strcmp(var->name, "day") == 0){
					output->day_of_year = (unsigned)var->stats.mean;
				} else if (strcmp(var->name, "msec") == 0){
					output->hour = var->stats.mean / 3600000.0;
				} else if (strcmp(var->name, "mside") == 0){
					output->mirror = (unsigned)var->stats.center_value;
				} else if (strcmp(var->name, "detnum") == 0){
					output->detector_number = (unsigned)var->stats.center_value;
				}
			} else {
				if (strcmp(var->name, "sst") == 0){
					output->sst = var->stats;
				} else if (strcmp(var->name, "sst4") == 0){
					output->sst4 = var->stats;
				} else if (strcmp(var->name, "BT_3750") == 0){
					output->bt_3750 = var->stats;
				} else if (strcmp(var->name, "BT_3959") == 0){
					output->bt_3959 = var->stats;
					output->bt_3959_var = var;
				} else if (strcmp(var->name, "BT_4050") == 0){
					output->bt_4050 = var->stats;
					output->bt_4050_var = var;
				} else if (strcmp(var->name, "BT_6715") == 0){
					output->bt_6715 = var->stats;
				} else if (strcmp(var->name, "BT_7325") == 0){
					output->bt_7325 = var->stats;
				} else if (strcmp(var->name, "BT_8550") == 0){
					output->bt_8550 = var->stats;
				} else if (strcmp(var->name, "BT_11000") == 0){
					output->bt_11000 = var->stats;
				} else if (strcmp(var->name, "BT_12000") == 0){
					output->bt_12000 = var->stats;
				} else if (strcmp(var->name, "rhot_678") == 0){
					output->rhot_678 = var->stats;
				} else if (strcmp(var->name, "rho_cirrus") == 0){
					output->rho_cirrus = var->stats;
				} else if (strcmp(var->name, "sstref") == 0){
					output->sstref = var->stats;
				} else if (strcmp(var->name, "qual_sst") == 0){
					output->quality_sst = var->stats.center_value;
				} else if (strcmp(var->name, "qual_sst4") == 0){
					output->quality_sst4 = var->stats.center_value;
				} else if (strcmp(var->name, "senz") == 0){
					output->senz = var->stats.center_value;
				} else if (strcmp(var->name, "solz") == 0){
					output->solz = var->stats.center_value;
				} else if (strcmp(var->name, "sena") == 0){
					output->sena = var->stats.center_value;
				} else if (strcmp(var->name, "sola") == 0){
					output->sola = var->stats.center_value;
				} else if (strcmp(var->name, "pixnum") == 0){
					output->pixel_number = (unsigned)var->stats.center_value;
					output->pixel_number_var = var;
				} else if (strcmp(var->name, "l2_flags") == 0){
					output->l2_flags = (unsigned)var->stats.center_value;
				} else if (strcmp(var->name, "flags_sst") == 0){
					output->flags_sst = (unsigned)var->stats.center_value;
				} else if (strcmp(var->name, "flags_sst4") == 0){
					output->flags_sst4 = (unsigned)var->stats.center_value;
				}
			}
		}
		break;
	}
	return 0;
}

int sstval_extract_clean(sstval_extract_arguments *arguments){
	val_extract_clean(&arguments->val_extract_arguments);
	return 0;
}

int sstval_extract(sstval_extract_arguments *arguments) {
	if (!arguments){
		return SSTVALEXTRACT_ERR_INPUT;
	}
	arguments->sstval_extract_parser(SSTVALEXTRACT_KEY_INIT, NULL, arguments->user_input);

	int ret = 0;
	sstval_extract_output *output = calloc(1, sizeof(sstval_extract_output));
	if (output == NULL){
		arguments->sstval_extract_parser(SSTVALEXTRACT_KEY_ERROR, NULL, arguments->user_input);
		ret = SSTVALEXTRACT_ERR_UNKNOWN;
	} else {
		output->print_header = arguments->header;
		if (arguments->val_extract_arguments.ifile != blank_string_for_ifile){
			output->buoy_id = arguments->buoy_id;
			output->buoy_time = arguments->buoy_time;

			if (arguments->val_extract_arguments.box_size <= 0){
				arguments->val_extract_arguments.box_size = 3;
			}

			arguments->val_extract_arguments.geospatial_only = 0;
			arguments->val_extract_arguments.geospatial_to_double = 1;
			arguments->val_extract_arguments.val_extract_parser = (val_extract_parser)read_val_extract;
			arguments->val_extract_arguments.user_input = output;

			ret = val_extract(&arguments->val_extract_arguments);
		} else {
			arguments->val_extract_arguments.ifile = NULL;
		}

		if (ret){
			arguments->sstval_extract_parser(SSTVALEXTRACT_KEY_ERROR, NULL, arguments->user_input);
		} else {
			arguments->sstval_extract_parser(SSTVALEXTRACT_KEY_OUTPUT, output, arguments->user_input);
			arguments->sstval_extract_parser(SSTVALEXTRACT_KEY_SUCCESS, NULL, arguments->user_input);
		}
		if (output->granule){
			free(output->granule);
		}
		if (output->reflective_lut_version){
			free(output->reflective_lut_version);
		}
		free(output);
	}
	arguments->sstval_extract_parser(SSTVALEXTRACT_KEY_FINI, NULL, arguments->user_input);
	// I'm not sure I like this error code paradigm.  It requires the parent and child parsers to have
	// complete different error code ranges with zero overlap.
	if (ret >= SSTVALEXTRACT_ERR_POINT_NOT_FOUND && ret <= SSTVALEXTRACT_ERR_L2QC){
		return ret;
	} else if (ret == VALEXTRACT_ERR_UNKNOWN){
		return SSTVALEXTRACT_ERR_UNKNOWN;
	} else if (ret == VALEXTRACT_ERR_POINT_NOT_FOUND){
		return SSTVALEXTRACT_ERR_POINT_NOT_FOUND;
	} else if (ret == VALEXTRACT_ERR_NCFILE_ERR){
		return SSTVALEXTRACT_ERR_NCFILE_ERR;
	} else if (ret == VALEXTRACT_ERR_NCFILE_INVALID){
		return SSTVALEXTRACT_ERR_NCFILE_INVALID;
	} else if (ret == VALEXTRACT_ERR_FLAG){
		return SSTVALEXTRACT_ERR_FLAG;
	} else if (ret == VALEXTRACT_ERR_VARIABLE){
		return SSTVALEXTRACT_ERR_VARIABLE;
	} else if (ret == VALEXTRACT_ERR_INPUT){
		return SSTVALEXTRACT_ERR_INPUT;
	} else if (ret == VALEXTRACT_ERR_L2QC){
		return SSTVALEXTRACT_ERR_L2QC;
	} else if (ret != VALEXTRACT_ERR_NONE){
		eprintf("val-extract returned an unknown error code\n");
		return ret;
	}

	return ret;
}
