
#include "argpar.h"
#include "val_extract.h"

#include <errno.h>
#include <string.h>

#define OFILE_DEFAULT_EXT ".qc"

#ifdef DEBUG
#define dprintf(...) do { printf(__VA_ARGS__); } while(0)
#else
#define dprintf(...) do {} while (0)
#endif

#define eprintf(...) do { fprintf(stderr, __VA_ARGS__); } while(0)

#define str(s) #s
#define expanded_str(s) str(s)

typedef struct val_extract_main_input {
	char *ofile;
	unsigned printed_products;
	FILE *base_ofile_h;
	val_extract_arguments *val_extract_arguments;
} val_extract_main_input;

static argpar_option options[] = {
	{ "ofile", 'o', "FILE", 0, "base output file path" },
	{ "help", 'h', "BOOL", OPTION_HIDDEN, "display usage information" },
	{ "version", 'v', "BOOL", OPTION_HIDDEN, "display version information" },
	{ 0,0,0,0, "Output files:", -2 },
	{ 0,0,0, OPTION_DOC,
		"Any number of products may be specified.  If none are given, every product that appears to be "
		"geospatial, having scan line and pixel as the only dimensions, are processed, outputting statistics "
		"such as min, max, and mean." },
	{ 0,0,0, OPTION_DOC, "Two sets of output files are created, one base file describing the extract as a whole, "
		"located at <ofile> or, if not set, <ifile>" OFILE_DEFAULT_EXT ", and one file for each product requested, "
		"located at <base>.<product>.  If the file is ignored by the L2QC, the product files are omitted." },
	{ 0,0,0, 0, "Return values:", -1 },
		{ str(VALEXTRACT_ERR_NONE) "=" expanded_str(VALEXTRACT_ERR_NONE), 0, 0, OPTION_DOC,
			"Extract successfully processed" },
		{ str(VALEXTRACT_ERR_NCFILE_ERR) "=" expanded_str(VALEXTRACT_ERR_NCFILE_ERR), 0, 0, OPTION_DOC,
			"NetCDF file is malformed" },
		{ str(VALEXTRACT_ERR_NCFILE_INVALID) "=" expanded_str(VALEXTRACT_ERR_NCFILE_INVALID), 0, 0, OPTION_DOC,
			"NetCDF file is fine but isn't in the format expected (doesn't have geospatial dimensions, etc)" },
		{ str(VALEXTRACT_ERR_FLAG) "=" expanded_str(VALEXTRACT_ERR_FLAG), 0, 0, OPTION_DOC,
			"Error processing/finding flags" },
		{ str(VALEXTRACT_ERR_VARIABLE) "=" expanded_str(VALEXTRACT_ERR_VARIABLE), 0, 0, OPTION_DOC,
			"Error processing/finding a product" },
		{ str(VALEXTRACT_ERR_INPUT) "=" expanded_str(VALEXTRACT_ERR_INPUT), 0, 0, OPTION_DOC,
			"Bad or no input given" },
		{ str(VALEXTRACT_ERR_L2QC) "=" expanded_str(VALEXTRACT_ERR_L2QC), 0, 0, OPTION_DOC,
			"L2QC flag over threshold" },
		{ str(VALEXTRACT_ERR_UNKNOWN) "=" expanded_str(VALEXTRACT_ERR_UNKNOWN), 0, 0, OPTION_DOC,
			"Unknown error, such as a malloc failure or permissions problem." },
	{ 0 } // tell argpar to stop checking options
};

static const char doc[] =
	"val_extract is a utility designed to process a small section of a Level-2 "
	"NetCDF file given either a bounding box or a center point and a box size.  "
	"A product list can be given, as well, to limit processing targets.";

static const char args_doc[] = "ifile=<file> <box definition> [products...]";

const char *argpar_program_name = "val_extract";
const char *argpar_program_version = "1.2.2";

static int parse_options(int key, char *argv, struct argpar_state *state) {
	int ret = 0;

	val_extract_main_input *arguments = state->input;
	switch (key){
		case 'h':
			return ARGPAR_ERR_USAGE;
			break;
		case 'v':
			printf("val_extract=%s\nval_extract-main=%s\n", val_extract_version(), argpar_program_version);
			return ARGPAR_ERR_ABORT;
			break;
		case 'o':
			arguments->ofile = (void*)argv;
			break;
		case ARGPAR_KEY_ARG:
			arguments->val_extract_arguments->products[arguments->val_extract_arguments->product_count++] = argv;
			break;
		case ARGPAR_KEY_INIT:
			state->child_inputs[0] = arguments->val_extract_arguments;
			arguments->val_extract_arguments->product_count = 0;

			arguments->val_extract_arguments->products = malloc(state->argc * sizeof(char*));
			if (arguments->val_extract_arguments->products == NULL){
				return ARGPAR_ERR_UNKNOWN;
			}
			break;
	}
	return ret;
}

static const argpar_child argpar_children[] = {{&val_extract_argpar}, {0}};
static argpar val_extract_main_argpar = { options, parse_options, args_doc, doc, argpar_children };

static const char *TOTAL_PIXEL_STR = "pixel_count";
static const char *TOTAL_PIXEL_NOT_FLAGGED_STR = "unflagged_pixel_count";
static const char *TOTAL_PIXEL_FLAGGED_STR = "flagged_pixel_count";
static const char *ASCII_TIME_STR = "time";
static const char *PROD_NAME_STR = "name";
static const char *MAX_PIXEL_STR = "max";
static const char *MIN_PIXEL_STR = "min";
static const char *MEAN_PIXEL_STR = "mean";
static const char *STDDEV_STR = "stddev";
static const char *RMS_STR = "rms";
static const char *VALID_PIXEL_STR = "valid_pixel_count";
static const char *MEDIAN_PIXEL_STR = "median";
static const char *CENTER_PIXEL_STR = "center_value";

static void print_stats_to_file(FILE *h, const char *prefix, nc_var_stats *stats){
	fprintf(h, "%s%s=%d\n", prefix, VALID_PIXEL_STR, stats->count);
	fprintf(h, "%s%s=%f\n", prefix, MAX_PIXEL_STR, stats->max);
	fprintf(h, "%s%s=%f\n", prefix, MIN_PIXEL_STR, stats->min);
	fprintf(h, "%s%s=%f\n", prefix, MEAN_PIXEL_STR, stats->mean);
	fprintf(h, "%s%s=%f\n", prefix, MEDIAN_PIXEL_STR, stats->median);
	fprintf(h, "%s%s=%f\n", prefix, STDDEV_STR, stats->stddev);
	fprintf(h, "%s%s=%f\n", prefix, RMS_STR, stats->rms);
}

static int save_extract(int key, void *nc_input, void *user_input){
	val_extract_main_input *input = user_input;
	switch(key){
	case VALEXTRACT_KEY_INIT:
		break;
	case VALEXTRACT_KEY_SUCCESS:
		break;
	case VALEXTRACT_KEY_ERROR:
		break;
	case VALEXTRACT_KEY_FINI:
		if (input->base_ofile_h){
			fprintf(input->base_ofile_h, "\n");
			fprintf(input->base_ofile_h, "versions=val_extract=%s val_extract-main=%s\n", val_extract_version(), argpar_program_version);
			fclose(input->base_ofile_h);
		}
		break;
	case VALEXTRACT_KEY_FILE:
		if (nc_input){ // useless check for dumb compiler problem
			nc_region *region = nc_input;
			nc_file *file = region->file;
			FILE *ofile_h = NULL;
			if ((ofile_h = fopen(input->ofile, "w")) == NULL){
				eprintf("Error opening output file %s, %s\n", input->ofile, strerror(errno));
				return -1;
			}

			fprintf(ofile_h, "%s=%s\n", ASCII_TIME_STR, region->ascii_time);
			fprintf(ofile_h, "utime=%zd\n", region->utime);
			fprintf(ofile_h, "dim_length_lines=%d\n", file->dim_lengths[file->line_dimid]);
			fprintf(ofile_h, "dim_length_pixels=%d\n", file->dim_lengths[file->pixel_dimid]);
			if (input->val_extract_arguments->start_lat != VALEXTRACT_UNSET){
				fprintf(ofile_h, "slat=%f\n", input->val_extract_arguments->start_lat);
				fprintf(ofile_h, "slon=%f\n", input->val_extract_arguments->start_lon);
				if (!input->val_extract_arguments->box_size){
					fprintf(ofile_h, "elat=%f\n", input->val_extract_arguments->end_lat);
					fprintf(ofile_h, "elon=%f\n", input->val_extract_arguments->end_lon);
				}
			} else {
				fprintf(ofile_h, "sline=%d\n", input->val_extract_arguments->start_line);
				fprintf(ofile_h, "spixl=%d\n", input->val_extract_arguments->start_pixel);
				if (!input->val_extract_arguments->box_size){
					fprintf(ofile_h, "eline=%d\n", input->val_extract_arguments->end_line);
					fprintf(ofile_h, "epixl=%d\n", input->val_extract_arguments->end_pixel);
				}
			}
			if (input->val_extract_arguments->box_size){
				fprintf(ofile_h, "boxsize=%d\n", input->val_extract_arguments->box_size);
			}
			fprintf(ofile_h, "center_lat=%f\n", region->lat);
			fprintf(ofile_h, "center_lon=%f\n", region->lon);
			fprintf(ofile_h, "center_line=%zd\n", region->center.line);
			fprintf(ofile_h, "center_pixel=%zd\n", region->center.pixel);
			fprintf(ofile_h, "%s=%d\n", TOTAL_PIXEL_STR, region->pixel_count);
			fprintf(ofile_h, "%s=%d\n", TOTAL_PIXEL_NOT_FLAGGED_STR, region->unflagged_pixel_count);
			fprintf(ofile_h, "%s=%d\n", TOTAL_PIXEL_FLAGGED_STR, region->flagged_pixel_count);
			fprintf(ofile_h, "platform=%s\n", region->file->platform);
			fprintf(ofile_h, "instrument=%s\n", region->file->instrument);
			fprintf(ofile_h, "flag_counts=");
			int flag_i;
			for (flag_i=0;flag_i<file->flag_count;flag_i++){
				if (strcmp("SPARE", file->flag_meanings[flag_i])){
					fprintf(ofile_h, "%*s%s=%d", flag_i ? 1 : 0, "", file->flag_meanings[flag_i], region->flag_counts[flag_i]);
				}
			}
			fprintf(ofile_h, "\nflag_percents=");
			for (flag_i=0;flag_i<file->flag_count;flag_i++){
				if (strcmp("SPARE", file->flag_meanings[flag_i])){
					double flag_percent = (double)region->flag_counts[flag_i]/region->pixel_count;
					fprintf(ofile_h, "%*s%s=%f", flag_i ? 1 : 0, "", file->flag_meanings[flag_i], flag_percent);
				}
			}
			fprintf(ofile_h, "\nvariables=");
			input->base_ofile_h = ofile_h;
			input->printed_products = 0;
		}
		break;
	case VALEXTRACT_KEY_VAR:
		if (nc_input){ // useless check for dumb compiler problem
			nc_var *var = nc_input;
			if (!var->is_geospatial){
				return 0;
			}
			nc_file *file = var->file;
			char output_filename[255];
			strcpy(output_filename, input->ofile);
			strcat(output_filename, ".");
			strcat(output_filename, var->name);

			FILE *output_h = NULL;
			if ((output_h = fopen(output_filename, "w")) == NULL){
				eprintf("Error opening output file %s, %s\n", output_filename, strerror(errno));
				return -1;
			}

			if (var->is_geospatial){
				fprintf(output_h, "%s=%s\n", PROD_NAME_STR, var->name);
				if (var->units != NULL){
					fprintf(output_h, "units=%s\n", var->units);
				}
				fprintf(output_h, "%s=%f\n", CENTER_PIXEL_STR, var->stats.center_value);
				print_stats_to_file(output_h, "", &var->stats);
				print_stats_to_file(output_h, "filtered_", &var->filtered_stats);
				print_stats_to_file(output_h, "iqr_", &var->iqr_stats);
				fclose(output_h);
			} else {
				int i, j, nc_ret;

				int ndims = var->ndims;
				int dim_lengths[ndims];
				for (i=0;i<ndims;i++){
					dim_lengths[i] = file->dim_lengths[var->dim_ids[i]];
				}

				int value_count = dim_lengths[ndims-1];
				for (i=ndims-2;i>=0;i--){
					value_count *= dim_lengths[i];
				}

				double values[value_count];
				if ((nc_ret = nc_get_var_double(var->gid, var->varid, values)) != NC_NOERR){
					eprintf("Error getting entire dataset for %s, error %d\n", var->name, nc_ret);
					return -1;
				}
				FILE *output_h = NULL;
				if ((output_h = fopen(output_filename, "w")) == NULL){
					eprintf("Error opening output file %s, %s\n", output_filename, strerror(errno));
					return -1;
				}
				for (i=0;i<value_count;i++){
					fprintf(output_h, "value");
					int multidim_index[ndims];
					unflatten_index(i, ndims, dim_lengths, multidim_index);
					for (j=0;j<ndims;j++){
						fprintf(output_h, "[%d]", multidim_index[j]);
					}
					if (var->has_fill && values[i] == var->fill){
						fprintf(output_h, "=%s\n", "NULL");
					} else {
						double v = values[i];
						if (var->has_scale){
							v *= var->scale;
						}
						if (var->has_offset){
							v += var->offset;
						}
						fprintf(output_h, "=%f\n",v);
					}
				}
				fclose(output_h);
			}
			fprintf(input->base_ofile_h, "%s%s", input->printed_products++ ? "," : "", var->name);
		}
		break;
	}
	return 0;
}

int main(int argc, char *argv[]) {
	if (argc <= 1){
		argpar_usage_default(&val_extract_main_argpar);
		return VALEXTRACT_ERR_INPUT;
	}
	val_extract_main_input input = {0};
	val_extract_arguments arguments = {
		.geospatial_only = 1, .geospatial_to_double = 1,
		.val_extract_parser = (val_extract_parser)save_extract,
		.user_input = (void*)&input
	};
	input.val_extract_arguments = &arguments;

	int ret = EXIT_SUCCESS;
	if (argpar_parse_args(&val_extract_main_argpar, argc, argv, 0, NULL, &input)){
		ret = VALEXTRACT_ERR_INPUT;
	} else {
	    char output_filename[255];
	    printf("%s\n", arguments.ifile);
		if (!input.ofile){
			strcpy(output_filename, arguments.ifile);
			strcat(output_filename, OFILE_DEFAULT_EXT);
			input.ofile = output_filename;
		}
		ret = val_extract(&arguments);
	}
	val_extract_clean(&arguments);
	argpar_clean(&val_extract_main_argpar);
	return ret;
}
