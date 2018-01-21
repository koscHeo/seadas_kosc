#include "argpar.h"
#include "sstval_extract.h"

#include <float.h>
#include <limits.h>

#define str(s) #s
#define expanded_str(s) str(s)

#ifdef DEBUG
#define dprintf(...) do { printf(__VA_ARGS__); } while(0)
#else
#define dprintf(...) do {} while (0)
#endif

#define eprintf(...) do { fprintf(stderr, __VA_ARGS__); } while(0)

typedef struct sstval_extract_main_input {
	char *ofile;
	unsigned printed_products;
	FILE *base_ofile_h;
	sstval_extract_arguments *sstval_extract_arguments;
} sstval_extract_main_input;

static argpar_option options[] = {
	{ "help", 'h', "1", OPTION_HIDDEN, "display usage information", 0 },
	{ 0,0,0, OPTION_PARENT, "Return values:", -1 },
		{ str(SSTVALEXTRACT_ERR_NONE) "=" expanded_str(SSTVALEXTRACT_ERR_NONE), 0, 0, OPTION_DOC,
			"Extract successfully processed" },
		{ str(SSTVALEXTRACT_ERR_NCFILE_ERR) "=" expanded_str(SSTVALEXTRACT_ERR_NCFILE_ERR), 0, 0, OPTION_DOC,
			"NetCDF file is malformed" },
		{ str(SSTVALEXTRACT_ERR_NCFILE_INVALID) "=" expanded_str(SSTVALEXTRACT_ERR_NCFILE_INVALID), 0, 0, OPTION_DOC,
			"NetCDF file is fine but isn't in the format expected (doesn't have geophysical dimensions, etc)" },
		{ str(SSTVALEXTRACT_ERR_FLAG) "=" expanded_str(SSTVALEXTRACT_ERR_FLAG), 0, 0, OPTION_DOC,
			"Error processing/finding flags" },
		{ str(SSTVALEXTRACT_ERR_VARIABLE) "=" expanded_str(SSTVALEXTRACT_ERR_VARIABLE), 0, 0, OPTION_DOC,
			"Error processing/finding a product" },
		{ str(SSTVALEXTRACT_ERR_INPUT) "=" expanded_str(SSTVALEXTRACT_ERR_INPUT), 0, 0, OPTION_DOC,
			"Bad or no input given" },
		{ str(SSTVALEXTRACT_ERR_L2QC) "=" expanded_str(SSTVALEXTRACT_ERR_L2QC), 0, 0, OPTION_DOC,
			"L2QC flag over threshold" },
		{ str(SSTVALEXTRACT_ERR_UNKNOWN) "=" expanded_str(SSTVALEXTRACT_ERR_UNKNOWN), 0, 0, OPTION_DOC,
			"Unknown error, such as a malloc failure or permissions problem." },
	{ 0 } // tell argpar to stop checking options
};


static int parse_options(int key, char *argv, struct argpar_state *state) {
	int ret = 0;

	sstval_extract_main_input *arguments = state->input;
	switch (key){
	case 'h':
		return ARGPAR_ERR_USAGE;
	case ARGPAR_KEY_ARG:
		eprintf("No arguments are accepted\n");
		return ARGPAR_ERR_USAGE;
	case ARGPAR_KEY_INIT:
		state->child_inputs[0] = arguments->sstval_extract_arguments;
		break;
	}
	return ret;
}

static const char doc[] = "sstval_extract blah blah blah.";
static const char args_doc[] = "[header=1] [<request ID> [<request file>]]";

static const argpar_child argpar_children[] = {{&sstval_extract_argpar}, {0}};
static argpar sstval_extract_main_argpar = { options, parse_options, args_doc, doc, argpar_children };


const char *argpar_program_name = "sstval_extract";
const char *argpar_program_version = "1.0.0";

static int print_sstval(int key, void *sstval_input, void *user_input){
	switch(key){
	case SSTVALEXTRACT_KEY_INIT:
		break;
	case SSTVALEXTRACT_KEY_SUCCESS:
		break;
	case SSTVALEXTRACT_KEY_ERROR:
		break;
	case SSTVALEXTRACT_KEY_FINI:
		break;
	case SSTVALEXTRACT_KEY_OUTPUT:
		if (sstval_input){ // useless check for dumb compiler problem
			sstval_extract_output *output = (sstval_extract_output*)sstval_input;

			#define output_row() \
				output_string("bid", buoy_id); \
				output_string("btime", buoy_time); \
				output_uint("year", year); \
				output_uint("day", day_of_year); \
				output_float("shour", hour); \
				output_float("clon", lon); \
				output_float("clat", lat); \
				output_float("csolz", solz); \
				output_float("csenz", senz); \
				output_float("cazim", relaz); \
				output_float("caoi", angle_of_incidence); \
				output_uint("cmirror", mirror); \
				output_uint("cdetnum", detector_number); \
				output_uint("pixnum", pixel_number); \
				output_uint("cl2comflg", l2_flags); \
				output_static_int("cqual", -999); \
				output_float("cqsst", quality_sst); \
				output_float("cqsst4", quality_sst4); \
				output_uint("cglntf", glint_flag); \
				output_stats("sst", sst); \
				output_stats("sst4", sst4); \
				output_stats("20", bt_3750); \
				output_stats("22", bt_3959); \
				output_stats("23", bt_4050); \
				output_stats("31", bt_11000); \
				output_stats("32", bt_12000); \
				output_stats("26", rho_cirrus); \
				output_stats("27", bt_6715); \
				output_stats("28", bt_7325); \
				output_stats("29", bt_8550); \
				output_stats_no_c("cd2223ref", bt3959minus4050); \
				output_stats_no_c("ch14l", rhot_678); \
				output_static_int_stats_no_c("ch14h", -999); \
				output_stats("creysst", sstref); \
				output_char("sat", sensor); \
				output_string("granule", granule); \
				output_uint("csstflg", flags_sst); \
				output_uint("csst4flg", flags_sst4); \
				output_last_string("l1blutver", reflective_lut_version);

			if (output->print_header){
				#define output_string(col_name, name) fprintf(stdout, "%s,", col_name)
				#define output_last_string(col_name, name) fprintf(stdout, "%s\n", col_name)
				#define output_char(col_name, name) output_string(col_name, name)
				#define output_uint(col_name, name) output_string(col_name, name)
				#define output_float(col_name, name) output_string(col_name, name)
				#define output_static_int(col_name, i) fprintf(stdout, "%s,", col_name)
				#define output_static_float(col_name, f) fprintf(stdout, "%s,", col_name)
				#define output_stats(col_name, name) fprintf(stdout, "%s,%s,%s,%s,%s,%s,", "c" col_name, "min" col_name, "max" col_name, "avg" col_name, "med" col_name, "sd" col_name)
				#define output_stats_no_c(col_name, name) fprintf(stdout, "%s,%s,%s,%s,%s,%s,", col_name, "min" col_name, "max" col_name, "avg" col_name, "med" col_name, "sd" col_name)
				#define output_static_int_stats_no_c(col_name, name) output_stats_no_c(col_name, name)
				fprintf(stdout, "/begin_header\n/missing=-999\n/delimiter=comma\n/fields=");
				output_row();
				fprintf(stdout, "/request_id=%s\n", output->request_id ? output->request_id : "");
				fprintf(stdout, "/request_file=%s\n", output->request_file ? output->request_file : "");
				fprintf(stdout, "/end_header\n");
			}

			if (output->buoy_id && output->buoy_time){
				#undef output_string
				#undef output_last_string
				#undef output_char
				#undef output_uint
				#undef output_float
				#undef output_static_int
				#undef output_static_float
				#undef output_stats
				#undef output_stats_no_c
				#undef output_static_int_stats_no_c

				#define output_string(col_name, name) fprintf(stdout, "%s,", output->name ? output->name : "")
				#define output_last_string(col_name, name) fprintf(stdout, "%s\n", output->name ? output->name : "")
				#define output_char(col_name, name) if (output->name){fprintf(stdout, "%c,", output->name);} else {fprintf(stdout, ",");}
				#define output_uint(col_name, name) if (output->name != UINT_MAX){fprintf(stdout, "%u,", output->name);} else {fprintf(stdout, "-999,");}
				#define output_float(col_name, name) if (output->name != FLT_MAX){fprintf(stdout, "%f,", output->name);} else {fprintf(stdout, "-999,");}
				#define output_static_int(col_name, i) fprintf(stdout, "%d,", i)
				#define output_static_float(col_name, f) fprintf(stdout, "%f,", f)
				#define output_stats(col_name, name) \
					if (output->name.initialized){ \
						fprintf(stdout, "%f,%f,%f,%f,%f,%f,", output->name.center_value, output->name.min, output->name.max, output->name.mean, output->name.median, output->name.stddev); \
					} else { \
						fprintf(stdout, "-999,-999,-999,-999,-999,-999,"); \
					}
				#define output_stats_no_c(col_name, name) output_stats(col_name, name)
				#define output_static_int_stats_no_c(col_name, name) fprintf(stdout, "-999,-999,-999,-999,-999,-999,")

				output_row();
			}
		}
		break;
	}
	return 0;
}

int main(int argc, char *argv[]) {
	if (argc <= 1){
		argpar_usage_default(&sstval_extract_main_argpar);
		return SSTVALEXTRACT_ERR_INPUT;
	}
	sstval_extract_main_input input = {0};
	sstval_extract_arguments arguments = {
		.sstval_extract_parser = (sstval_extract_parser)print_sstval,
		.user_input = &input
	};
	input.sstval_extract_arguments = &arguments;
	int ret = EXIT_SUCCESS;
	if (argpar_parse_args(&sstval_extract_main_argpar, argc, argv, 0, NULL, &input)){
		ret = SSTVALEXTRACT_ERR_INPUT;
	} else {
		ret = sstval_extract(&arguments);
	}
	sstval_extract_clean(&arguments);
	argpar_clean(&sstval_extract_main_argpar);
	return ret;
}
