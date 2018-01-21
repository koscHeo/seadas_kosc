
// implemented version numbers
#define ARGPAR_IMPLEMENTATION 2000002
#define ARGPAR_IMPLEMENTATION_STR "2.0.2"
#define ARGPAR_API_VERSION 2000002
#define ARGPAR_API_VERSION_STR "2.0.2"
#include "argpar.h"

#include <ctype.h>
#include <errno.h>
#include <libgen.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct parser;
typedef struct parser parser;
struct parser_chain;
typedef struct parser_chain parser_chain;

struct parser {
	argpar_parser parser;
	argpar *argpar;
	unsigned args_parsed;
	parser *parent;
	unsigned parent_index;
	void *input, **child_inputs;
	void *hook;
};

struct parser_chain {
	argpar *root_argpar;
	parser *parsers;
	parser *last_parser;
	void **child_inputs;
	void *last_child_input;
	argpar_state state;
	void *storage;
};

typedef struct parser_counts {
	size_t parsers, child_inputs;
} parser_counts;

static long strtol_strict(const char *s){
	char *endptr;
	errno = 0;
	long ret = strtol(s, &endptr, 0);
	if (!errno && *endptr != '\0'){
		errno = 1;
	}
	return ret;
}

static double strtod_strict(const char *s){
	char *endptr;
	errno = 0;
	double ret = strtod(s, &endptr);
	if (!errno && *endptr != '\0'){
		errno = 1;
	}
	return ret;
}
static void count_argpars(argpar *argpar, parser_counts *counts){
	if (argpar->options || argpar->parser) {
		counts->parsers++;
	}
	const argpar_child *children = argpar->children;
	if (children){
		while (children->argpar) {
			count_argpars((children++)->argpar, counts);
			counts->child_inputs++;
		}
	}
}

static void add_parser(parser_chain *parser_chain, argpar *argpar, parser *parent, unsigned parent_index){
	parser *p = parser_chain->last_parser++;
	p->parser = argpar->parser;
	p->argpar = argpar;
	p->args_parsed = 0;
	p->parent = parent;
	p->parent_index = parent_index;
	p->input = NULL;
	p->hook = NULL;

	const argpar_child *children = argpar->children;
	if (children){
		p->child_inputs = parser_chain->last_child_input;
		unsigned child_count = 0;
		while (children++->argpar) {
			child_count++;
		}
		parser_chain->last_child_input += child_count;

		children = argpar->children;
		unsigned child_index = 0;
		while (children->argpar) {
			add_parser(parser_chain, children++->argpar, p, child_index++);
		}
	}
}
static int call_parser(parser *p, argpar_state *state, int key, char *arg){
	if (p->parser) {
		int err;
		state->hook = p->hook;
		state->input = p->input;
		state->child_inputs = p->child_inputs;
		state->arg_num = p->args_parsed;
		err = (*p->parser)(key, arg, state);
		p->hook = state->hook;
		return err;
	} else {
		return ARGPAR_ERR_UNKNOWN;
	}
}

static int parser_chain_init(parser_chain *parser_chain, argpar *argpar, unsigned argc, char **argv, unsigned flags, void *input) {
	parser_counts counts = {0, 0};

	if (argpar){
		count_argpars(argpar, &counts);
	}

	size_t parser_size = (counts.parsers + 1) * sizeof(parser);
	size_t input_size = counts.child_inputs * sizeof(void*);
	parser_chain->storage = malloc(parser_size + input_size);
	if (!parser_chain->storage){
		return ENOMEM;
	}
	parser_chain->root_argpar = argpar;
	parser_chain->parsers = parser_chain->storage;
	parser_chain->last_parser = parser_chain->parsers;
	parser_chain->child_inputs = parser_chain->storage + parser_size;
	parser_chain->last_child_input = parser_chain->child_inputs;

	memset(parser_chain->child_inputs, 0, input_size);
	memset(&parser_chain->state, 0, sizeof(argpar_state));

	parser_chain->state.argpar = parser_chain->root_argpar;
	parser_chain->state.argc = argc;
	parser_chain->state.argv = argv;
	parser_chain->state.flags = flags;
//	parser_chain->state.err_stream = stderr;
//	parser_chain->state.out_stream = stdout;
	parser_chain->state.next = 0;
//	parser_chain->state.pstate = parser;

	add_parser(parser_chain, argpar, NULL, 0);

	if (parser_chain->parsers < parser_chain->last_parser){
		parser_chain->parsers->input = input;
	}
	int err = 0;
	parser *p;
	for (p = parser_chain->parsers; p < parser_chain->last_parser && (!err || err == ARGPAR_ERR_UNKNOWN); p++) {
		if (p->parent != NULL){
			p->input = p->parent->child_inputs[p->parent_index];
		}
		if (!p->parser && p->argpar->children && p->argpar->children->argpar){
			// without a parser, every child gets the parent's input
//			p->child_inputs[0] = p->input;
			const argpar_child *children = p->argpar->children;
			unsigned child_index = 0;
			while (children++->argpar) {
				p->child_inputs[child_index++] = p->input;
			}
		}
		err = call_parser(p, &parser_chain->state, ARGPAR_KEY_INIT, 0);
	}
	if (err == ARGPAR_ERR_UNKNOWN){
		err = 0;
	}

	if (err){
		return err;
	}

	return 0;
}
static int parser_chain_finalize(parser_chain *p_chain, int err, unsigned *end_index) {
	struct parser *parser;

	if (!err) {
		if (p_chain->state.next == p_chain->state.argc){
			for (parser = p_chain->parsers; parser < p_chain->last_parser && (!err || err == ARGPAR_ERR_UNKNOWN); parser++){
				if (parser->args_parsed == 0){
					err = call_parser(parser, &p_chain->state, ARGPAR_KEY_NO_ARGS, 0);
				}
			}
			for (parser = p_chain->last_parser - 1; parser >= p_chain->parsers && (!err || err == ARGPAR_ERR_UNKNOWN); parser--){
				err = call_parser(parser, &p_chain->state, ARGPAR_KEY_END, 0);
			}

			if (err == ARGPAR_ERR_UNKNOWN){
				err = 0;
			}
			if (end_index){
				*end_index = p_chain->state.next;
			}
		} else if (end_index){
			*end_index = p_chain->state.next;
		} else {
//			if (!(p_chain->state.flags & ARGP_NO_ERRS) && p_chain->state.err_stream){
//				fprintf(p_chain->state.err_stream, dgettext(p_chain->argp->argp_domain, "%s: Too many arguments\n"), p_chain->state.name);
//			}
			fprintf(stderr, "%s: Too many arguments\n", p_chain->state.name);
			err = ARGPAR_ERR_UNKNOWN;
		}
	}

	if (err) {
		if (err == ARGPAR_ERR_USAGE){
//			__argp_state_help(&p_chain->state, p_chain->state.err_stream, ARGP_HELP_STD_ERR);
			argpar_usage(&p_chain->state);
		}

		for (parser = p_chain->parsers; parser < p_chain->last_parser; parser++){
			call_parser(parser, &p_chain->state, ARGPAR_KEY_ERROR, 0);
		}
	} else {
		for (parser = p_chain->last_parser - 1; parser >= p_chain->parsers && (!err || err == ARGPAR_ERR_UNKNOWN); parser--){
			err = call_parser(parser, &p_chain->state, ARGPAR_KEY_SUCCESS, 0);
		}
		if (err == ARGPAR_ERR_UNKNOWN){
			err = 0;
		}
	}

	for (parser = p_chain->last_parser - 1; parser >= p_chain->parsers; parser--){
		call_parser(parser, &p_chain->state, ARGPAR_KEY_FINI, 0);
	}

	free(p_chain->storage);

	return err;
}

static int parse_option(parser_chain *parser_chain, char *key, char *value){
	parser *parser;
	for (parser = parser_chain->parsers; parser < parser_chain->last_parser; parser++){
		const argpar_option *o = parser->argpar->options;
		bool is_parent = (parser->argpar == parser_chain->root_argpar);

		int i = -1;

		while (o[++i].name || o[i].doc){
			if (o[i].name && !((o[i].flags & OPTION_PARENT) && !is_parent) && !((o[i].flags & OPTION_CHILD) && is_parent) && !strcmp(o[i].name, key)){
				argpar_state s = parser_chain->state;
				s.arg_name = key;
				s.arg_value = value;

				if (o[i].flags & OPTION_INT){
					s.argv_as_int = strtol_strict(value);
					s.argv_as_int_err = errno;
				} else {
					s.argv_as_int = 0;
					s.argv_as_int_err = 1;
				}
				if (o[i].flags & OPTION_DBL){
					s.argv_as_dbl = strtod_strict(value);
					s.argv_as_dbl_err = errno;
				} else {
					s.argv_as_dbl = 0;
					s.argv_as_dbl_err = 1;
				}
				return call_parser(parser, &s, o[i].key, value);
			}
		}
	}
	return ARGPAR_ERR_UNKNOWN;
}

static int parse_arg(parser_chain *parser_chain, char *arg){
	parser *parser;
	int err = ARGPAR_ERR_UNKNOWN;

	argpar_state s = parser_chain->state;
	if (s.flags & ARGPAR_CAST_ARGS){
		s.argv_as_int = strtol_strict(arg);
		s.argv_as_int_err = errno;
		s.argv_as_dbl = strtod_strict(arg);
		s.argv_as_dbl_err = errno;
	} else {
		s.argv_as_int = 0;
		s.argv_as_int_err = 1;
		s.argv_as_dbl = 0;
		s.argv_as_dbl_err = 1;
	}

	for (parser = parser_chain->parsers; parser < parser_chain->last_parser && err == ARGPAR_ERR_UNKNOWN; parser++){
		err = call_parser(parser, &s, ARGPAR_KEY_ARG, arg);
	}
	s.next++;
	return err;
}

static int parse_file(parser_chain *p_chain, const char *path, const char *cwd){
	argpar *p = NULL;

	parser *parser;
	for (parser = p_chain->parsers; parser < p_chain->last_parser && !p; parser++){
		p = parser->argpar;
	}
	if (!p){
		return ARGPAR_ERR_PARFILE;
	}
	if (p->parfile_buffer_count >= MAX_PARFILES){
		return ARGPAR_LIMIT_REACHED;
	}

	FILE *file_h = NULL;
	char *this_path, *this_cwd;
	this_path = this_cwd = NULL;

	if (cwd != NULL && *path != '/'){
		this_path = malloc((strlen(cwd)+strlen(path)+2)*sizeof(char));
		if (this_path == NULL){
			return ENOMEM;
		}
		strcpy(this_path, cwd);
		strcat(this_path, "/");
		strcat(this_path, path);
		file_h = fopen(this_path, "rb");
	} else {
		file_h = fopen(path, "rb");
	}
	if (file_h == NULL){
		if (this_path != NULL){
			free(this_path);
		}
		return ARGPAR_ERR_PARFILE;
	}

	fseek(file_h, 0, SEEK_END);
	long file_size = ftell(file_h);
	fseek(file_h, 0, SEEK_SET);
	unsigned long this_buffer_index = p->parfile_buffer_count++;
	if (p->parfile_buffer == NULL){
		p->parfile_buffer = calloc(MAX_PARFILES, sizeof(char*));
		if (p->parfile_buffer == NULL){
			if (this_path != NULL){
				free(this_path);
			}
			return ENOMEM;
		}
	}
	size_t ufile_size = (size_t)file_size;
	p->parfile_buffer[this_buffer_index] = malloc((ufile_size + 1) * sizeof(char));
	if (p->parfile_buffer[this_buffer_index] == NULL){
		if (this_path != NULL){
			free(this_path);
		}
		return ENOMEM;
	}
	size_t fread_count = fread(p->parfile_buffer[this_buffer_index], sizeof(char), ufile_size+1, file_h);
	if (fread_count != ufile_size || !feof(file_h)){
		if (this_path != NULL){
			free(this_path);
		}
		fclose(file_h);
		return ARGPAR_ERR_PARFILE;
	}
	fclose(file_h);

	char *buffer = p->parfile_buffer[this_buffer_index];
	char *buffer_end = &p->parfile_buffer[this_buffer_index][file_size];
	char *v, *k;
	int err = 0;
	do {
		k = strsep(&buffer, "=");
		v = strsep(&buffer, "\n");
		if (k == NULL || v == NULL){
			break;
		}
		while(isspace(*k)) k++;
		char *k_end = v-2;
		while(isspace(*k_end)){
			*k_end = '\0';
			k_end--;
		}
		while(isspace(*v)) v++;
		if (*k == '#'){
			continue;
		} else if (strlen(k)){
			if (strcmp(k, PARFILE_STR)){
				err = parse_option(p_chain, k, v);
				if (err == ARGPAR_ERR_UNKNOWN){
					fprintf(argpar_ostream ? argpar_ostream : stderr, "Unknown option: %s\n\n", k);
				}
			} else {
				if (this_cwd == NULL){
					if (this_path == NULL){
						this_path = strdup(path);
						if (this_path == NULL){
							return ENOMEM;
						}
					}
					this_cwd = dirname(this_path);
				}
				err = parse_file(p_chain, v, this_cwd);
			}
		}
	} while (buffer != buffer_end && !err);
	if (this_path){
		free(this_path);
	}
	return err;
}

int argpar_clean(argpar *p){
	if (p){
		unsigned i;
		for (i=0;i<p->parfile_buffer_count;i++){
			free(p->parfile_buffer[i]);
		}
		free(p->parfile_buffer);
		p->parfile_buffer_count = 0;

		// in case parsers manually loaded parfiles
		const argpar_child *children = p->children;
		if (children){
			while (children++->argpar) {
				argpar_clean(children->argpar);
			}
		}
	}

	return 0;
}

// stolen directly from argp
static char *nondestructive_basename(char *name){
	char *short_name = strrchr(name, '/');
	return short_name ? short_name + 1 : name;
}

int argpar_parse_args(argpar *p, unsigned argc, char *argv[], unsigned flags, unsigned *end_index, void *input){
	if (argc < 1 || argv == NULL){
		argpar_help(p, stdout, flags, (argv == NULL ? "cmd" : nondestructive_basename(argv[0])));
		return ARGPAR_ERR_USAGE;
	} else {
		parser_chain p_chain;
		if (p->parfile_buffer == NULL){
			p->parfile_buffer_count = 0;
		}
		int err = parser_chain_init(&p_chain, p, argc, argv, flags, input);
		if (!err){
			for (p_chain.state.next=1;p_chain.state.next<argc && !err;p_chain.state.next++){
				if (p_chain.state.quoted && p_chain.state.next < p_chain.state.quoted){
					p_chain.state.quoted = 0;
				}
				if (p_chain.state.quoted){
					err = parse_arg(&p_chain, argv[p_chain.state.next]);
					p_chain.state.arg_num++;
				} else if (strcmp(argv[p_chain.state.next], "--") == 0){
					p_chain.state.quoted = p_chain.state.next + 1;
				} else if (strchr(argv[p_chain.state.next], '=') != NULL){
					char *v = (char*)argv[p_chain.state.next];
					char *k = strsep(&v, "=");
					while(isspace(*k)) k++;
					if (strcmp(k, PARFILE_STR)){
						err = parse_option(&p_chain, k, v);
						if (err == ARGPAR_ERR_UNKNOWN){
							fprintf(argpar_ostream ? argpar_ostream : stderr, "Unknown option: %s\n\n", k);
						}
					} else {
						err = parse_file(&p_chain, v, NULL);
					}
				} else {
					err = parse_arg(&p_chain, argv[p_chain.state.next]);
					p_chain.state.arg_num++;
				}
			}
			err = parser_chain_finalize(&p_chain, err, end_index);
		}
		return err;
	}
}

int argpar_parse_file(argpar *p, const char *path, unsigned flags, void *input){
	parser_chain p_chain;
	int err = parser_chain_init(&p_chain, p, 0, NULL, flags, input);
	if (!err){
		err = parse_file(&p_chain, path, NULL);
		err = parser_chain_finalize(&p_chain, err, NULL);
	}
	return err;
}

const char *argpar_version(){
	return "argpar " ARGPAR_IMPLEMENTATION_STR ", API " ARGPAR_API_VERSION_STR;
}
