/** @file argpar.h
	@brief Library for reading command-line arguments in the form of key=value.

	This library parses a pointer of null-terminated strings (usually argv) and
	may modify the array in the process.  All pairs are stored for later access.
	Any argument that does not contain an equal sign is considered a "key argument"
	and is stored as-is in a non-associative array in the order in which they are
	found.  An argument of "--" will cause argpar to treat all remaining arguments
	as key arguments.

	Files containing newline separated key=value entries (called a par file) can also
	be parsed, triggered by a special key.  Any argument matching parfile=file
	(configurable via the PARFILE_STR macro) will cause the parser to immediately
	open the file (relative to the current working directory) and parse through the
	entire file.  Par files can also contain parfile=file entries, which will
	interrupt the current file parsing; these entries are relative to the current
	file being parsed instead of the current working directory, unless specified as
	absolute paths.  All repeat keys will overwrite previous entries.

	This library can also generate a usage statement based on the expected arguments.

	This library should not be used in conjunction with argp as some macros may have
	been re-used in an attempt to make this as much like argp as possible.

	argpar is not expected to free any strings passed in as arguments and they must
	remain in-memory and accessible until the argpar object is destroyed.
*/

#ifndef __ARGPAR_H__
#define __ARGPAR_H__

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>


#ifndef ARGPAR_API_VERSION
/** @brief Current API version.  If already defined, use API defined first. */
#define ARGPAR_API_VERSION 2000002
/** @brief Current API version as a string.  If already defined, use API defined first. */
#define ARGPAR_API_VERSION_STR "2.0.2"
#endif

/** @brief Maximum amount of parfiles to be able to load into one object. */
#ifndef MAX_PARFILES
#define MAX_PARFILES 64
#endif

/** @brief Used to override the par file trigger. */
#ifndef PARFILE_STR
#define PARFILE_STR "parfile"
#endif

extern const char *argpar_program_name;
extern FILE *argpar_ostream;

/** @typedef typedef struct argpar_option
	@brief For the lazy. */
/** @struct argpar_option
	@brief Stores the configuration for a par argument.

	The easiest way to use this is to create a static array at the top of your program
	to be passed to the argpar constructor, as shown below.

	If the first four fields are 0, the option is considered a new section header and is
	only used for the automated usage statement.

	Terminate the list with an empty struct.

	@code
		static argpar_option options[] = {
			{ "ifile", 'i', "FILE", 0, "input file" },
			{ "ofile", 'o', "FILE", 0, "output file" },
			{ 0,0,0,0, "This is another section:" },
			{ "efile", 'e', "FILE", 0, "where to output stderr" },
			{ 0 } // tell argpar to stop checking options
		};
	@endcode

*/
typedef struct argpar_option {
	/** @brief The key to search for. */
	const char *name;
	/** @brief The value to pass to the parser callback.  Must be a positive and non-zero. */
	int key;
	/** @brief What to display as the right side of the argument in the usage statement.

		If non-zero, this is the name of an argument associated with this option, which must
		be provided.  If zero, it will be automatically determined based on the desired type,
		which may be modified with the flag argument.
	*/
	const char *arg;
	/** @brief Modify the behavior of the argument by OR'ing one or more OPTION_ macros. */
	int flags;
	/** @brief A (generally short) description of an option.  In case of a documentation option,
			this may be as lengthy as required and will be wrapped to a reasonable size. */
	const char *doc;
	/** @brief The group associated with this option.  If 0 and a documentation option, it is
			auto-incremented.

			The group number is used for display order in the usage statement.  The order is
			0, 1, 2, ..., n, -m, ..., -2, -1 and options in each group are displayed in the
			order they are given in the options specified.
	*/
	int group;
} argpar_option;

/** @brief Passed to the argpar_parse_ functions, this tells the parser to cast all key arguments
		to long and double, as if OPTION_INT|OPTION_DBL were given to each. */
#define ARGPAR_CAST_ARGS 0x1

/** @brief Don't throw an error if not given and don't display at the top of the usage summary
		in the command line example. (Error not currently implemented.) */
#define OPTION_ARG_OPTIONAL 0x1

/** @brief Do not display this option anywhere in the usage summary. */
#define OPTION_HIDDEN 0x2

/** @brief This option isn't actually an option, merely text for the usage summary. */
#define OPTION_DOC 0x4

/** @brief Cast this option as a double.  The value and any error will be reflected in the
 	 	 argpar_state struct during the parser callback. */
#define OPTION_DBL 0x10

/** @brief Cast this option as a long.  The value and any error will be reflected in the
 	 	 argpar_state struct during the parser callback. */
#define OPTION_INT 0x20

/** @brief This option only applies if this parser is the top-most parser. Useful for return
 	 	 value documentation. */
#define OPTION_PARENT 0x40

/** @brief This option only applies if this parser is not the top-most parser. */
#define OPTION_CHILD 0x80

/** @brief Do not add an extra newline after this documentation string. Useful for lists and
 	 	 manual formatting. */
#define OPTION_DOC_NO_BREAK 0x100

/** @brief How to display the program name.  If not given, it will be derived during a call
		to argpar_parse_args. */
const char *argpar_program_name;

/** @brief Where to print errors and usage summaries to, defaults to STDOUT. */
FILE *argpar_ostream;

/** @typedef typedef struct argpar_state
	@brief For the lazy. */
/** @struct argpar_state
	@brief State variable to be filled before each call to the parser callback. */
typedef struct argpar_state {
	/** @brief Pointer to the parent argpar object. */
	struct argpar *argpar;

	/** @brief Total number of arguments being parsed by a call of argpar_parse_args. */
	unsigned argc;
	/** @brief Pointer to arguments passed to argpar_parse_args. */
	char **argv;

	/** @brief The index in ARGV of the next arg that to be parsed.  May be modified. */
	unsigned next;

	/** @brief Flags passed to parser functions to modify the behavior. */
	unsigned int flags;

	/** @brief The number of key arguments processed so far. */
	unsigned arg_num;

	/** @brief The index of the first argument after --, if found. */
	unsigned quoted;

	/** @brief Arbitrary pointer given by the user, generally to store parsed arguments. */
	void *input;

	/** @brief Values to pass to child parsers, same length as the number of child parsers
			for the current parser. */
	void **child_inputs;

	/** @brief For the parser's use.  Initialized to NULL. */
	void *hook;

	/** @brief The name used when printing messages.  This is initialized to ARGV[0],
			or PROGRAM_INVOCATION_NAME if that is unavailable. */
	char *name;

	/** @brief key part of key=value pair. */
	const char *arg_name;

	/** @brief value part of key=value pair. */
	const char *arg_value;

	/** @brief value part of key=value pair, parsed as a long, if OPTION_INT was used. */
	long argv_as_int;
	/** @brief Error state of integer casting, if OPTION_INT was used; 1 otherwise. */
	int argv_as_int_err;

	/** @brief value part of key=value pair, parsed as a double, if OPTION_DBL was used. */
	double argv_as_dbl;
	/** @brief Error state of double casting, if OPTION_DBL was used; 1 otherwise. */
	int argv_as_dbl_err;
} argpar_state;

/** @brief This is not an option at all, but rather a command line argument.  If a
   parser receiving this key returns success, the fact is recorded, and the
   ARGP_KEY_NO_ARGS case won't be used.  HOWEVER, if while processing the
   argument, a parser function decrements the NEXT field of the state it is
   passed, the option won't be considered processed; this is to allow you to
   actually modify the argument (perhaps into an option), and have it
   processed again.  */
#define ARGPAR_KEY_ARG		0

/** @brief There are remaining arguments not parsed by any parser, which may be found
   starting at (STATE->argv + STATE->next).  If success is returned, but
   STATE->next left untouched, it's assumed that all arguments were consume,
   otherwise, the parser should adjust STATE->next to reflect any arguments
   consumed.  */
#define ARGPAR_KEY_ARGS		-8

/** @brief Passed as the key to the parser callback function when there are no more
		arguments left to parse. */
#define ARGPAR_KEY_END -1

/** @brief Passed as the key to each parser callback function before any parsing
		occurs.  For most cases, this is the proper time to initialize child parsers'
		input pointers.  This is guaranteed to be called on every parser. */
#define ARGPAR_KEY_INIT -3

/** @brief Passed as the key to each parser callback function at the very end of
 	 	 the process.  This is guaranteed to be called on every parser. */
#define ARGPAR_KEY_FINI -4

/** @brief Passed as the key to each parser callback function when parsing
		has finished and no errors were returned. */
#define ARGPAR_KEY_SUCCESS -5

/** @brief Passed as the key to each parser callback function when parsing
		has finished and an error was returned. */
#define ARGPAR_KEY_ERROR -6

/** @brief Passed as the key to each parser callback function when parsing
		has finished and no key arguments were found. */
#define ARGPAR_KEY_NO_ARGS -7

/** @brief Passed as the key to each help filter when printing the args_doc. */
#define ARGPAR_KEY_HELP_ARGS_DOC -1

/** @brief Passed as the key to each help filter when printing the doc string before
		options are listed.  The text given is the portion preceding any vertical tabs
		('\v'), if any, or the entire string otherwise. */
#define ARGPAR_KEY_HELP_PRE_DOC -2

/** @brief Passed as the key to each help filter when printing the doc string after all
 	 	 options are listed, if a vertical tab ('\v') was found in doc.  Only the text
 	 	 after the tab is given to the filter. */
#define ARGPAR_KEY_HELP_POST_DOC -3

/** @brief Passed as the key to each help filter when printing a group header. */
#define ARGPAR_KEY_HELP_HEADER -4

/** @brief Passed as the key to each help filter when printing a group's trailing documentation. */
#define ARGPAR_KEY_HELP_OPTION_DOC -5

/** @brief Passed as the key to each help filter after all other documentation is printed, to
 	 	 return extra information not contained in the options list. */
#define ARGPAR_KEY_HELP_EXTRA -7 /* After all other documentation; text is NULL for this key. */

/** @brief Returned from the parser callback to signal argpar to stop parsing and return to
		the caller. */
#define ARGPAR_ERR_ABORT -1

/** @brief Returned from the parser callback to signal argpar to stop parsing, print a usage
		summary, and return to the caller. */
#define ARGPAR_ERR_USAGE -2

/** @brief What to return for unrecognized keys within an argpar_parser function. */
#define ARGPAR_ERR_UNKNOWN -5

/** @brief Returned from argpar_parse_args and argpar_parse_file when an error has occurred
		in reading a parfile. */
#define ARGPAR_ERR_PARFILE -3

/** @brief Returned from argpar_parse_args and argpar_parse_file when attempting to load more
		than the maximum amount of parfiles. */
#define ARGPAR_LIMIT_REACHED -4

/** @brief Pointer to a callback function to call for each argument and option parsed.

	@param[in] key argpar_option.key of recognized options or one of the ARGPAR_KEY_ macros.
	@param[in] argv value part of key=value parameter, or the entire argument in case of a key
		argument.
	@param[in] state Structure containing all the information needed to process an argument.

	@return 0 to continue or any of the ARGPAR_ERR_ macros to fail.
*/
typedef int (*argpar_parser) (int key, char *argv, argpar_state *state);

/** @brief Pointer to a callback function to call for each argument parsed.

	@param[in] key argpar_option.key of recognized options or one of the ARGPAR_KEY_HELP_ macros.
	@param[in] text The text about to be printed.  May be NULL for some ARGPAR_KEY_HELP_ keys.
	@param[in] input Currently unused, I think.

	@return NULL to print nothing, the input text pointer to print it unchanged, or a malloc'd
		string to print which is free'd by the parser after use.
*/
typedef char *(*help_filter) (int key, const char *text, void *input);


/** @typedef typedef struct argpar_child
	@brief For the lazy. */
/** @struct argpar_child
	@brief Child parser for nesting argpars. */
typedef struct argpar_child {
	/** @brief Underlying argpar struct. */
	struct argpar *argpar;

	/** @brief This should probably be the flags passed to the parser and help functions, but
	 	 	 is currently unused. Argp uses it to determine which sections to print. */
	unsigned flags;

	/** @brief If non-NULL, don't merge this child's options with the the parent's. If header
	 	 	 length is non-zero, it is printed before the child's options. */
	const char *header;

	/** @brief The group after which to print the child's arguments. */
	int group;
} argpar_child;


/** @typedef typedef struct argpar
	@brief For the lazy. */
/** @struct argpar
	@brief Master structure containing options, document strings, child parsers, and text
 	 	 filters. Everything needed to parse options and arguments and print help documentation. */
typedef struct argpar {
	/** @brief Array of option structures, terminated by an empty entry (IE: {0}). */
	const argpar_option *options;
	/** @brief Parser function handed all options and arguments, as well as various ARGPAR_KEY_
	 	 	 keys for events such as initialization, successful parsing, etc. */
	argpar_parser parser;
	/** @brief A string defining the desired arguments to the program, printed as part of the
	 	 	 usage documentation. */
	const char *args_doc;
	/** @brief A string printed before and/or after the options.  To display text after options,
	 	 	 use a vertical-tab (\v) to separate the text. */
	const char *doc;
	/** @brief Array of argpar_child structures containing parsers to nest within this one. All
	 	 	 options and arguments are passed to parent parsers first, then to children in the
	 	 	 order in which they're listed. */
	const argpar_child *children;
	/** @brief A function passed all documentation to print. */
	help_filter help_filter;
	/** @brief Used internally to store the number of parfiles parsed. Initialize to zero. */
	unsigned parfile_buffer_count;
	/** @brief Used internally to store the parfiles. Initialize to NULL. */
	char **parfile_buffer;
} argpar;

/** @brief Free any space consumed by argpar for parfiles.

	@param[in] p argpar structure to clean up.

	@return 0 on success.
*/
int argpar_clean(argpar *p);

/** @brief Parse a key=value store file.

	@param[in] p argpar object to use to parse file
	@param[in] path Path of file, absolute or relative to the current directory.
	@param[in] flags Flags to modify the behavior of the parser.
	@param[in] input Pointer to pass to the configured parser callback.

	@return 0 when the parsing is completed or the return value of the callback function
		if it bailed out early.

	@see ARGPAR_CAST_ARGS
*/
int argpar_parse_file(argpar *p, const char *path, unsigned flags, void *input);

/** @brief Parse an array of key=value pairs and/or key arguments.

	@param[in] p argpar object to use to parse file
	@param[in] argc Number of arguments to parse.
	@param[in] argv Array of key=value pairs and/or key arguments, starting at index 1.
		The first element is assumed to be the calling program's name.
	@param[in] flags Flags to modify the behavior of the parser.
	@param[out] end_index Set to the last index parsed.
	@param[in] input Pointer to pass to the configured parser callback.

	@return 0 when the parsing is completed or the return value of the callback function
		if it bailed out early.

	@see ARGPAR_CAST_ARGS
*/
int argpar_parse_args(argpar *p, unsigned argc, char *argv[], unsigned flags, unsigned *end_index, void *input);

/** @brief Print usage summary, called within a argpar_parser.

	@param[in] state argpar_state from which to generate usage summary.

	@return 0 on success, non-0 on error.
*/
int argpar_usage(argpar_state *state);

/** @brief Print the default usage summary with all available sections.

	@param[in] argpar argpar object from which to generate the usage summary.

	@return 0 on success, non-0 on error.
*/
int argpar_usage_default(argpar *argpar);

/** @brief Print the default usage summary with all available sections.

	@param[in] argpar argpar object from which to generate the usage summary.
	@param[in] stream Stream to which the help gets printed.
	@param[in] flags Flags controlling which sections to print.  Currently unused.
	@param[in] name Program name to use.

	@return 0 on success, non-0 on error.
*/
int argpar_help(argpar *argpar, FILE *stream, unsigned flags, char *name);

/** @brief Returns the source code version and the implemented API version.

	No assumptions are made about the format of the return value.
*/
const char *argpar_version();

#ifdef __cplusplus
}
#endif

#endif // __ARGPAR_H__
