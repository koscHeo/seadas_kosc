#include <stdio.h>
#include <stdlib.h>
#include "grb2.h"
#include "wgrib2.h"
#include "fnlist.h"

#ifdef USE_REGEX
#include <sys/types.h>
#include <regex.h>

/*
 * match option
 *   requires POSIX-2
 *    if regcomp and regexec are not available,  don't define USE_REGEX in makefile
 *
 * Only process messages that match the "-s inventory" using a regular expression
 *  This short circuits having to read the file twice.
 *  Can call f_match a multiple number of times.
 *  ex. -match (dates)  -match (levels)  -match (variables)
 *
 */

/*
 * HEADER:100:match:setup:1:process data that matches X (regular expression)
 */

extern int match;			/* run matching routines */

static int match_count = 0;
static regex_t preg[MATCH_MAX];
static int type[MATCH_MAX];
static int match_val[MATCH_MAX];

int f_match(ARG1)  {
    if (mode == -1) {
        if (match_count >= MATCH_MAX) fatal_error("too many -match, -not options","");
        match = 1;
        if (regcomp(&(preg[match_count]), arg1, REG_EXTENDED | REG_NOSUB | REG_NEWLINE)) 
            fatal_error("bad regular expression \"%s\"", arg1);
        type[match_count] = 1;
        match_count++;
    }
    return 0;
}

/*
 * HEADER:100:not:setup:1:process data that does not match X (regular expression)
 */
int f_not(ARG1)  {
    if (mode == -1) {
        if (match_count >= MATCH_MAX) fatal_error("too many -match, -not options","");
        match = 1;
        if (regcomp(&(preg[match_count]), arg1, REG_EXTENDED | REG_NOSUB | REG_NEWLINE)) 
            fatal_error("bad regular expression \"%s\"", arg1);
        type[match_count] = 0;
        match_count++;
    }
    return 0;
}

/*
 * HEADER:100:if:misc:1:process data that matches X (regular expression)
 */

int f_if(ARG1)  {
    struct local_struct {
        int match_cnt;
    };
    struct local_struct *save;


    if (mode == -1) {
        if (match_count >= MATCH_MAX) fatal_error("too many -match, -not -if options","");
        match = 1;
        if (regcomp(&(preg[match_count]), arg1, REG_EXTENDED | REG_NOSUB | REG_NEWLINE)) 
            fatal_error("bad regular expression \"%s\"", arg1);
        type[match_count] = 2;

        save = (struct local_struct *) malloc( sizeof(struct local_struct));
        if (save == NULL) fatal_error("memory allocation f_if","");
        save -> match_cnt = match_count;
	*local = save;
        match_count++;
    }
    else if (mode >= 0) {
	save = *local;
	printf("if --  %d = %d \n",save->match_cnt,match_val[save->match_cnt]);
    }
    return 0;
}

/*
 * is_match
 *
 * return codes
 *     1                      ignore
 *     0                      process
 */


int is_match(char *s) {
    int i, j;

    /* process  if-tests */

    for (i = 0; i < match_count; i++) {
        if (type[i] == 2) match_val[i] = regexec(&(preg[i]), s, (size_t) 0, NULL, 0);
    }

    /* process match and not tests */

    for (i = 0; i < match_count; i++) {
        if (type[i] == 2) continue;
	j = regexec(&(preg[i]), s, (size_t) 0, NULL, 0) != 0;
	if (j == type[i]) return 1;
    }
    return 0;
}
#else
int f_match(ARG1)  {
   if (mode == -1) {fprintf(stderr,"MATCH package not installed\n"); return 1;}
   return 1;
}
int f_not(ARG1)  {
   if (mode == -1) {fprintf(stderr,"MATCH package not installed\n"); return 1;}
   return 1;
}
#endif

/*
 * HEADER:100:end:misc:0:stop after first (sub)message (save time)
 */
extern int last_message;		/* last message to process */
int f_end(ARG0) {
    if (mode >= 0) last_message = 1;
    return 0;
}

