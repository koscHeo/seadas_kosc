#include <stdio.h>
#include <stdint.h>

/*
 *  (C) Copyright by Ariel Faigon, 1996
 *  Released under the GNU GPL (General Public License) version 2
 *  or any later version (http://www.gnu.org/licenses/licenses.html)
 */
/*---------------            sort.h              --------------*/
/*--------------- sort library customizable file --------------*/

/*
 * This is the key TYPE.
 * Replace this typedef by YOUR key type.
 * e.g. if you're sorting an array of pointers to strings
 * You should do:
 *
 *	typedef char * KEY_T
 *
 * The keys are the items in the array that you're moving
 * around using the SWAP macro.
 * Note: the comparison function may compare any "function"
 * of this key, it doesn't necessarily need to compare the
 * key itself. example: you compare the strings pointed to
 * by the key itself.
 */ 
//typedef int32_t KEY_T;
//typedef uint64_t	KEY_T;
typedef uint32_t	KEY_T;

/*
 * These are the COMPARISON macros
 * Replace these macros by YOUR comparison operations.
 * e.g. if you are sorting an array of pointers to strings
 * you should define:
 *
 *	GT(x, y)  as   (strcmp((x),(y)) > 0) 	Greater than
 *	LT(x, y)  as   (strcmp((x),(y)) < 0) 	Less than
 *	GE(x, y)  as   (strcmp((x),(y)) >= 0) 	Greater or equal
 *	LE(x, y)  as   (strcmp((x),(y)) <= 0) 	Less or equal
 *	EQ(x, y)  as   (strcmp((x),(y)) == 0) 	Equal
 *	NE(x, y)  as   (strcmp((x),(y)) != 0) 	Not Equal
 */
#define GT(x, y) ((x) > (y))
#define LT(x, y) ((x) < (y))
#define GE(x, y) ((x) >= (y))
#define LE(x, y) ((x) <= (y))
#define EQ(x, y) ((x) == (y))
#define NE(x, y) ((x) != (y))

/*
 * This is the SWAP macro to swap between two keys.
 * Replace these macros by YOUR swap macro.
 * e.g. if you are sorting an array of pointers to strings
 * You can define it as:
 *
 *	#define SWAP(x, y) temp = (x); (x) = (y); (y) = temp
 *
 * Bug: 'insort()' doesn't use the SWAP macro.
 */
#define SWAP(x, y) temp = (x); (x) = (y); (y) = temp

/*-------------------- End of customizable part -----------------------*/
/*-------------------- DON'T TOUCH BEYOND THIS POINT ------------------*/
extern void    insort ();
extern void    quicksort ();
extern void    quickersort ();
extern void    partial_quickersort ();
extern void    sedgesort ();
extern void    shellsort ();
extern void    heapsort ();
extern void    gamasort ();
extern void    siftdown ();
extern int     sorted ();
extern int     array_eq ();

//#include  "sort.h"

/*
 *  A library of sorting functions
 *
 *  Written by:  Ariel Faigon,  1987
 *
 *  (C) Copyright by Ariel Faigon, 1996
 *  Released under the GNU GPL (General Public License) version 2
 *  or any later version (http://www.gnu.org/licenses/licenses.html)
 */

/* http://www.yendor.com/programming/sort */

/*
 |  void  insort (array, idx, len)
 |  KEY_T  array[];
 |  int    len;
 |
 |  Abstract:	Sort array[0..len-1] into increasing order.
 |
 |  Method:	Optimized insertion-sort (ala Jon Bentley)
 */

void  insort (array, idx, len)
register KEY_T  array[];
register uint32_t   idx[];
register int    len;
{
	register int	i, j;
	register KEY_T	temp;
	register int32_t temp_idx;

	for (i = 1; i < len; i++) {
		/* invariant:  array[0..i-1] is sorted */
		j = i;
		/* customization bug: SWAP is not used here */
		temp = array[j];
		temp_idx = idx[j];
		while (j > 0 && GT(array[j-1], temp)) {
			array[j] = array[j-1];
			idx[j] = idx[j-1];
			j--;
		}
		array[j] = temp;
		idx[j] = temp_idx;
	}
}


/* 15 has been found empirically as the optimal cutoff value in 1996 */
/* In 2006, with computers very different and much faster it was found
 * to be a close tie between 15 and 16 */
#ifndef CUTOFF
#  define CUTOFF 15
#endif

/*
 |  void  partial_quickersort (array, lower, upper)
 |  KEY_T  array[];
 |  int    lower, upper;
 |
 |  Abstract:
 |	Sort array[lower..upper] into a partial order
 |     	leaving segments which are CUTOFF elements long
 |     	unsorted internally.
 |
 |  Efficiency:
 |	Could be made faster for _worst_ cases by selecting
 |	a pivot using median-of-3. I don't do it because
 |	in practical cases my pivot selection is arbitrary and
 |	thus pretty random, your mileage may vary.
 |
 |  Method:
 |	Partial Quicker-sort using a sentinel (ala Robert Sedgewick)
 |
 |  BIG NOTE:
 |	Precondition: array[upper+1] holds the maximum possible key.
 |	with a cutoff value of CUTOFF.
 */

void  partial_quickersort (array, idx, lower, upper)
register KEY_T  array[];
register uint32_t   idx[];
register int    lower, upper;
{
    register int	i, j;
    register KEY_T	temp, pivot;

    if (upper - lower > CUTOFF) {
	SWAP(array[lower], array[(upper+lower)/2]);
	SWAP(idx[lower], idx[(upper+lower)/2]);
	i = lower;  j = upper + 1;  pivot = array[lower];
	while (1) {
	    /*
	     * ------------------------- NOTE --------------------------
	     * ignoring BIG NOTE above may lead to an infinite loop here
	     * ---------------------------------------------------------
	     */
	    do i++; while (LT(array[i], pivot));
	    do j--; while (GT(array[j], pivot));
	    if (j < i) break;
	    SWAP(array[i], array[j]);
	    SWAP(idx[i], idx[j]);
	}
	SWAP(array[lower], array[j]);
	SWAP(idx[lower], idx[j]);
	partial_quickersort (array, idx, lower, j - 1);
	partial_quickersort (array, idx, i, upper);
    }
}


/*
 |  void  sedgesort (array, idx, len)
 |  KEY_T  array[];
 |  int    len;
 |
 |  Abstract:
 |	Sort array[0..len-1] into increasing order.
 |
 |  Method:
 |	Use partial_quickersort() with a sentinel (ala Sedgewick)
 |	to reach a partial order, leave the unsorted segments of
 |	length == CUTOFF to a simpler low-overhead, insertion sort.
 |
 |	This method seems to me the ultimative sort method in terms
 |	of average efficiency (Skeptic ? try to beat it).
 |
 |  BIG NOTE:
 |	precondition: array[len] must hold a sentinel (largest
 |	possible value) in order for this to work correctly.
 |	An easy way to do this is to declare an array that has
 | 	len+1 elements [0..len], and assign MAXINT or some such
 |	to the last location before starting the sort (see sorttest.c)
 */
void  sedgesort (array, idx, len)
register KEY_T  array[];
register uint32_t   idx[];
register int    len;
{
    /*
     * ------------------------- NOTE --------------------------
     * ignoring BIG NOTE above may lead to an infinite loop here
     * ---------------------------------------------------------
     */
    partial_quickersort (array, idx, 0, len - 1);
    insort (array, idx, len);
}

