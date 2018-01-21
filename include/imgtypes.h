/*
	$Id: imgtypes.h,v 1.3 1992/12/09 19:19:20 angel Exp $

	$Log: imgtypes.h,v $
 * Revision 1.3  1992/12/09  19:19:20  angel
 * Add support for Alpha.
 *
 * Revision 1.2  1992/07/13  13:24:17  angel
 * Change from int32 to INT32, etc.  This is the start of support for the
 * Alpha chip.
 *
 * Revision 1.1  1992/06/25  19:53:02  angel
 * Initial release.
 *
*/
#ifndef _IMGTYPES_H
#define _IMGTYPES_H

#ifndef UINT32
#define	UINT32	unsigned int
#else
#define	UINT32	uint32_t
#endif /* UINT32 */

#ifndef INT32
#if defined(__alpha) || (_MIPS_SZLONG == 64)
#define	INT32	int
#else	/* !__alpha */
#define	INT32	int32_t
#endif	/* __alpha */
#endif /* INT32 */

#ifndef UINT16
#define	UINT16	unsigned short
#endif /* INT16 */

#ifndef INT16
#define	INT16	short
#endif /* INT16 */

#ifndef FLOAT32
#define	FLOAT32	float
#endif /* FLOAT32 */

#ifndef FLOAT64
#define	FLOAT64	double
#endif /* FLOAT64 */

#endif /* _IMGTYPES_H */
