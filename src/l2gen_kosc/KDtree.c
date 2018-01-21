#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include "KDvector.h" 


/* ------------------------------------------------------------------ 

// Nearest Neighbor algorithm after Kalantari and McDonald,
// (IEEE Transactions on Software Engineering, v. SE-9, pp.
//    631-634,1983)

//  converted from C++ to C by Ewa Kwiatkowska, NASA GSFC, SAIC
//  November 2003

//  modified to use recursion instead of a double-linked tree
//  and simplified so that it does a bit less checking for
//  things like is the vdistance to the right less than the
//  vdistance to the left; it was found that these checks little
//  to no difference.

// copyright by Larry Andrews, 2001
//  may be freely distributed or used as long as this copyright notice
//  is included

// This template is used to contain a collection of objects. After the
// collection has been loaded into this structure, it can be quickly
// queried for which object is "closest" to some probe object of the
// same type. The major restriction on applicability of the near-tree
// is that the algorithm only works if the objects obey the triangle
// inequality. The triangle rule states that the length of any side of
// a triangle cannot exceed the sum of the lengths of the other two sides.

// The user of this class needs to provide at least the following
// functionality for the template to work. For the built-in
// numerics of C++, they are provided by the system.

//    operator vdistance( );   // conversion constructor from the templated class to double
//                            (usually will return a "length")
//    operator- ( );          // geometrical (vector) difference of two objects
//    a copy constructor
//    a constructor would be nice
//    a destructor would be nice

// The provided interface is:
//
//    #include "TNear.h"
//
//    CNearTree( void )   // constructor
//       instantiated by something like:      CNearTree <v> vTree;
//       for some type v
//
//    void m_fnInsert( T& t )
//       where t is an object of the type v
//
//    bool m_bfnNearestNeighbor ( const double& dRadius,  T& tClosest,   const T& t ) const
//       dRadius is the largest radius within which to search; make it
//          very large if you want to include every point that was loaded; dRadius
//          is returned as the closest vdistance to the probe (or the search radius
//          if nothing is found)
//       tClosest is returned as the object that was found closest to the probe
//          point (if any were within radius dRadius of the probe)
//       t is the probe point, used to search in the group of points m_fnInsert'ed
//       return value is true if some object was found within the search radius, false otherwise
//
//    bool m_bfnFarthestNeighbor ( T& tFarthest,   const T& t ) const
//       tFarthest is returned as the object that was found farthest to the probe
//          point
//       t is the probe point, used to search in the group of points m_fnInsert'ed
//       return value is true if some object was found, false otherwise
//
//    int32_t m_lfnFindInSphere ( const double& dRadius,  std::vector<  T >& tClosest,   const T& t ) const
//       dRadius is the radius within which to search; make it very large if you want to
//           include every point that was loaded;
//       tClosest is returned as the vector of objects that were found within a radius dRadius
//          of the probe point
//       t is the probe point, used to search in the group of points m_fnInsert'ed
//       return value is the number of objects found within the search radius
//
//    ~CNearTree( void )  // destructor
//       invoked by  vTree.CNeartree<v>::~CNearTree
//       for an object vTree of some type v

// So a complete program is:
//
// #include "TNear.h"
// #include <cstdio>
// void main()
// {
//   CNearTree< double > dT;
//   double dNear;
//   dT.m_fnInsert( 1.5 );
//   if ( dT.m_bfnNearestNeighbor( 10000.0,   dNear,  2.0 )) printf( "%f\n",dRad );
// }
//
// and it should print 0.5 (that's how for 2.0 is from 1.5)
//
//
//-------------------------------------------------------------------------  */

/*=======================================================================
//  CNearTree ( )
//
//  Default constructor for class CNearTree
//  creates an empty tree with no right or left node and with the dMax-below
//  set to negative values so that any match found will be stored since it will
//  greater than the negative value
//
//=======================================================================   */
void newKDtree( KDtree **Tree )
   {
   
      if ( (*Tree = (KDtree *) malloc ( sizeof(KDtree) )) == NULL ) {
         printf("Cannot allocate memory to a KDtree\n");
	 exit(1);
      }
      (*Tree)->m_ptLeft       = NULL;
      (*Tree)->m_ptRight      = NULL;
      (*Tree)->m_pLeftBranch  = NULL;
      (*Tree)->m_pRightBranch = NULL;
      (*Tree)->m_dMaxLeft     = FLT_MIN;
      (*Tree)->m_dMaxRight    = FLT_MIN;
      
   }  /*  CNearTree constructor   */

/*=======================================================================
//  ~CNearTree ( )
//
//  Destructor for class CNearTree
//
//=======================================================================   */
void freeKDtree( KDtree *Tree )
   {

      if (Tree->m_pLeftBranch != NULL) freeKDtree(Tree->m_pLeftBranch);
      if (Tree->m_pLeftBranch != NULL) free(Tree->m_pLeftBranch);
      Tree->m_pLeftBranch  = NULL;
      if (Tree->m_pRightBranch != NULL) freeKDtree(Tree->m_pRightBranch);
      if (Tree->m_pRightBranch != NULL) free(Tree->m_pRightBranch);
      Tree->m_pRightBranch = NULL;
      if (Tree->m_ptLeft != NULL)  free(Tree->m_ptLeft);
      Tree->m_ptLeft       = NULL;
      if (Tree->m_ptRight != NULL) free(Tree->m_ptRight);
      Tree->m_ptRight      = NULL;
      Tree->m_dMaxLeft     = FLT_MIN;
      Tree->m_dMaxRight    = FLT_MIN;
          
   }   /*  CNearTree destructor   */

/*=======================================================================
//  void m_fnInsert ( const T& t )
//
//  Function to insert some "point" as an object into a CNearTree for
//  later searching
//
//     t is an object of the templated type which is to be inserted into a
//     Neartree
//
//  Three possibilities exist: put the datum into the left
//  postion (first test),into the right position, or else
//  into a node descending from the nearer of those positions
//  when they are both already used.
//
//=======================================================================     */
   void m_fnInsert( KDtree **Tree, vector *t )
   {
      /* do a bit of precomputing if it is possible so that we can            */
      /* reduce the number of calls to operator 'double' as much as possible; */
      /* 'double' might use square roots in some cases                        */
      float  dTempRight =  0.0;
      float  dTempLeft  =  0.0;

      if ( (*Tree)->m_ptRight != NULL )
      {
         dTempRight  = fabs( vdistance( t, (*Tree)->m_ptRight ) );
         dTempLeft   = fabs( vdistance( t, (*Tree)->m_ptLeft  ) );
      }

      if ( (*Tree)->m_ptLeft == NULL )
      {
         v_insert( &((*Tree)->m_ptLeft), t );
      }
      else if ( (*Tree)->m_ptRight == NULL )
      {
         v_insert( &((*Tree)->m_ptRight), t );
      }
      else if ( dTempLeft > dTempRight )
      {
         if ( (*Tree)->m_pRightBranch == NULL ) newKDtree(&((*Tree)->m_pRightBranch));
	 
         /* note that the next line assumes that m_dMaxRight is negative for a new node */
         if ( (*Tree)->m_dMaxRight < dTempRight ) (*Tree)->m_dMaxRight = dTempRight;
	 
         m_fnInsert( &((*Tree)->m_pRightBranch), t );
      }
      else  /* ((double)(t - *m_tLeft) <= (double)(t - *m_tRight) )     */
      {
         if ( (*Tree)->m_pLeftBranch  == NULL ) newKDtree(&((*Tree)->m_pLeftBranch));
	 
         /* note that the next line assumes that m_dMaxLeft is negative for a new node  */
         if ( (*Tree)->m_dMaxLeft < dTempLeft ) (*Tree)->m_dMaxLeft = dTempLeft;
	 
         m_fnInsert( &((*Tree)->m_pLeftBranch), t );
      }

   }  /*  m_fnInsert     */
   
   void v_insert( vector **m_pt, vector *t )
   {
       if ( (*m_pt = (vector *)malloc(sizeof(vector))) == NULL ) {
         printf("Cannot allocate memory to the vector\n");
	 exit(1);
       }
       memcpy( (void *)(*m_pt), (const void *)t, sizeof(vector) );
   }

/*=======================================================================
//  bool m_bfnNearestNeighbor ( const double& dRadius,  T& tClosest,   const T& t ) const
//
//  Function to search a Neartree for the object closest to some probe point, t. This function
//  is only here so that the function m_bfnNearest can be called without having dRadius const
//
//    dRadius is the maximum search radius - any point farther than dRadius from the probe
//             point will be ignored
//    tClosest is an object of the templated type and is the returned nearest point
//             to the probe point that can be found in the Neartree
//    t  is the probe point
//    the return value is true only if a point was found
//
//=======================================================================    */
   int m_bfnNearestNeighbor ( KDtree  * Tree, float  dRadius,  vector  * tClosest, vector  * t )
   {
      float  dSearchRadius = dRadius;
      
      return ( m_bfnNearest ( Tree, &dSearchRadius, tClosest, t ) );
      
   }  /*  m_bfnNearestNeighbor    */

/*=======================================================================
//  bool m_bfnFarthestNeighbor ( const double& dRadius,  T& tClosest,   const T& t ) const
//
//  Function to search a Neartree for the object closest to some probe point, t. This function
//  is only here so that the function m_bfnFarthestNeighbor can be called without the user
//  having to input a search radius and so the search radius can be guaranteed to be
//  negative at the start.
//
//    tFarthest is an object of the templated type and is the returned farthest point
//             from the probe point that can be found in the Neartree
//    t  is the probe point
//    the return value is true only if a point was found (should only be false for
//             an empty tree)
//
//=======================================================================    */
   int m_bfnFarthestNeighbor ( KDtree  * Tree, vector  ** tFarthest, vector  * t )
   {
     float  dSearchRadius = FLT_MIN;
     
     return ( m_bfnFindFarthest ( Tree, &dSearchRadius, tFarthest, t ) );
     
   }  /*  m_bfnFarthestNeighbor    */

/*=======================================================================
//  int32_t m_lfnFindInSphere ( const double& dRadius,  std::vector<  T >& tClosest,   const T& t ) const
//
//  Function to search a Neartree for the set of objects closer to some probe point, t,
//  than dRadius. This is only here so that tClosest can be cleared before starting the work.
//
//    dRadius is the maximum search radius - any point farther than dRadius from the probe
//             point will be ignored
//    tClosest is a vector of objects of the templated type and is the returned set of nearest points
//             to the probe point that can be found in the Neartree
//    t  is the probe point
//    return value is the number of points found within dRadius of the probe point
//
//=======================================================================    */
   int32_t m_lfnFindInSphere ( KDtree  * Tree, float dRadius, vector  ** tClosest, vector  * t )
   {
   
      int32_t lReturn = 0;
      
      m_lfnInSphere( Tree, dRadius, t, tClosest, &lReturn );
      
      return ( lReturn );
   }  /*  m_lfnFindInSphere    */


/*=======================================================================
//  int32_t m_lfnInSphere ( const double& dRadius,  std::vector<  T >& tClosest,   const T& t ) const
//
//  Private function to search a Neartree for the object closest to some probe point, t.
//  This function is only called by m_lfnFindInSphere.
//
//    dRadius is the search radius
//    tClosest is a vector of objects of the templated type found within dRadius of the
//         probe point
//    t  is the probe point
//    the return value is the number of points found within dRadius of the probe
//
//=======================================================================    */
   void m_lfnInSphere ( KDtree  * Tree, float dRadius, vector  * t, vector ** tClosest, int32_t  * lReturn )
   {
   
      /* first test each of the left and right positions to see if        */
      /* one holds a point nearer than the search radius.    	          */
      
      if ( *lReturn >= 10000 ) return;
      
      if (( Tree->m_ptLeft != NULL ) && (( fabs( vdistance( t, Tree->m_ptLeft  ))) <= dRadius ))
      {
/*	 v_insert( &(tClosest[*lReturn++]), Tree->m_ptLeft );   */
	 memcpy( (void *)tClosest[*lReturn], (const void *)Tree->m_ptLeft, sizeof(vector) );
	 ++(*lReturn);
      }
      if (( Tree->m_ptRight != NULL ) && (( fabs( vdistance( t, Tree->m_ptRight ))) <= dRadius ))
      {
/*         v_insert( &(tClosest[*lReturn++]), Tree->m_ptRight );   */
	 memcpy( (void *)tClosest[*lReturn], (const void *)Tree->m_ptRight, sizeof(vector) );
	 ++(*lReturn);
      }
      /*   								  */
      /* Now we test to see if the branches below might hold an object    */
      /* nearer than the search radius. The triangle rule is used         */
      /* to test whether it's even necessary to descend.                  */
      /*    								  */
      if (( Tree->m_pLeftBranch  != NULL )  && (( dRadius + Tree->m_dMaxLeft  ) >= fabs( vdistance( t, Tree->m_ptLeft  ))))
      {
         m_lfnInSphere( Tree->m_pLeftBranch, dRadius, t, tClosest, lReturn );
      }

      if (( Tree->m_pRightBranch != NULL )  && (( dRadius + Tree->m_dMaxRight ) >= fabs( vdistance( t, Tree->m_ptRight ))))
      {
         m_lfnInSphere( Tree->m_pRightBranch, dRadius, t, tClosest, lReturn );
      }

   }  /*  m_lfnInSphere    */

/*=======================================================================
//  bool m_bfnNearest ( double& dRadius,  T& tClosest,   const T& t ) const
//
//  Private function to search a Neartree for the object closest to some probe point, t.
//  This function is only called by m_bfnNearestNeighbor.
//
//    dRadius is the smallest currently known vdistance of an object from the probe point.
//    tClosest is an object of the templated type and is the returned closest point
//             to the probe point that can be found in the Neartree
//    t  is the probe point
//    the return value is true only if a point was found within dRadius
//
//=======================================================================    */
   int m_bfnNearest ( KDtree  * Tree, float  * dRadius, vector  * tClosest, vector  * t )
   {
      float    dTempRadius;
      int      bRet = 0;

      /* first test each of the left and right positions to see if       */
      /* one holds a point nearer than the nearest so far discovered.    */
      if (( Tree->m_ptLeft != NULL ) && (( dTempRadius = fabs( vdistance( t, Tree->m_ptLeft ))) <= *dRadius ))
      {
         *dRadius  = dTempRadius;
         memcpy( (void *)(tClosest), (const void *)Tree->m_ptLeft, sizeof(vector) );
         bRet     = 1;
      }
      if (( Tree->m_ptRight != NULL ) && (( dTempRadius = fabs( vdistance( t, Tree->m_ptRight ))) <= *dRadius ))
      {
         *dRadius  = dTempRadius;
         memcpy( (void *)(tClosest), (const void *)Tree->m_ptRight, sizeof(vector) );
         bRet     = 1;
      }

      /*  								  */
      /* Now we test to see if the branches below might hold an object    */
      /* nearer than the best so far found. The triangle rule is used     */
      /* to test whether it's even necessary to descend.                  */
      /* 								  */
      if (( Tree->m_pLeftBranch  != NULL )  && (( *dRadius + Tree->m_dMaxLeft  ) >= fabs( vdistance( t, Tree->m_ptLeft ))))
      {
         bRet |= m_bfnNearest( Tree->m_pLeftBranch, dRadius, tClosest, t );
      }

      if (( Tree->m_pRightBranch != NULL )  && (( *dRadius + Tree->m_dMaxRight ) >= fabs( vdistance( t, Tree->m_ptRight ))))
      {
         bRet |= m_bfnNearest( Tree->m_pRightBranch, dRadius, tClosest, t );
      }

      return ( bRet );
   }   /* m_bfnNearest    */

/*=======================================================================
//  bool m_bfnFindFarthest ( double& dRadius,  T& tFarthest,   const T& t ) const
//
//  Private function to search a Neartree for the object farthest from some probe point, t.
//  This function is only called by m_bfnFarthestNeighbor.
//
//    dRadius is the largest currently known vdistance of an object from the probe point.
//    tFarthest is an object of the templated type and is the returned farthest point
//             from the probe point that can be found in the Neartree
//    t  is the probe point
//    the return value is true only if a point was found (should only be false for
//             an empty tree)
//
//=======================================================================    */
   int m_bfnFindFarthest ( KDtree  * Tree, float  * dRadius, vector  ** tFarthest, vector  * t )
   {
      float    dTempRadius;
      int      bRet    = 0;

      /* first test each of the left and right positions to see if              */
      /* one holds a point farther than the farthest so far discovered.         */
      /* the calling function is presumed initially to have set dRadius to a    */
      /* negative value before the recursive calls to m_bfnFindFarthestNeighbor */

      if (( Tree->m_ptLeft != NULL  ) && (( dTempRadius = fabs( vdistance( t, Tree->m_ptLeft ))) >= *dRadius ))
      {
         *dRadius   = dTempRadius;
         memcpy( (void *)(*tFarthest), (const void *)Tree->m_ptLeft, sizeof(vector) );
         bRet      = 1;
      }
      if (( Tree->m_ptRight != NULL ) && (( dTempRadius = fabs( vdistance( t, Tree->m_ptRight))) >= *dRadius ))
      {
         *dRadius   = dTempRadius;
         memcpy( (void *)(*tFarthest), (const void *)Tree->m_ptRight, sizeof(vector) );
         bRet      = 1;
      }
      /*                                                                         */
      /* Now we test to see if the branches below might hold an object           */
      /* farther than the best so far found. The triangle rule is used           */
      /* to test whether it's even necessary to descend.                         */
      /*                                                                         */
      if (( Tree->m_pLeftBranch  != NULL )  && (( *dRadius - Tree->m_dMaxLeft  ) <= fabs( vdistance( t, Tree->m_ptLeft  ))))
      {
         bRet |=  m_bfnFindFarthest( Tree->m_pLeftBranch, dRadius, tFarthest, t );
      }

      if (( Tree->m_pRightBranch != NULL )  && (( *dRadius - Tree->m_dMaxRight ) <= fabs( vdistance( t, Tree->m_ptRight ))))
      {
         bRet |=  m_bfnFindFarthest( Tree->m_pRightBranch, dRadius, tFarthest, t );
      }

      return ( bRet );
   }   /* m_bfnFindFarthest     */



void newVector ( vector  * v )
{
	v->m_dvec[0]		= FLT_MAX;
	v->m_dvec[1]		= FLT_MAX;
	v->m_dvec[2]		= FLT_MAX;
	v->m_pos		= -1L;
	
}



void alloc_Vector ( vector *v, float da, float db, float dc, int32_t dd )
{
	v->m_dvec[0]		= da;
	v->m_dvec[1]		= db;
	v->m_dvec[2]		= dc;
	v->m_pos		= dd;
}





float vdistance( vector *v1, vector *v2 )
{
   /* return the Euclidean (L2) length of the vector   */

   return ( sqrt(vlength(v1->m_dvec[0]-v2->m_dvec[0],v1->m_dvec[1]-v2->m_dvec[1],v1->m_dvec[2]-v2->m_dvec[2])) );
   
   /* city block measure (L1) 	return ( fabs(m_dvec[0]) + fabs(m_dvec[1]) + fabs(m_dvec[2]) );
   // max value (L-infinity) return the largest of fabs of any element
   // return ( fabs(m_dvec[0])>=fabs(m_dvec[1]) && fabs(m_dvec[0])>=fabs(m_dvec[2]) ? fabs(m_dvec[0]) : fabs(m_dvec[1])>=fabs(m_dvec[2]) ? fabs(m_dvec[1]) : fabs(m_dvec[2]) );
   // Hamming measure (if this is a difference vector) return ( (m_dvec[0]==0 ? 0 : 1) + (m_dvec[1]==0 ? 0 : 1) + (m_dvec[2]==0 ? 0 : 1) )   */
}

