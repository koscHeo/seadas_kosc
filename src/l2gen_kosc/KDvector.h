/* v.cpp: implementation of the v class.
//
///////////////////////////////////////////////////////////////////// */

/*/////////////////////////////////////////////////////////////////////
// Construction/Destruction
///////////////////////////////////////////////////////////////////// */

#define vdiff(x1, x2, y1, y2, z1, z2) ((((x1)-(x2))*((x1)-(x2)))+(((y1)-(y2))*((y1)-(y2)))+(((z1)-(z2))*((z1)-(z2))))
#define minreach(x,xstart) (((x)>=(xstart)) ? (x) : (xstart))
#define maxreach(x,xend) (((x)<=(xend)) ? (x) : (xend))
#define vlength(x,y,z) ((x)*(x)+(y)*(y)+(z)*(z))

typedef struct v {

	float m_dvec[3];
	int32_t  m_pos;
	
} vector;


typedef struct KDtree KDtree;

struct KDtree {

   /* m_fnInsert copies the input objects into a binary NEAR tree. When a node has
   // two entries, a descending node is used or created. The current datum is
   // put into the branch descending from the nearer of the two
   // objects in the current node.

   // m_bfnNearestNeighbor retrieves the object nearest to some probe by descending
   // the tree to search out the appropriate object. Speed is gained
   // by pruning the tree if there can be no data below that are
   // nearer than the best so far found.

   // The tree is built in time O(n log n), and retrievals take place in
   // time O(log n).							     */


    vector *   m_ptLeft;            /* left object (of type T) stored in this node  */
    vector *   m_ptRight;           /* right object (of type T) stored in this node */
    float      m_dMaxLeft;          /* longest distance from the left object to     */
                                    /* anything below it in the tree 		    */
    float      m_dMaxRight;         /* longest distance from the right object to    */
                                    /* anything below it in the tree 		    */
    KDtree *   m_pLeftBranch;       /* tree descending from the left object	    */
    KDtree *   m_pRightBranch;      /* tree descending from the right object	    */

};


void v_insert( vector **m_pt, vector *t );
void newKDtree( KDtree **Tree );
void freeKDtree( KDtree *Tree );
void m_fnInsert( KDtree **Tree, vector *t );
int m_bfnNearestNeighbor ( KDtree  * Tree, float  dRadius,  vector  * tClosest, vector  * t );
int m_bfnFarthestNeighbor ( KDtree  * Tree, vector  ** tFarthest, vector  * t );
int32_t m_lfnFindInSphere ( KDtree  * Tree, float dRadius, vector ** tClosest, vector * t );
void m_lfnInSphere ( KDtree  * Tree, float dRadius, vector  * t, vector ** tClosest, int32_t  * lReturn );
int m_bfnNearest ( KDtree  * Tree, float  * dRadius, vector  * tClosest, vector  * t );
int m_bfnFindFarthest ( KDtree  * Tree, float  * dRadius, vector  ** tFarthest, vector  * t );
void alloc_Vector ( vector *v, float da, float db, float dc, int32_t dd );
float vdistance( vector *v1, vector *v2 );






