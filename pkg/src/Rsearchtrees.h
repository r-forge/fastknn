#include <stdio.h>
#include <R.h>
#include <Rdefines.h>


typedef struct qtree {
  double upper;
  double lower;
  double left;
  double right;
  int numdata;
  char pos; //to indicate which sector in the parent node this node is. 1 = uleft 2 = uright 3=lright 4 = lleft;
  unsigned char depth;
  int *data;
  char datatype;//indicates the type of objects stored in the tree. 1=points 2=rectangle
  struct qtree *parent;
  struct qtree *uleft;
  struct qtree *uright;
  struct qtree *lleft;
  struct qtree *lright;
} qtree_t;


typedef struct point {
  double x;
  double y;
  int index;
} point_t;

typedef struct rect {
  //parallel to the axes only, other rectangles will need to wait until we do the general polygon case;
  double left; //lower x bound
  double right; //upper x bound
  double low; //lower y bound
  double high; //lower x bound
  struct rect *parent;
  int index;
} rect_t;


typedef struct kdtree {
  int depth;
  double split;
  int data;
  struct kdtree *parent;
  struct kdtree *left;
  struct kdtree *right;
} kdtree_t;

double eucl_Dist(double *pt1, double *pt2, int d);
// function prototypes
qtree_t *Build_Branch(double up, double low, double left, double right, qtree_t *parent, char pos, unsigned char depth);
void Free_Quad_Tree( qtree_t *tree);
qtree_t *Descend_To_Container(qtree_t *node, double newx, double newy);
int Add_Bucket_Point(qtree_t *node, double *x, double *y, int *newdat, int newdatsize, unsigned char maxDepth);
void Add_Depth(qtree_t *curnode);
void Insert_Dist(double *dists, double newdist, int *inds, int newind, int k, int start);
void Harvest_Data_KNN(qtree_t *node, int excludepos, double leftbound, double rightbound, double lowbound, double highbound, int *inds, double *dists, double newx, double newy, double *x, double *y, int k, int start);
int Check_Bounds(qtree_t *node, double left, double right, double down, double up);
void Check_Children(qtree_t *tree, double left, double right, double down, double up, int *ret);
void Get_Data_In_Box(qtree_t *tree, double left, double right, double down, double up, int **found, int *pos, int *cursize, double *x, double *y);
void Get_Data_For_Test( qtree_t *node, int *pos, SEXP ans);
void Free_Quad_Tree( qtree_t *tree);
int Find_Max_Depth(qtree_t *tree, unsigned char curdepth);
int * Grow_Return_Array(int **ar, int *size);
void Get_XY_Ptrs( SEXP dattype, SEXP cols, SEXP Rx, SEXP Ry, double **ptrs);
void Get_Tree_Attributes(qtree_t *tree, int *curattr);

SEXP R_Build_Quadtree(SEXP Rx, SEXP Ry, SEXP RxMax, SEXP RxMin, SEXP RyMax, SEXP RyMin, SEXP RmaxDepth);
SEXP R_Find_KNN(SEXP Rtree, SEXP Rnewdat, SEXP Rfulldat, SEXP Rnewcols, SEXP Rfullcols, SEXP Rnewtyp, SEXP Rfulltyp, SEXP Rk, SEXP Rnewn, SEXP Rfulln);

SEXP R_Get_Points_In_Box(SEXP Rtree, SEXP Rleft, SEXP Rright, SEXP Rdown, SEXP Rup, SEXP Rx, SEXP Ry, SEXP dattype, SEXP cols);
SEXP R_Get_Top_Left_Order(SEXP Rtree, SEXP Rlen);
void R_Free_Quad_Tree( SEXP treeptr);
SEXP R_Find_Max_Depth(SEXP Rtree);

//kdtree


//typedef struct qtree qtree_t;
  
//Samet spatial data structures
