#include "Rsearchtrees.h"
#include <time.h>

static long alloctimetaken = 0;
static long getDatatime = 0;
static long checkBoundtime = 0;
static long checkDatatime = 0;
static long descendtime = 0;
static int clockerrs = 0;
static int childrenChecked = 0;
static int dataChecked = 0;
static int dataInBoxCount = 0;
static struct timespec outside, inside;
static long recursiontime = 0;
SEXP 
R_Build_Quadtree(SEXP Rx, SEXP Ry, SEXP RxMax, SEXP RxMin, SEXP RyMax, SEXP RyMin, SEXP RmaxDepth)
{
  //fprintf(stderr, "Starting R_Build_Quadtree\n"); fflush(stderr);
  double *x = REAL(Rx);
  double *y = REAL(Ry);
  int len = LENGTH(Rx);
  unsigned char maxDepth = (unsigned char) INTEGER( RmaxDepth)[ 0 ];
  
  //qtree_t tree;
  qtree_t *tree;
  double upper = REAL(RyMax)[0];
  double lower = REAL(RyMin)[0];
  double left = REAL(RxMin)[0];
  double right = REAL(RxMax)[0];
  tree = Build_Branch(upper, lower, left, right, NULL, 0, 0 );
  
  int res;
  for (int i=0; i<len; i++)
    {
      //res = Add_Point(&tree, x, y, i);
      Add_Bucket_Point(tree, x, y, &i, 1, maxDepth);
    }

  //get info about the tree;
  
  //first int is number of nodes, second int is number of nodes with data, third int is maxDepth, fourth int is max bucketsize

  int *attr = calloc(4, sizeof(int));
  Get_Tree_Attributes(tree, attr);
  
  SEXP klass, ans, ptr;
  PROTECT( klass = MAKE_CLASS( "QuadTree" ) );
  PROTECT( ans = NEW( klass ) );
  PROTECT( ptr = R_MakeExternalPtr( tree ,
				    Rf_install( "QuadTree" ),
				    R_NilValue ) );
  R_RegisterCFinalizerEx(ptr, &R_Free_Quad_Tree, 1);
  SET_SLOT( ans, Rf_install( "ref" ), ptr );
  SET_SLOT( ans, Rf_install( "points" ), ScalarInteger( len ) );
  SET_SLOT( ans, Rf_install( "numNodes" ), ScalarInteger(attr[0] ));
  SET_SLOT( ans, Rf_install( "dataNodes" ), ScalarInteger(attr[1] ) );
  SET_SLOT( ans, Rf_install( "maxDepth" ), ScalarInteger(attr[2] ) );

  SET_SLOT( ans, Rf_install( "maxBucketSize" ), ScalarInteger( attr[3] ) );
  UNPROTECT(3);
  free(attr);
  return ans;  
}

int Add_Bucket_Point(qtree_t *node, double *x, double *y, int *newdat, int newdatsize, unsigned char maxDepth)
{
  int toret = 0;
  int *olddat, olddatsize;
  qtree_t *curnode = node;
  double newx, newy;
  for ( int i = 0; i < newdatsize; i++)
    {
      newx = x[ newdat[ i ] ];
      newy = y[ newdat[ i ] ];
      //descend to leaf node whose area contains the new cordinates
      curnode = Descend_To_Container(node, newx, newy);

      //if the leaf is empty, insert the data
      if (curnode -> numdata == 0)
	{
	  olddat = calloc(1, sizeof(int));
	  olddat[0] = newdat[i];
	  curnode -> data = olddat;
	  curnode -> numdata ++;
	  toret = 1;
	}
      //else if (curnode -> depth == maxDepth || (x[ *curnode -> data  ] == newx && y[ *curnode -> data ] == newy ) ) 
      else if (curnode -> depth == maxDepth)
      {
	  
      //realloc preserves the values of the previous elements
	  curnode -> data = realloc(curnode -> data, sizeof(int) * (curnode -> numdata + 1));
	  // for (int j = 0; j < newdatsize; j++)
	  //  curnode -> data[curnode -> numdata + j] = newdat[j];
	  curnode -> data [curnode -> numdata ] = newdat[i];
	  
	  curnode -> numdata = curnode -> numdata + 1;
	  //curnode -> data = tmpdat;
	  toret = 1;
	} else {
	//leaf had data and we are not yet at max depth. Thus we add depth and try again
	olddat = curnode -> data;
	olddatsize = curnode -> numdata;
	Add_Depth(curnode);
	curnode -> numdata = 0;
	curnode -> data = NULL;
	Add_Bucket_Point(curnode, x, y, olddat, olddatsize, maxDepth);
	toret = Add_Bucket_Point(curnode, x, y, newdat, newdatsize, maxDepth);
      }
    }
    return toret;
}

  
void Add_Depth(qtree_t *curnode)
{
  
  double up = curnode -> upper;
  double low = curnode -> lower;
  double left = curnode -> left;
  double right = curnode -> right;
  double midvert = (up + low) / 2.0;
  double midhoriz = (left + right) / 2.0;
  unsigned char depth = curnode -> depth + 1;
  curnode -> uleft = Build_Branch(up, midvert, left, midhoriz, curnode, 1, depth);
  curnode -> uright  = Build_Branch(up, midvert, midhoriz, right, curnode, 2, depth);
  curnode -> lleft =  Build_Branch(midvert, low, left, midhoriz, curnode, 4, depth) ;
  curnode -> lright =  Build_Branch(midvert, low, midhoriz, right, curnode, 3, depth);

  return;
}
  
qtree_t *Build_Branch(double up, double low, double left, double right, qtree_t *parent, char pos, unsigned char depth)
{
  qtree_t *toret = malloc(sizeof(qtree_t));
  toret->upper = up;
  toret->left = left;
  toret->right = right;
  toret->lower = low;
  toret->parent = parent;
  //malloc doesn't initialize the allocated memory to 0s.
  toret->uleft = NULL;
  toret->uright = NULL;
  toret->lleft = NULL;
  toret->lright = NULL;
  toret->numdata = 0;
  toret->data = NULL;
  toret->pos = pos;
  toret->depth = depth;
  return toret;
}

static checked = 0;
static int checkpassed = 0;
SEXP
//returns candidate indices for k nearest neighbors
R_Find_KNN(SEXP Rtree, SEXP Rnewdat, SEXP Rfulldat, SEXP Rnewcols, SEXP Rfullcols, SEXP Rnewtyp, SEXP Rfulltyp, SEXP Rk, SEXP Rnewn, SEXP Rfulln)
{
  checkpassed = 0;
  qtree_t *tree = (qtree_t *) R_ExternalPtrAddr( GET_SLOT( Rtree, Rf_install( "ref" ) ) );

  double *x, *y, *newx, *newy;
  
  double **ptrs = calloc(2, sizeof(double*));
  Get_XY_Ptrs(Rnewtyp, Rnewcols, Rnewdat, Rnewn, ptrs );
  newx = ptrs[0];
  newy = ptrs[1];
  ptrs = calloc(2, sizeof(double*));
  Get_XY_Ptrs(Rfulltyp, Rfullcols, Rfulldat, Rfulln, ptrs );
  x = ptrs[0];
  y = ptrs[1];
  
  int n = INTEGER(Rfulln)[0];
  int newn = INTEGER(Rnewn)[0];
  
  qtree_t *curnode;
  int k = INTEGER( Rk ) [ 0 ];
  double dists[ newn * k ];
  int chosen[ newn * k ];
    int ind = 0;
  double newpt[2] = { 0 , 0 };
  double oldpt[2] = { 0 , 0 };
  double tmpdist; 
  int filled;
  int j, pos, tmp, cnt = 0;
  double tmpnewx, tmpnewy;
  for(int i= 0; i < newn*k; i++)
    {
      dists[i] = -1.0;
      chosen[i] = -1;
    }
  for(int l = 0; l < newn; l++)
    {
      cnt = 0;
      tmpnewx = newx[ l ];
      tmpnewy = newy[ l ];
      newpt[0] = tmpnewx; newpt[1] = tmpnewy;

      //descend to leaf closest to our x, y coordinates
      //XXX but sometimes that leaf doesn't have data!!
      curnode = Descend_To_Container(tree, tmpnewx, tmpnewy);
      if( curnode -> numdata > 0)
	{
	  
	  for (int i =0; i < curnode -> numdata ; i ++)
	    {
	      ind = curnode -> data[ i ];
	      oldpt[0] = x[ ind ]; oldpt[1] = y[ ind ];
	      
	      tmpdist = eucl_Dist(newpt, oldpt, 2);
	      //chosen[ l * k  + i ] = ind;
	      //dists[ l * k + i ] = tmpdist;
	      Insert_Dist(&dists, tmpdist, &chosen, ind, k, l*k);    
	      cnt = cnt + 1;
	    }    
	}
      //initialize the rest with bad guesses just to have numbers to compare to.
      
      
      if (ind > n - k)
	ind = 0;
      for (int j = ind + 1; j <= ind + 1 + (k - cnt) ; j++)
	{
	  oldpt[0] = x[j]; oldpt[1] = y[j];
	  tmpdist = eucl_Dist(newpt, oldpt, 2);
	  
	  Insert_Dist(&dists, tmpdist, &chosen, j, k, l*k);
	}
      
      //maxdist, ie distance to beat to get into the list, is in dists[ (l + 1) *  k - 1 ]
      tmp = ( l + 1 ) * k - 1;
      while(curnode -> parent != NULL)
	{
	  pos = curnode -> pos;
	  curnode = curnode -> parent;
	  Harvest_Data_KNN(curnode, pos, tmpnewx - dists[ tmp ], tmpnewx + dists[ tmp ], tmpnewy - dists[ tmp ], tmpnewy + dists[ tmp ], &chosen, &dists, tmpnewx, tmpnewy, x, y, k, l * k);
	}
    }
  SEXP ans;
  PROTECT( ans = NEW_INTEGER( newn * k ) );
  for (int i = 0; i < newn * k; i++)
    INTEGER( ans )[ i ] = chosen[ i ] + 1;
		    //INTEGER(ans) = & chosen;
  UNPROTECT(1);
  return ans;
}

qtree_t *Descend_To_Container(qtree_t *node, double newx, double newy)
{
  qtree_t *curnode = node;
  while(curnode -> uleft != NULL)
    {
      if (newx < ( curnode -> lleft ) -> right)
	{
	  if ( newy < (curnode -> lleft ) -> upper)
	    curnode = curnode -> lleft;
	  else
	    curnode = curnode -> uleft;
	} else 
	{
	  if ( newy < (curnode -> lleft ) -> upper)
	    curnode = curnode -> lright;
	  else
	    curnode = curnode -> uright;
	}
    }
  return curnode;

}

void Insert_Dist(double *dists, double newdist, int *inds, int newind, int k, int start)
{

  int pos = start;
  int done = 0;
  while (!done && pos < start + k)
    {
      //if the index is already in the list, we're done;
      if (inds[pos] == newind)
	return;
      //tie goes to the new index.
      if (dists[ pos ] >= newdist || dists[ pos ] == -1.0)
	done = 1;
      else
	pos = pos + 1;
    }
  
  //if we found a place for it;
 
  if (done) 
    {
      int tmpind;
      double tmpdist;
      for(int curpos = pos; curpos < start +  k; curpos ++)
	{
	  
	  tmpind = inds[ curpos ]; //old value at curpos
	  inds[ curpos ] = newind; //newvalue at curpos
	  tmpdist = dists[ curpos ];
	  dists[ curpos ] = newdist;
	  if (curpos + 1 < start + k)
	    {
	      newind = tmpind;
	      newdist = tmpdist;
	    }    
	}
    }
  return;
}  

void
Harvest_Data_KNN(qtree_t *node, int excludepos, double leftbound, double rightbound, double lowbound, double highbound, int *inds, double *dists, double newx, double newy, double *x, double *y, int k, int start)
{
  int datind;
  if (node -> numdata > 0)
    {
      
      double newpt[2] = {newx, newy};
      double oldpt[2];
      double dist;
      for (int i =0; i < node -> numdata; i++)
	{
	  datind = node -> data[ i ]; 
	  oldpt[0] = x[ datind ]; 
	  oldpt[1] = y[ datind ];
	  dist = eucl_Dist(newpt, oldpt, 2 );
	  Insert_Dist(dists, dist, inds, datind, k, start);
	}
      checked += node -> numdata;
    } else {
    
    double midhoriz = ( node -> left + node -> right ) / 2.0;
    double midvert = (node -> upper + node -> lower ) / 2.0;
    //each valid node has 0 or 4 children.
    if (node -> uleft != NULL)
      {
	if (excludepos != 1)
	  {
	    if ( Check_Bounds(node -> uleft, leftbound, rightbound, lowbound, highbound) )
	      Harvest_Data_KNN(node -> uleft, 0, leftbound, rightbound, lowbound, highbound, inds, dists, newx, newy, x, y, k, start);
	  }
	if (excludepos != 2)
	  {
	    if ( Check_Bounds(node -> uright, leftbound, rightbound, lowbound, highbound) )
	      Harvest_Data_KNN(node -> uright, 0, leftbound, rightbound, lowbound, highbound, inds, dists, newx, newy, x, y, k, start);
	  }
	if (excludepos != 3)
	  {
	    if ( Check_Bounds(node -> lright, leftbound, rightbound, lowbound, highbound) )
	      Harvest_Data_KNN(node -> lright, 0, leftbound, rightbound, lowbound, highbound, inds, dists, newx, newy, x, y, k, start);
	  }
	if (excludepos != 4)
	  {
	    if ( Check_Bounds(node -> lleft, leftbound, rightbound, lowbound, highbound) )
	      Harvest_Data_KNN(node -> lleft, 0, leftbound, rightbound, lowbound, highbound, inds, dists, newx, newy, x, y, k, start);
	  }
      }
  }
  return;
}

int 
Check_Bounds(qtree_t *node, double left, double right, double down, double up)
{
  int res1, res2;
  struct timespec before, after;
  //res1 = clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &before);
  //clockerrs += res1 ? 1 : 0;
  int toret = 0;
  if ( ! (node -> left > right || node -> right < left || node -> upper < down || node -> lower > up) )
    {
      checkpassed++;
      toret = 1;
    }
  //clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &after);
  //clockerrs += res2 ? 1 : 0;
  //checkBoundtime += after.tv_nsec - before.tv_nsec;
  return toret;
}

SEXP
R_Get_Points_In_Box(SEXP Rtree, SEXP Rleft, SEXP Rright, SEXP Rdown, SEXP Rup, SEXP Rx, SEXP Ry, SEXP dattype, SEXP cols)
{
  struct timespec before, after;
  alloctimetaken = 0;
  dataChecked = 0;
  checkDatatime = 0;
  getDatatime = 0;
  clockerrs = 0;
  childrenChecked = 0;
  dataInBoxCount = 0;
  
  double up, down, left, right;
  double *x, *y;
  
  double **ptrs = calloc(2, sizeof(double*));
  Get_XY_Ptrs(dattype, cols, Rx, Ry, ptrs );
  x = ptrs[0];
  y = ptrs[1];
  /*  
  if (INTEGER(dattype)[0] == 1)
    {
      int len = LENGTH(Rx)/INTEGER(Ry)[0];
      int xcol = INTEGER(cols)[ 0 ];
      int ycol = INTEGER(cols)[ 1 ];
	
      x = &(REAL( Rx )[ len * ( xcol - 1 ) ] );
      y = &(REAL( Rx )[ len * ( ycol - 1 ) ] );
    } else if (INTEGER(dattype)[0] == 2) 
    {
      int xcol = INTEGER(cols)[ 0 ];
      int ycol = INTEGER(cols)[ 1 ];
      
      x = REAL(VECTOR_ELT(Rx, xcol - 1) );
      y = REAL(VECTOR_ELT(Rx, ycol - 1 ) ) ;
    
  } else {
    x = REAL(Rx);
    y = REAL(Ry);
  }
  */
  qtree_t *tree = (qtree_t *) R_ExternalPtrAddr( GET_SLOT( Rtree, Rf_install( "ref" ) ) );
  up = REAL(Rup)[0];
  down = REAL(Rdown)[0];
  left = REAL(Rleft)[0];
  right = REAL(Rright)[0];
  int pos=0;
  int size = 100;
  //XXX how many spots should we allocate? currently it is as many as there are points but this is probably way too many!
  int *found = calloc(size, sizeof(int)) ; //start at 100, double when necessary
  //int *found = calloc( INTEGER( GET_SLOT( Rtree, Rf_install( "points" ) ) )[ 0 ], sizeof(int));
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &before);
  Get_Data_In_Box(tree, left, right, down, up, &found, &pos, &size, x, y );
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &after);
  getDatatime += (after.tv_sec - before.tv_sec) * 1000000000 + after.tv_nsec - before.tv_nsec;
  SEXP ans;
  fprintf(stderr, "Total time taken in Grow_Return_Array (nanosecs): %ld \nTotal time taken checking data %ld\nTotal time spent in Get_Data_In_Box function: %ld\nTime spent build recursion frames: %ld\nNumber of clock errors detected: %d\nTotal data points chedcked: %d\nTimes Get_Data_In_Box was called: %d", alloctimetaken, checkDatatime, getDatatime, recursiontime, clockerrs, dataChecked, dataInBoxCount);fflush(stderr);
  PROTECT( ans = NEW_INTEGER( pos ) );
  for (int i = 0; i < pos; i ++)
    INTEGER( ans )[ i ] = found[i] + 1; 
  UNPROTECT(1);
  free(found);
  free(ptrs);
  return ans;
}

void Get_XY_Ptrs( SEXP dattype, SEXP cols, SEXP Rx, SEXP Ry, double **ptrs)
{
  
  if (INTEGER(dattype)[0] == 1)
    {
      //in the case of a matrix/dataframe the number of rows is passed in the Ry slot.
      int len = INTEGER(Ry)[0];
      int xcol = INTEGER(cols)[ 0 ];
      int ycol = INTEGER(cols)[ 1 ];
	
      ptrs[0] = &(REAL( Rx )[ len * ( xcol - 1 ) ] );
      ptrs[1] = &(REAL( Rx )[ len * ( ycol - 1 ) ] );
    } else if (INTEGER(dattype)[0] == 2) 
    {
      //Ry is ignored int he case of the data frame because we can just select the columns we want from the list.
      int xcol = INTEGER(cols)[ 0 ];
      int ycol = INTEGER(cols)[ 1 ];
      
      ptrs[0] = REAL(VECTOR_ELT(Rx, xcol - 1) );
      ptrs[1] = REAL(VECTOR_ELT(Rx, ycol - 1 ) ) ;
    
  } else {
    ptrs[0] = REAL(Rx);
    ptrs[1] = REAL(Ry);
  }

  return;
}

void
Get_Data_In_Box(qtree_t *tree, double left, double right, double down, double up, int **found, int *pos, int *cursize,  double *x, double *y)
{
  int res, res2;
  dataInBoxCount++;
  struct timespec before1, before2;
  struct timespec after1, after2;
  int ind;
  double tmpx, tmpy;
  
  if (tree -> numdata > 0)
    {
      res = clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &before1);
      clockerrs += res ? 1 : 0;
      //if (res)
      //{fprintf(stderr, "clock_gettime returned non-zero: %d", res); fflush(stderr);}
      for(int i =0; i < tree -> numdata; i ++)
	{ 
	  
	  ind = tree -> data[i];
	  tmpx = x[ind];
	  tmpy = y[ind];
	  //if ( x[ ind ] >= left && x[ind] <= right && y[ ind ] >= down && y[ ind ] <= up )	
	    
	    if (tmpx >= left)
	    {
	      if (tmpx <= right)
		{
		  if (tmpy >= down)
		    {
		    if(tmpy <= up)
		      {
			
			if (*pos >= *cursize - 1)
			  {
			    
			    res = clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &before2);
			    clockerrs += res ? 1 : 0;
			    *found = Grow_Return_Array(found, cursize);
			    res = clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &after2);
			    clockerrs += res ? 1 : 0;
			    alloctimetaken += after2.tv_nsec - before2.tv_nsec;
			  }
			
			( *found )[ *pos ] = ind;
     			*pos += 1 ;
			
		      }
		    
		    }
		}
	    }
	}
      
      res2 = clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &after1);
      clockerrs += res2 ? 1 : 0;
      checkDatatime += after1.tv_nsec - before1.tv_nsec;
      dataChecked += tree -> numdata;
    }
  else if (tree -> uleft != NULL)
    {
      /*
      childrenChecked++;
      //double midx = (tree -> left + tree -> right ) / 2.0;
      //double midy = (tree -> down + tree -> up ) / 2.0;
      
      int possible[4] = {1, 1, 1, 1};
      Check_Children(tree, left, right, down, up, possible);

      if (possible[0]) //uleft
	Get_Data_In_Box(tree -> uleft, left, right, down, up, found, pos, cursize, x, y);
      if (possible[1]) //uright
	Get_Data_In_Box(tree -> uright, left, right, down, up, found, pos, cursize, x, y);
      if (possible[2]) //lright
	Get_Data_In_Box(tree -> lright, left, right, down, up, found, pos, cursize, x, y);
      if(possible[3]) //lleft
	Get_Data_In_Box(tree -> lleft, left, right, down, up, found, pos, cursize, x, y);

      */
      qtree_t *tmptree = tree -> uleft;
      if (Check_Bounds(tmptree, left, right, down, up))
	{
	  Get_Data_In_Box(tmptree, left, right, down, up, found, pos, cursize, x, y);
	}
      tmptree = tree -> uright;
      if (Check_Bounds(tmptree, left, right, down, up))
	{
	  Get_Data_In_Box(tmptree, left, right, down, up, found, pos, cursize, x, y);
	}
      tmptree = tree ->lleft;
      if (Check_Bounds(tmptree, left, right, down, up))
	{
	  Get_Data_In_Box(tmptree, left, right, down, up, found, pos, cursize, x, y);
	}
      tmptree = tree -> lright;
      if (Check_Bounds(tmptree, left, right, down, up))
	{
	  Get_Data_In_Box(tmptree, left, right, down, up, found, pos, cursize, x, y);
	}
    }
  return;
}
  
void Check_Children(qtree_t *tree, double left, double right, double down, double up, int *ret)
//ret is an array of length 4 indicating children in this order: uleft, uright, lright, lleft (clockwise) Should be initialized to {1, 1, 1,1 }
{
  double midx = ( tree -> left + tree -> right ) / 2.0;
  double midy = ( tree -> lower + tree -> upper ) / 2.0;

  //for (int i =0; i < 4; i++)
  //  (*ret)[i] = 1;
  if (midx < left)
    {
      ret[0] = 0; //uleft
      ret[3] = 0; //lleft
    }
  if (midx > right)
    {
      ret[1] = 0; //uright
      ret[2] = 0; //lright
    }
  if (midy < down)
    {
      ret[2] = 0; //lright			
      ret[3] = 0; //lleft
    }
  if (midy > up)
    {
      ret[0] = 0; //uleft
      ret[0] = 0; //uright
    }
  return;
}


int * Grow_Return_Array(int **ar, int *size)
{
  int *toret = NULL;
  int oldsize = *size;
  *size = 2 * oldsize;
  toret = realloc( *ar, 2 * oldsize * sizeof(int));
  //toret = memcpy(topret, ar, oldsize * sizeof(int));
  
  return toret;
}

SEXP
R_Get_Top_Left_Order(SEXP Rtree, SEXP Rlen)
{
  int curpos = 0;
  SEXP ans;
  int len = INTEGER( Rlen ) [ 0 ];
  qtree_t *tree = (qtree_t *) R_ExternalPtrAddr( GET_SLOT( Rtree, Rf_install( "ref" ) ) );
  PROTECT( ans = NEW_INTEGER( len ) );
  Get_Data_For_Test( tree, &curpos, ans );
  UNPROTECT(1);
  return ans;
}

void Get_Data_For_Test( qtree_t *node, int *pos, SEXP ans)
{
  if (node -> numdata > 0)
    {
      for (int i =0; i < node -> numdata; i ++)
	{
	  INTEGER(ans) [ *pos ] = node -> data[ i ] + 1;
	  *pos = *pos + 1;
	}
    } else if (node -> uleft != NULL ) {
    Get_Data_For_Test(node -> uleft, pos, ans);
    Get_Data_For_Test(node -> uright, pos, ans);
    Get_Data_For_Test(node -> lright, pos, ans);
    Get_Data_For_Test(node -> lleft, pos, ans);
  }
    
  return;
}

void
R_Free_Quad_Tree( SEXP treeptr)
{
  qtree_t *tree  = R_ExternalPtrAddr(treeptr);
  Free_Quad_Tree(tree);
  return;
}

void Free_Quad_Tree( qtree_t *tree)
{

  if (tree -> parent != NULL)
    {
      qtree_t *par= tree -> parent;
      switch (tree -> pos)
	{
	case 1:
	  par -> uleft = NULL;
	  break;
	case 2:
	  par -> uright = NULL;
	  break;
	case 3:
	  par -> lright = NULL;
	  break;
	case 4:
	  par -> lleft = NULL;
	  break;
	}
    }
  
  if( tree -> uleft != NULL)
    {
      Free_Quad_Tree( tree -> uleft);
      tree -> uleft = NULL;
      Free_Quad_Tree( tree -> uright);
      tree -> uright = NULL;
      Free_Quad_Tree( tree -> lright);
      tree -> lright = NULL;
      Free_Quad_Tree( tree -> lleft);
      tree -> lleft = NULL;
    }
  if (tree -> data != NULL)
    free(tree -> data);

  free( tree );
  return;
}


int
Find_Max_Depth(qtree_t *tree, unsigned char curdepth)
{
  if (tree -> uleft != NULL)
    {
      curdepth = Find_Max_Depth( tree -> uleft, curdepth);
      curdepth = Find_Max_Depth( tree -> uright, curdepth);
      curdepth = Find_Max_Depth( tree -> lright, curdepth);
      curdepth = Find_Max_Depth( tree -> lleft, curdepth);
    } else {
    curdepth = tree->depth > curdepth ? tree->depth : curdepth;
  }
  return curdepth;
}

SEXP
R_Find_Max_Depth(SEXP Rtree)
{
  qtree_t *tree = (qtree_t *) R_ExternalPtrAddr( GET_SLOT( Rtree, Rf_install( "ref" ) ) );
  int toret = Find_Max_Depth( tree, 0) ;
  return ScalarInteger(toret);
}


double 
eucl_Dist(double *pt1, double *pt2, int d)
{
  double tmp = 0;
  double toret = 0;
  for(int i =0; i < d; i++)
    {
      tmp = pt1[i] - pt2[i];
      toret = toret + tmp * tmp;
    }
  return sqrt(toret);
}

void Get_Tree_Attributes(qtree_t *tree, int *curattr)
{
  //first int is number of nodes, second int is number of nodes with data, third int is maxDepth, fourth int is max bucketsize
  curattr[0]++;
  
  
  if(tree -> uleft == NULL)
    {
      if(tree -> numdata > 0)
	{
	  curattr[1]++;
	  curattr[3] = (curattr[3] >= tree -> numdata) ? curattr[3] : tree -> numdata;
	  curattr[2] = (curattr[2] >= tree -> depth) ? curattr[2] : (int) tree -> depth;
	}
    }
  else
    {
      Get_Tree_Attributes(tree -> uleft, curattr);
      Get_Tree_Attributes(tree -> uright, curattr);
      Get_Tree_Attributes(tree -> lright, curattr);
      Get_Tree_Attributes(tree -> lleft, curattr);
    }
     
  return; 
}


/*
double
call_R_Dist(double *pt1, double *pt2, int d, SEXP fun)
{
  int err = 0;
  SEXP Rpt1, Rpt2, call, p, ans;
  PROTECT(Rpt1 = NEW_NUMERIC( d ) );
  PROTECT(Rpt2 = NEW_NUMERIC( d ) );
  for (int i = 0; i < d; i++)
    {
      REAL( Rpt1 ) [ i ] = pt1[ i ];
      REAL( Rpt2 ) [ i ] = pt2[ i ]; 
    }
  
  PROTECT( call = allocVector( LANGSXP, 3 ) );
  SETCAR( call, fun); p = CDR(call);
  SETCAR( p , Rpt1); p = CDR(p);
  SETCAR( p , Rpt2); 
  
  PROTECT( ans = R_tryEval( call, R_GlobalEnv, &err) );
  
  if (err)
    {
      fprintf(stderr, "There was an error evaluating the distance.\n"); fflush(stderr);
    }
  return REAL(ans)[0];
}
  

*/

SEXP
R_LoopTest1()
{
  int *thing = calloc(10000000, sizeof(int));
  int pos = 0;
  for(int i = 0; i < 333333; i++)
    {
      for (int j = 0; j < 33; j++)
	{
	  pos ++;
	  thing[pos] = i;
	}
    }
  free(thing);
  return ScalarLogical(1);
}

    
SEXP
R_LoopTest2()
{
  int pos = 0;
  int *thing = calloc(10000000, sizeof(int));
  for(int i = 0; i < 10000; i++)
    {
      for (int j = 0; j < 1000; j++)
	{
	  
	  thing[pos] = i;
	  pos ++;
	}
    }
  free(thing);
  return ScalarLogical(1);
} 



SEXP
R_LoopTest3()
{
  int statj;
  int *thing = calloc(50*5300000, sizeof(int));
  int pos = 0;
  for(int i = 0; i < 5300000; i++)
    {
      for (statj = 0; statj < 50; statj++)
	{
	  pos ++;
	  thing[pos] = i;
	}
    }
  free(thing);
  return ScalarLogical(1);
}

void Find_Node(qtree_t *node, int val, qtree_t **ret)
{
  if(node -> numdata > 0 )
    {
      for (int i =0; i < node -> numdata; i++)
	{
	  if (node -> data[i] == val)
	    {
	      *ret = node;
	      break;
	    }
	}
    } else if (node -> uleft != NULL) {
    Find_Node(node -> uleft, val, ret);
    Find_Node(node -> uright, val, ret);
    Find_Node(node -> lright, val, ret);
    Find_Node(node -> lleft, val, ret);
  }
  return;
}
