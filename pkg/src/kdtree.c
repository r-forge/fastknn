#include "Rsearchtrees.h"

void R_Free_KD_Tree(SEXP Rtree);
kdtree_t *Build_kdTree( int *oldorder, int npts, double *data, kdtree_t *parent, int k, int n, int *allorders);
void Get_Left_Order(kdtree_t *tree, int *dats, int *pos);
SEXP R_Get_Left_Order(SEXP Rtree);

SEXP
R_Build_kdTree(SEXP Rdata, SEXP Rorders, SEXP Rk, SEXP Rn)
{
  int k = INTEGER( Rk ) [ 0 ];
  int *orders = INTEGER( Rorders );
  double *data = REAL( Rdata );
  int n = INTEGER( Rn ) [ 0 ];
  
  //kdtree_t *tree = Build_kdTree(&orders[0], data, NULL, k, n, 0, n - 1, orders);
  kdtree_t *tree = Build_kdTree(orders, n, data, NULL, k, n,  orders);

  SEXP klass, ans, ptr;
  PROTECT( klass = MAKE_CLASS( "KDTree" ) );
  PROTECT( ans = NEW( klass ) );
  PROTECT( ptr = R_MakeExternalPtr( tree ,
				    Rf_install( "KDTree" ),
				    R_NilValue ) );
    R_RegisterCFinalizerEx(ptr, &R_Free_KD_Tree, 1);
  SET_SLOT( ans, Rf_install( "ref" ), ptr );
  //SET_SLOT( ans, Rf_install( "maxDepth" ), ScalarInteger( Find_Max_Depth( tree, 0) ) );
  UNPROTECT(3);
  return ans;

}

 
//kdtree_t *Build_kdTree( int *oldorder, double *data, kdtree_t *parent, int k, int n, int startpos, int endpos, int *allorders)
kdtree_t *Build_kdTree( int *oldorder, int npts, double *data, kdtree_t *parent, int k, int n, int *allorders)
{
  
  //just incase;
  if (npts < 1)
    return NULL;
  int depth;
  if (!parent)
    depth = 0;
  else
    depth = parent -> depth + 1;

  int col = depth % k;
  //int npoints = endpos - startpos + 1; //0-4: 0 1 2 3 4: 5pts;
       
  int neworder[npts];
  int datinds[npts];
  int found = 0, pos=0;

  //XXX this seems cumbersome, is there a better way?
  /* 
 while ( found < npts && pos < n)
    {
      for(int j = 0; j< npts; j++)
	{
	  if(allorders[ n * ((col) % k) + pos ] ==  oldorder [ j ])
	    {
	      datinds[found] = pos;
	      //neworder[found] = oldorder[j];
	      //neworder[found] = allorders[ n * col + pos ];
	      neworder[found] = allorders[ n * ((col + 1) % k ) + pos];
	      found++;
	      break;
	    }
	}
      pos++;
    }
  */
  for (int i =0; i < npts; i++)
    neworder[i] = allorders[ n * ( ( col +1 )% k ) + oldorder[ i ] - 1 ] ;
     
  kdtree_t *tree = malloc(sizeof(kdtree_t));
  tree -> depth = depth;
  tree -> parent = parent;
  if (npts == 1)
    {
      tree -> data = oldorder[0];
      tree -> left = NULL;
      tree -> right = NULL;
    } else {

    int med = npts /2 ; //integer division!!
    double medval;
    //the orders, oldorder, neworder, etc were caluclated in R!! this means they start at 1, not 0!! 
    if (npts % 2 == 0)
      medval = (data[ n * col + oldorder[ med ] - 1 ] + data[ n* col + oldorder[ med - 1 ] - 1] ) / 2.0;
    else
      medval = data[ n * col + oldorder[ med ] - 1 ];
	
    int forleft = 0;
    for (int i = 0; i < found; i++)
      {
	if (data[ n * col + oldorder[ i ] - 1 ] <= medval)
	  forleft++;
      }

  
	
    tree -> split = medval;
    
    //tree -> left = Build_kdTree(&nextcol[0], data, tree, k, n,
    //				startpos, startpos + med - 1, allorders);
    tree -> left = Build_kdTree(neworder, forleft, data, tree, k, n, allorders);
    
    tree -> right = Build_kdTree( &neworder[ forleft + 1 ], found - forleft, data, tree, k, n, allorders); 
    
    //tree -> right = Build_kdTree(&nextcol[ med + 1], data, tree, k, n, startpos + med + 1, endpos, allorders);
  }
  return tree;
    
}

SEXP 
R_Get_Left_Order(SEXP Rtree)
{
  kdtree_t *tree = R_ExternalPtrAddr( GET_SLOT( Rtree, Rf_install( "ref" ) ) );
  int npts = INTEGER( GET_SLOT( Rtree , Rf_install( "npoints" ) ) )[0];
  int toret[npts];

  int pos = 0;
  if (npts)
    Get_Left_Order(tree, toret, &pos);
  SEXP ans;
  PROTECT( ans = NEW_INTEGER( npts ) );
  for(int i =0; i < npts; i++)
    INTEGER(ans)[i] = toret[i] + 1;
  UNPROTECT(1);
  return ans;
}

void 
Get_Left_Order(kdtree_t *tree, int *dats, int *pos)
{

  //unlike quadtree, ever node has data.
  dats[ *pos ] = tree -> data + 1;
  *pos = *pos + 1;
  if (tree -> left != NULL)
    Get_Left_Order(tree -> left, dats, pos);
  
  if (tree -> right != NULL)
    Get_Left_Order(tree -> right, dats, pos);

  return;
}

void
R_Free_KD_Tree(SEXP Rtree)
{
  return;
}
