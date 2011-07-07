\name{getKNNIndices}
\alias{getKNNIndices}
\title{
Lookup k Nearest Neighbors using a QuadTree
}
\description{
This function looks up the indices of the nearest neighbors to a give set of new points.
}
\usage{
getKNNIndices(tree, newx, newy = NULL, allx, ally, k = 5)
}
\arguments{
  \item{tree}{
  A QuadTree object created using the training, or existing, dataset.
}
  \item{newx}{
  A vector of the x coordinates of points to look up neighbors of. Alternatively, a matrix containing both x and y coordinates of the new points, see details.
}
  \item{newy}{
  A vector of y coordinates for the points to look up neighbors for, or NULL
}
  \item{allx}{
The x coordinates of the points used to create the QuadTree.
}
  \item{ally}{
  The y coordinates of the points used to create the QuadTree.
}
  \item{k}{
  Number of neighbors to find for each new point. Defaults to 5.
}
}
\details{

  If newy is null, then newx is checked for having at least 2 columns. If it does, newy is set to be the second column, while newx is set to the first. All other columns are ignored.
}
\value{
  If a single point is specified the function returns an integer vector indicating the indices (row numbers) of the original data (specified by allx, and ally) which represent the nearest neighbors to the new point.
  
  If multiple new points are specified, an n x k (n being the number of new points) integer matrix is returned. Each row represents the k nearest neighbors of the corresponding point.

}
\references{
}
\author{
Gabriel Becker
}
\note{
  As with other methods using QuadTrees, this method requires the original data to be passed to it as well as the new data because the quadtree stores only row numbers, as opposed to the actual data, in order to conserve memory for large data sets.
}
\seealso{
\code{\link{getPointsInBox}}
}
\examples{
  #generate data
  x = runif(100, 0, 10)
  y = runif(100, 0, 10)
  
  #create tree
  tree = quadTree(x, y)
  
  #get the 6 nearest neighbors to the point (4, 5)
  inds = getKNNIndices(tree, 4, 5, x, y, 6)
  
  #the proof is in the plot
  ch = rep(1, times=100)
  ch[inds] = 3
  plot(x, y, pch = ch)
  abline(h = 5, v = 4)
}
\keyword{KNN}
\keyword{lookup}