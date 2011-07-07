\name{getPointsInBox}
\alias{getPointsInBox}
\title{
Lookup Points Falling in a Rectangular Region in 2 Dimensions
}
\description{
This function determines which points used to create a QuadTree object fall within a specified rectangular region.
}
\usage{
getPointsInBox(tree, left, right, down, up, x, y)
}
\arguments{
  \item{tree}{
A QuadTree object created using all possible points.
}
  \item{left}{
  The leftmost (lower limit in the x direction) boundary of the rectangle.
}
  \item{right}{
  The rightmost (upper limit in the x direction) boundary of the rectangle.
}
  \item{down}{
  The lower (in the y direction) boundary of the rectangle.
}
  \item{up}{
  The upper (in the y direction) boundary of the rectangle.
}
  \item{x}{
  The x data used to create the QuadTree object. See note for why this is necessary.
}
  \item{y}{
  The y data used to create the QuadTree object. See note for why this is necessary.
}
}
\details{
}
\value{
  The function returns an integer vector containing the indices (row numbers) of the points which fall within the specified rectangle.
}
\references{
}
\author{
Gabriel Becker
}
\note{
  Passing in the original data is necessary because to save on storage, the tree is built to store only row indices, not the actual (x,y) data.,
}
\seealso{
}
\examples{
  #generate the data.
  x = rnorm(100)
  y = rnorm(100)
  
  #build the tree
  tree = quadTree(x,y)

  #our rectangle has corners at (0, 0) and (1, 2)
  inds = getPointsInBox(tree, 0, 1, 0, 2, x, y)
  
  #now we plot it
  ch = rep(1, times=100)
  ch[inds] = 3
  plot(x, y, pch = ch)
  lines(c(0, 0, 1, 1, 0), c(0, 2, 2, 0, 0))

}
\keyword{lookup}