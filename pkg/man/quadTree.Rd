\name{quadTree}
\alias{quadTree}
\title{
Create a Quad Tree for 2-Dimensional Data
}
\description{
This function creates a Quad Tree recursive partitioning of the data. This is useful for spatial lookups, etc.
}
\usage{
quadTree(x, y)
}
\arguments{
  \item{x}{
  Vector of x coordinates in 2 dimensional space. Will be converted to numeric.
}
  \item{y}{
  Vector of y coordinates in 2 dimensional space. Will be converted to numeric.
}
}
\details{
}
\value{
  The function returns an object of S4 class QuadTree, which has the following slots:
  \item{name}{Unused.}
  \item{ref}{External pointer to the underlying C data structure representing the QuadTree (ie the root node).}
  \item{maxDepth}{Integer representing the maximum depth of the tree.}
  \item{points}{Integer indicating the number of points stored in the tree, ie the number of points passed into quadTree}
}
\references{
}
\author{
Gabriel Becker
}
\note{
}
\seealso{
}
\examples{
x = rnorm(1000)
y = rnorm(1000)
tree = quadTree( x , y )
}
\keyword{spatial}
\keyword{storage}