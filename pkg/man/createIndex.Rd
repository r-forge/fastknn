\name{createIndex}
\alias{createIndex}
\title{
Create a Search Tree Index
}
\description{
  Create a search tree from the supplied data for use
  in during futurelookups.
}
\usage{
createIndex(data, type = "quad", columns = c(1, 2), ...)
}
\arguments{
  \item{data}{
    Data to create the tree from. Must be a matrix or data.frame like object.
  }
  \item{type}{
    Type of tree structure to create. Currently only "quad" is supported
  }
  \item{columns}{
    Which columns in \code{data} to use to create the tree. Defaults to
    columns 1 and 2.
  }
  \item{\dots}{
    Any additional/type specific parameters to be passed to the tree
    creation function.
  }
}
\value{
  The class of the returned object depends on the tree type created,
  though all will inherit from the \code{SearchTree} S4 class and have the
  following slots:
  \item{ref}{An external pointer to the C level data structure.}
  \item{numNodes}{Total number of nodes comprising the tree.}
  \item{dataNodes}{Number of terminal/leaf nodes which store at least
    one data point.}
  \item{maxBucket}{Maximum number of data points stored in a single
    terminal/leaf node.}
  \item{points}{Number of points used to create the tree.}
}
\references{
  Finkel, R. A. and Bentley, J. L. "Quad Trees, a Data Structure for Retrieval on Composite Keys." Acta Informatica 4, 1-9, 1974.}
\author{
  Gabriel Becker
}
\seealso{
  \code{\link{findKNN}}, \code{\link{getPointsInRect}}
}
\examples{
 x = rnorm(100)
 y = rnorm(100)
 dat = cbind(x,y)
 #maxDepth is a parameter settable for quad trees
 tree = createIndex(dat, maxDepth = 10, type="quad")
}
\keyword{storage}
\keyword{lookup}
\keyword{indexing}