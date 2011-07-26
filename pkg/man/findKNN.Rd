\name{findKNN}
\docType{methods}
\alias{findKNN}
\alias{findKNN,QuadTree-method}
\title{
  Perform k-Nearest Neighbors Lookup Using a Search Tree
}
\description{
  This function performs fast Nearest Neighbors lookup on an object
  created using \code{createIndex}
}
\usage{
findKNN(tree, newdat, fulldat, newcols = 1:2, fullcols = 1:2, k = 5)
}
\arguments{
  \item{tree}{
    An object which inherits from the \code{SearchTree} S4 class.
  }
  \item{newdat}{
    A matrix or data.frame containing the data for the points to lookup
    neighbors for.
  }
  \item{fulldat}{
    A matrix or data.frame containing the data used to create \code{tree}  
  }
  \item{newcols}{
    The columns within \code{newdat} where the relevant
    data can be found. Defaults to columns 1 and 2.
  }
  \item{fullcols}{
    The columns within \code{fulldat} where the relevant data can be
    found. Defaults to columns 1 and 2.
  }
  \item{k}{
    Number of Nearest Neighbors to find for each new point in
    \code{newdat}
  }
}
\value{
  The return value is an integer matrix indicating the indices in
  \code{fulldat} where the nearest neighbors were found. Row indicates
  point in \code{newdat}, while column indicates the order of the k
  neighbors.
}
\author{
  Gabriel Becker
}
\note{
  If given a point in \code{newdat} that is in the original dataset
  \code{fulldat}, the first neighbor will always be the point itself.
}
\seealso{
  \code{\link{createIndex}}, \code{\link{getPointsInRect}}
}
\examples{
x= rnorm(100)
y = rnorm(100)
dat = cbind(x,y)
newdat = cbind(rnorm(2), rnorm(2))
tree = createIndex(dat)
inds = findKNN(tree, newdat, dat)

ch = rep(1, times=100)
ch[inds[1:5]] = 3
ch[inds[6:10]] = 5
cls = rep("black", times=100)
cls[inds[1:5]] = "red"
cls[inds[6:10]] ="blue"

plot(x,y, pch=ch, col = cls)
abline(v=newdat[1,1], h = newdat[1, 2], col="red")
abline(v=newdat[2,1], h = newdat[2,2], col = "blue")
}
\keyword{lookup}
\keyword{knn}
