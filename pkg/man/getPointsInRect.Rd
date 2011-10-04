\name{getPointsInRect}
\docType{methods}
\alias{getPointsInRect}
\alias{getPointsInRect,QuadTree-method}
\title{
Perform Rectangular Range Lookup on a Tree Object
}
\description{
  Determine which two dimensional points (which have already
  been stored in a search tree) fall within a specified
  rectangle.
}
\usage{
getPointsInRect(tree, ptOne, ptTwo, data, columns)
}
\arguments{
  \item{tree}{
    S4 object representing the tree to be searched over.
  }
  \item{ptOne}{
    Numeric vector of length 2 indicating x and y coordinates of first
    point defining a rectangle.
  }
  \item{ptTwo}{
    Numeric vector of length 2 indicating x and y coordinates of second
    point defining a rectangle. Must be non-colinear with ptOne in both
    the x and y directions. 
  }
  
  \item{data}{
    Data used to create the tree.
  }
  \item{columns}{
    Columns in the \code{data} object which contain
    the relevant data. 
  }
}
\value{
  An integer vector is returned containing the indices (row numbers) of
  points which fall within the specified rectangle.
}
\author{
  Gabriel Becker
}
\seealso{
  \code{\link{findKNN}},\code{\link{createIndex}} 
}
\examples{
x = rnorm(100)
y = rnorm(100)
dat = cbind(x,y)
tree = createIndex(dat) #default tree type is "quad"
#our rectangle has corners at (0,0), (0,2), (2, 2), and (2,0)
inds = getPointsInRect(tree, c(0, 0), c(2, 2), dat)

#prepare vector of plotting characters
ch = rep(1, times = 100)
ch[inds] = 3

#plot points indicating which ones were selected as in our rectangle
plot(x,y, pch = ch)
#add rectangle to plot
xpts = c(0, 0, 2, 2, 0)
ypts = c(0, 2, 2, 0, 0)
lines(xpts, ypts)
}
\keyword{lookup}
\keyword{search}
