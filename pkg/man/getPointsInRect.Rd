\name{getPointsInRect}
\alias{getPointsInRect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform Rectangular Range Lookup on a Tree Object
}
\description{
Determine which two dimensional points (which have already been stored
in a search tree) fall within a specified rectangle.
}
\usage{
getPointsInRect(tree, left, right, down, up, data)
}
\arguments{
  \item{tree}{
    S4 object representing the tree to be searched over.
}
\item{left}{
  Left bound of the rectangular range. Will be converted to numeric.
}
\item{right}{
Right bound of the rectungular range. Will be converted to numeric.
}
\item{down}{
  Lower bound of the rectangular range. Will be converted to numeric.
}
\item{up}{
Upper bound of the rectangular range. Will be converted to numeric.
}
  \item{data}{
Data used to create the tree.
}
\item{columns}{
  Columns in the \code{data} object which contain the relevant
  data. Defaults to columns 1 and 2.
}
}
\details{
}
\value{
  An integer vector is returned containing the indices (row numbers) of
  points which fall within the specified rectangle.
}
\references{
}
\author{
Gabriel Becker
}
\seealso{
  \link{\code{findKNN} } 
}
\examples{
x = rnorm(100)
y = rnorm(100)
dat = cbind(x,y)
tree = createIndex(dat) #default tree type is "quad"
#our rectangle has corners at (0,0), (0,2), (2, 2), and (2,0)
inds = getPointsInRect(tree, 0, 2, 0, 2, dat)

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
