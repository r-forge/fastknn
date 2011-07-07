\name{getPointsInRect}
\alias{getPointsInRect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getPointsInRect(tree, left, right, down, up, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
%%     ~~Describe \code{tree} here~~
}
  \item{left}{
%%     ~~Describe \code{left} here~~
}
  \item{right}{
%%     ~~Describe \code{right} here~~
}
  \item{down}{
%%     ~~Describe \code{down} here~~
}
  \item{up}{
%%     ~~Describe \code{up} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
structure(function (tree, left, right, down, up, data) 
standardGeneric("getPointsInRect"), generic = structure("getPointsInRect", package = "QuadTreeKNN"), package = "QuadTreeKNN", group = list(), valueClass = character(0), signature = c("tree", 
"left", "right", "down", "up", "data"), default = `NULL`, skeleton = function (tree, 
    left, right, down, up, data) 
stop("invalid call in method dispatch to \"getPointsInRect\" (no default method)", 
    domain = NA)(tree, left, right, down, up, data), class = structure("standardGeneric", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
