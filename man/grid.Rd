\name{grid}
\alias{grid}
\title{Creation of a grid}
\description{
Creation of an object that defines the map structure of a som object.
}
\usage{
grid(xdim, ydim, type = "hexagonal", current.index = 1 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xdim}{number of rows of the grid}
  \item{ydim}{number of cols of the grid}
  \item{type}{a character specifying the type of the grid, 'rectangular' or 'hexagonal'.}
  \item{current.index}{cluster index used as reference in the distance plot.}
}
\value{
a grid object
}
\author{David Gohel}
\seealso{
\code{\link{som}}
\code{\link{learn}}
\code{\link{learnBatch}}
\code{\link{plot.som}}
\code{\link{biplot.som}}
}
\examples{
library(MASS)
lcrabs <- log(crabs[, 4:8])
x <- grid ( xdim = 20 , ydim = 13 , type = "hexagonal" )
lcrabs.som <- som ( formula = ~ . , data = lcrabs
	, neighborhood = "uniform"
	, grid =  x
	, weights.min = 0 , weights.max = 4
	)
}
