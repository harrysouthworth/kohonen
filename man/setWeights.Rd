\name{setWeights}
\alias{setWeights}
\title{Set a matrix as the prototypes' matrix}
\description{
Set a matrix as the prototypes' matrix of a som object.
}
\usage{
setWeights(x, new)
}
\arguments{
  \item{x}{a som object}
  \item{new}{a numeric matrix that have the same dimension that the prototypes' matrix of the som object}
}
\value{
a som object
}
\author{David Gohel}
\seealso{
\code{\link{som}}
\code{\link{learn}}
\code{\link{getWeights}}
}
\examples{
library(MASS)
lcrabs <- log(crabs[, 4:8])

lcrabs.som <- som ( formula = ~ . , data = lcrabs
	, neighborhood = "uniform"
	, grid = grid ( xdim = 10 , ydim = 10 , type = "hexagonal" ) 
	, weights.min = min (lcrabs), weights.max = max (lcrabs)
	)
lcrabs.som <- learn( lcrabs.som , number.iter = 500 , max.alpha = 0.5, min.alpha = .001, max.rayon = 3 , step.eval.si = 50)
plot( lcrabs.som , "effectif" )

x <- getWeights ( lcrabs.som )
# change values of the first prototype
x[1,]=0
lcrabs.som <- setWeights ( lcrabs.som , x )
plot( lcrabs.som , "effectif" )
}
