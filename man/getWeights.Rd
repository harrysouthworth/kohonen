\name{getWeights}
\alias{getWeights}
\title{Extraction of the prototypes' matrix}
\description{
Extract the matrix of prototypes
}
\usage{
getWeights(x)
}
\arguments{
  \item{x}{a som object}
}
\value{
the matrix of prototypes
}
\author{David Gohel}
\seealso{
\code{\link{som}}
\code{\link{learn}}
\code{\link{setWeights}}
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
getWeights ( lcrabs.som )
}
