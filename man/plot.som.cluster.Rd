\name{plot.som.cluster}
\alias{plot.som.cluster}
\title{plot a som.cluster object}
\description{plot a som.cluster object}
\usage{plot.som.cluster( object )}
\arguments{
  \item{object}{a som object}
}
\author{David Gohel}
\note{Values are not scaled before plotting. Be sure that all data have the same range because, if not, graphics won't allow to see correctly everything.}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
library(MASS)
lcrabs <- log(crabs[, 4:8])

lcrabs.som <- som ( formula = ~ . , data = lcrabs
	, neighborhood = "uniform"
	, grid = grid ( xdim = 10 , ydim = 10 , type = "hexagonal" ) 
	, weights.min = min (lcrabs), weights.max = max (lcrabs)
	)
lcrabs.som <- learn( lcrabs.som , number.iter = 500 , max.alpha = 0.5, min.alpha = .001, max.rayon = 3 , step.eval.si = 50)

plot ( get.cluster( lcrabs.som , 6 ) )
}
