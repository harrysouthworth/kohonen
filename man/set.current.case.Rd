\name{set.current.case}
\alias{set.current.case}
\title{Set a reference prototype}
\description{
Set a reference prototype for the function plot.som with type = 'distances'. 
This cluster will be used as reference for the representation of its distances of all other clusters.
}
\usage{set.current.case( object , current.index )}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a som object}
  \item{current.index}{index of a cluster of the som object.}
}
\value{a som object}
\author{David Gohel}
\seealso{
\code{\link{som}}
\code{\link{plot.som}}
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

lcrabs.som <- set.current.case(lcrabs.som, 4)
plot(lcrabs.som, "distance", cex = .75, cex.label = 0.3 )
}
