\name{predict.som}
\alias{predict.som}
\title{Predictions from a som Object}
\description{Predictions are obtained from the minimum distance between a prototype and a row of the data.}
\usage{predict.som(object, newdata, type = "class")}
\arguments{
  \item{object}{a som object}
  \item{newdata}{a data frame in which to interpret the variables of the som object.}
  \item{type}{Output type (character), [ 'class' | 'meta' | 'dist'].}
}
\value{
type = 'class' gives the class id for each row.
type = 'meta' gives the meta-class id for each row.
type = 'dist' gives the entire distances matrix.
}
\author{David Gohel}
\seealso{
\code{\link{som}}
\code{\link{learn}}
\code{\link{update.som}}
\code{\link{summary.som}}
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

predict( lcrabs.som)
predict( lcrabs.som , type = "dist" )
}
