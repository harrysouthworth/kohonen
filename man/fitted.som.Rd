\name{fitted.som}
\alias{fitted.som}
\title{Extract som Fitted Values}
\description{The fitted values for the som object are extracted. This is
a method for the fitted generic function for objects of class som.}
\usage{
fitted.som(object)
}
\arguments{
  \item{object}{a som object}
}
\value{a vector with the fitted values.}
\references{ ~put references to the literature/web site here ~ }
\author{David Gohel}
\seealso{
\code{\link{som}}
\code{\link{learn}}
\code{\link{learnBatch}}
\code{\link{plot.som}}
\code{\link{summary.som}}
\code{\link{predict.som}}
\code{\link{biplot.som}}
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

fitted ( lcrabs.som )
}
