\name{update.som}
\alias{update.som}
\title{Update a som object}
\description{
Update a som object, compute statitics and can also project a new dataset
}
\usage{
update.som(object, newdata)
}
\arguments{
  \item{object}{a som object}
  \item{newdata}{a data frame in which to interpret the variables of the som object.}
}
\details{\code{groups}, \code{size}, \code{withinss} are re-computed.}
\value{a som object.}
\author{David Gohel}
\note{
This function is useful if you want to project a new dataset with the som object or if you have modified by the hand the \code{weights} matrix.
}
\seealso{
\code{\link{som}}
\code{\link{learn}}
\code{\link{learnBatch}}
\code{\link{plot.som}}
\code{\link{summary.som}}
\code{\link{predict.som}}
\code{\link{biplot.som}}
\code{\link{getWeights}}
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

# see the counts' map with learning data
plot( lcrabs.som , "effectif", cex.label = 0)

# see the counts' map with new data
lcrabs.som <- update( lcrabs.som , newdata = lcrabs[1:20,] )
plot( lcrabs.som , "effectif", cex.label = 0)

}
