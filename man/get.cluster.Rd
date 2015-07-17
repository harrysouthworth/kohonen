\name{get.cluster}
\alias{get.cluster}
\title{Extract data and prototype for a specified class.}
\description{Extract data and prototype for a specified class. The function return an object of class 'som.cluster'.}
\usage{
get.cluster(object, class.number)
}
\arguments{
  \item{object}{a som object}
  \item{class.number}{integer representing the number id of the class.}
}
\details{
The function allow to browse a specific class.
}
\value{
object of class \code{{som.cluster}} : a list with several components:
	\item{data}{subset of the model.matrix of the model that belong to the class}
	\item{weights}{prototype of the class}
	\item{class.number}{integer representing the number id of the class}
	\item{withinss}{mean of within sum of squares in the class}
	\item{ss}{within sum of squares for each data row in the class}
}
\author{David Gohel}
\seealso{
\code{\link{som}}
\code{\link{plot.som.cluster}}
\code{\link{summary.som.cluster}}
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

x <- get.cluster( lcrabs.som , 1)
x
plot( x )
}
