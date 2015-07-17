\name{summary.som.cluster}
\alias{summary.som.cluster}
\title{Summary Method for som.cluster object.}
\description{Returns a summary list for a som.cluster object.}
\usage{summary.som.cluster(object)}
\arguments{
  \item{object}{a som object}
}
\value{
object of class \code{{summary.som.cluster}} : a list with several components:
\item{data}{subset of the model.matrix of the model that belong to the class}
\item{weights}{prototype of the class}
\item{class.number}{integer representing the number id of the class}
\item{withinss}{mean of within sum of squares in the class}
\item{ss}{within sum of squares for each data row in the class}
\item{outliers}{identified outliers (boxplot method on the ss values)}
\item{stdevss}{Standart deviation of the sum of squares)}
}
\author{David Gohel}
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

summary ( get.cluster( lcrabs.som , 6 ) )
}
