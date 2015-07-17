\name{summary.som}
\alias{summary.som}
\title{Summary Method for som objects}
\description{Returns a summary list for a som object. This is a summary method
for objects of class som.}
\usage{
summary.som(object)
}
\arguments{
  \item{object}{a som object}
}
\details{
intra inertia + inter inertia gives the total inertia. 
Inter inertia is maximized and intra inertia is minimized.
}
\value{
\item{inertia}{inter and intra inertia}
\item{weights}{weights matrix of the network}
\item{withinss}{sums of squares for each cluster}
\item{size}{vector of length \code{nrow(weights)}. The ith value gives the number of data points in cluster i.}
\item{class}{vector of length \code{nrow(data)}. The ith value gives the cluster of the data points i.}
\item{metaclass}{vector of length \code{nrow(data)}. The ith value gives the meta-cluster of the data points i.}
}
\author{David Gohel}
\seealso{
\code{\link{som}}
\code{\link{learn}}
\code{\link{learnBatch}}
}
\examples{
library(MASS)
lcrabs <- log(crabs[, 4:8])

lcrabs.som <- som ( formula = ~ . , data = lcrabs
	, neighborhood = "uniform"
	, grid = grid ( xdim = 10 , ydim = 10 , type = "hexagonal" ) 
	, weights.min = min (lcrabs), weights.max = max (lcrabs)
	)
lcrabs.som <- learn( lcrabs.som , number.iter = 500 , max.alpha = 0.3, step.eval.si = 50)

summary( lcrabs.som )
}
