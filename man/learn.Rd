\name{learn}
\alias{learn}
\title{Train the network (Online algorithm).}
\description{
Perform the learning process of the network with online algorithm.
}
\usage{
learn( object, number.iter, order, max.alpha = 0.5, min.alpha = 0.0, max.rayon = 3, min.rayon = 0, step.eval.si = 100 )
}
\arguments{
  \item{object}{a som object}
  \item{number.iter}{number of learning process' step. It won't have effect if \code{order} is specified.}
  \item{order}{vector of length \code{number.iter}, specifying which row is presented to the learning process at the n-ith step.}
	\item{max.alpha}{the maximum amount of change during learning process.}
	\item{min.alpha}{the minimum amount of change during learning process.}
	\item{max.rayon}{the maximum radius during learning process.}
	\item{min.rayon}{the minimum radii during learning process.}
	\item{step.eval.si}{integer. A vector will be created, it will have values of evenly spaced numbers from 1 to \code{number.iter} with spacing \code{step.eval.si}. At each value of the vector, an evaluation of energy function will be performed.}
}
\details{
Be carefull with \code{step.eval.si}, the smallest is the value, the longer it will take time to process.
alpha decrease linearly from \code{max.alpha} to \code{min.alpha} over \code{number.iter} updates.
rayon decrease linearly from \code{max.rayon} to \code{min.rayon} over \code{number.iter} updates.
}
\value{
a som object.
}
\references{
	Kohonen, T. (1995).
	\emph{Self-Organizing Maps}
}
\author{David Gohel}
\seealso{
\code{\link{som}}
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
# train the network
lcrabs.som <- learn( lcrabs.som , number.iter = 500 , max.alpha = 0.3, step.eval.si = 50)
# plot energy function and print a summary
plot( lcrabs.som , "energy")
summary( lcrabs.som )

}
\keyword{learn}
