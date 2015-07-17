\name{learnBatch}
\alias{learnBatch}
\title{Train the network (Batch algorithm).}
\description{
Perform the learning process of the network with Batch algorithm.
}
\usage{
learnBatch( object , number.iter , max.rayon = 3 , min.rayon = 0 )
}
\arguments{
  \item{object}{a som object}
  \item{number.iter}{number of learning process' step. It won't have effect if \code{order} is specified.}
  \item{max.rayon}{the maximum radius during learning process.}
  \item{min.rayon}{the minimum radii during learning process.}
}
\details{
Be carefull with \code{number.iter}, it should be small.
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
\code{\link{learn}}
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
lcrabs.som <- learnBatch( lcrabs.som , number.iter = 10 )
# plot energy function and print a summary
plot( lcrabs.som , "energy")
summary( lcrabs.som )

}
