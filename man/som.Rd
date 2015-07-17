\name{som}
\alias{som}
\title{Creation of a som object}
\description{Create and initialize a som object.}
\usage{
som(formula, data, grid = grid(9, 9, type = "hexagonal"), neighborhood = "gaussian", weights.min = -1, weights.max = 1)
}
\arguments{
  \item{formula}{a formula object, no variables may appear on the left (response) side.}
  \item{data}{a data frame in which to interpret the variables named in the formula.}
  \item{grid}{a grid specifying the structure of the projection space.}
  \item{neighborhood}{a character, possible values are 'uniform' and 'gaussian'.}
  \item{weights.min}{numeric value, minimal value to used during the creation of the initial weights.}
  \item{weights.max}{numeric value, maximal value to used during the creation of the initial weights.}
}
\details{
Scale data to mean = 0 and std = 1 is recommanded.
If the model contains factors, they will be transformed with the contr.treatment function.
}
\value{
object of class \code{{som}}: a list with several components:
\item{data }{model.matrix of the model}
\item{withinss }{within sums of squares for each row}
\item{weights }{weights of the network}
\item{grid }{argument \code{{grid}}}
\item{size }{vector of length \code{{nrow(weights)}}. The ith value gives the number of data points in cluster i.}
\item{neighborhood }{argument \code{{neighborhood}}}
\item{groups }{vector of length \code{{nrow(data)}}. The ith value gives the cluster of the data points i}
\item{metagroups }{vector of length \code{{nrow(data)}}. The ith value gives the meta-cluster of the data points i}
}

\references{
	Kohonen, T. (1995).
	\emph{Self-Organizing Maps}
}
\author{David Gohel}
\note{
\section{side effects}{
Meta-groups are all equals to 1, the same meta-group.
}
\section{Warning}{
NA values are not allowed.
}
}
\seealso{
\code{\link{grid}}
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
lcrabs.som
}
