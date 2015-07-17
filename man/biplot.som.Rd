\name{biplot.som}
\alias{biplot.som}
\title{Biplot of a som object}
\description{Produces a biplot corresponding to the som object and its data.}
\usage{biplot.som(object, type = "lines", cex = 0.75)}
\arguments{
  \item{object}{a som object}
  \item{type}{character specifying the type of biplot : [lines | bar | boxplot | solid]}
  \item{cex}{character expansion of titles and legend.}
}
\details{
type = 'lines' : for each class, plot values and the associated prototype,
type = 'solid' : for each class, plot the prototype as a solid polygon.
type = 'stars' : for each class, plot the prototype as a star.
type = 'boxplot' : plot a boxplot of the sums of squares (ss of the data in the class).
}
\note{When type = 'lines', the plotting procedure is slow.}
\author{David Gohel}
\seealso{
\code{\link{som}}
\code{\link{learn}}
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

#--- profil for each class ---#
biplot( lcrabs.som , type = "stars")
biplot( lcrabs.som , type = "lines")
biplot( lcrabs.som , type = "solid")
biplot( lcrabs.som , type = "boxplot")
}
\keyword{multivariate}
\keyword{hplot}
