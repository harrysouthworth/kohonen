\name{plot.som}
\alias{plot.som}
\title{Plot a som object}
\description{
Plot a som object, several plot can be drawn, see below.
}
\usage{
plot.som( object, type = "effectif", cex = 0.75, nom.variable , precision = 4, cex.label = 0.75 )
}
\arguments{
  \item{object}{a som object}
  \item{type}{character specifying the type of the graph. Possible values are : [distance | mean.ss | effectif | meta | variable | energy | dxdy] (see details)}
  \item{cex}{character expansion of titles and legend.}
  \item{nom.variable}{It has no effect if type is different from 'variable'. Name of the variable of weights to plot.}
  \item{precision}{Precision of the displayed label if it is numeric}
  \item{cex.label}{character expansion of the displayed info in each polygon of the network.}
}
\details{
\item{distance}{euclidean distances between a class and all the other classes. Change the reference class with <s-expression>set.current.case</s-expression>.}
\item{mean.ss}{mean sums of squares for each prototype.}
\item{effectif}{number of rows classified in each prototype.}
\item{meta}{meta-groups (see <s-expression>plot.clust</s-expression> for precisions)}
\item{variable}{weights values for a variable in each prototype. It allows to see if a variable if used by the network or not (if not, values are quite equals).}
\item{energy}{plot the energy function during the learning process. This represent the intra-inertia extended to the neighborhood at a stage of the learning process. This is the objective function to be minimised during learning}
\item{dxdy}{The quality of a projection can be evaluated by the 'dy-dx' representation. It is a plot of all the possible distances in the input space, dx's, versus their respective distances in the output space, dy's. For a linear projection the 'dy-dx' plot should be linear.}

}
\references{
	Demartines, P. and J. Herault.
	\emph{Curvilinear component analysis: A self-organizing neural network for nonlinear mapping of data sets.}
	Kohonen, T. (1995).
	\emph{Self-Organizing Maps}
}
\author{David Gohel}
\seealso{
\code{\link{som}}
\code{\link{learn}}
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



plot(lcrabs.som, "energy")
#--- quality of learn, points must be higly correlated and in a line from c(0, 0) -> c(1, 1)
plot ( lcrabs.som , "dxdy" )
#--- plot number of cases per neuron/class
plot(lcrabs.som, "effectif", cex.label = 0)
#--- plot means sums of squares
plot(lcrabs.som, "mean.ss", cex.label = 0.3)
#--- change current case
lcrabs.som <- set.current.case(lcrabs.som, 4)
plot(lcrabs.som, "distance", cex = .75, cex.label = 0)



#--- plot values for a variable... useful for selecting active varaible
plot(lcrabs.som, "variable")
plot(lcrabs.som, "variable", nom.variable = "FL")

\dontrun{
#--- construct meta class , click on the graph, the y will be used to collapse classes in meta-classes ---#
lcrabs.som <- plot.clust(lcrabs.som, interactive = T)
}

#--- plot meta class ---#
plot(lcrabs.som, "meta" , cex.label = .5 )

}
\keyword{hplot}
