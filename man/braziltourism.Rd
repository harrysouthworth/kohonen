\name{braziltourism}
\alias{braziltourism}
\docType{data}
\title{ ~~ data name/kind ... ~~}
\description{
  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(braziltourism)}
\format{
  A data frame with 412 observations on the following 9 variables.
  \describe{
    \item{\code{Age}}{a numeric vector}
    \item{\code{Sex}}{a factor with levels \code{Female} \code{Male}}
    \item{\code{Income}}{a numeric vector}
    \item{\code{Travel.cost}}{a numeric vector}
    \item{\code{Access.road}}{a numeric vector}
    \item{\code{Active}}{a numeric vector}
    \item{\code{Passive}}{a numeric vector}
    \item{\code{Logged.income}}{a numeric vector}
    \item{\code{Trips}}{a numeric vector}
  }
}
\details{
This dataset is used to illustrate SOM with categorical and numerical data.
}
\source{
	Jeffrey S. Simonoff
	\emph{Analyzing Categorical Data}
}

\examples{
braziltourism2 <- braziltourism
braziltourism2[, sapply( braziltourism2 , is.numeric ) ] <- sapply ( braziltourism2 [, sapply( braziltourism2 , is.numeric ) ] , scale )
bt.som.init  <- som ( ~ Age + Sex + Income + Travel.cost + Active + Passive + Logged.income + Trips
	, data = braziltourism2
	, grid = grid ( xdim = 15 , ydim = 15 , type = "rectangular" ) 
	)
bt.som <- learn( bt.som.init , number.iter = 5000, max.alpha = 0.5, min.alpha = .001, step.eval.si = 100)
plot(bt.som, "energy")
plot(bt.som, "effectif", cex = .75 , cex.label = .0 )}
\keyword{sysdata}
