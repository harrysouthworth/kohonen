\name{authorship}
\alias{authorship}
\docType{data}
\title{Word counts for books by various author}
\description{
The authorship data set is a data.frame of 69 numeric columns and 2 factor columns.
Each row contains data giving counts of each word. The last two columns are providing the id of the book and its author. 
}
\usage{data(authorship)}
\format{
  A data frame with 840 observations on the following 71 variables.
  \describe{
    \item{\code{str.a}}{a numeric vector that represent the number of times that the string 'a' appeared}
    \item{\code{str.all}}{a numeric vector that represent the number of times that the string 'all' appeared}
    \item{\code{str.also}}{a numeric vector that represent the number of times that the string 'also' appeared}
    \item{\code{str.an}}{a numeric vector that represent the number of times that the string 'an' appeared}
    \item{\code{str.and}}{...}
    \item{\code{str.any}}{...}
    \item{\code{str.are}}{...}
    \item{\code{str.as}}{...}
    \item{\code{str.at}}{...}
    \item{\code{str.be}}{...}
    \item{\code{str.been}}{...}
    \item{\code{str.but}}{...}
    \item{\code{str.by}}{...}
    \item{\code{str.can}}{...}
    \item{\code{str.do}}{...}
    \item{\code{str.down}}{...}
    \item{\code{str.even}}{...}
    \item{\code{str.every}}{...}
    \item{\code{str.for}}{...}
    \item{\code{str.from}}{...}
    \item{\code{str.had}}{...}
    \item{\code{str.has}}{...}
    \item{\code{str.have}}{...}
    \item{\code{str.her}}{...}
    \item{\code{str.his}}{...}
    \item{\code{str.if}}{...}
    \item{\code{str.in}}{...}
    \item{\code{str.into}}{...}
    \item{\code{str.is}}{...}
    \item{\code{str.it}}{...}
    \item{\code{str.its}}{...}
    \item{\code{str.may}}{...}
    \item{\code{str.more}}{...}
    \item{\code{str.must}}{...}
    \item{\code{str.my}}{...}
    \item{\code{str.no}}{...}
    \item{\code{str.not}}{...}
    \item{\code{str.now}}{...}
    \item{\code{str.of}}{...}
    \item{\code{str.on}}{...}
    \item{\code{str.one}}{...}
    \item{\code{str.only}}{...}
    \item{\code{str.or}}{...}
    \item{\code{str.our}}{...}
    \item{\code{str.should}}{...}
    \item{\code{str.so}}{...}
    \item{\code{str.some}}{...}
    \item{\code{str.such}}{...}
    \item{\code{str.than}}{...}
    \item{\code{str.that}}{...}
    \item{\code{str.the}}{...}
    \item{\code{str.their}}{...}
    \item{\code{str.then}}{...}
    \item{\code{str.there}}{...}
    \item{\code{str.things}}{...}
    \item{\code{str.this}}{...}
    \item{\code{str.to}}{...}
    \item{\code{str.up}}{...}
    \item{\code{str.upon}}{...}
    \item{\code{str.was}}{...}
    \item{\code{str.were}}{...}
    \item{\code{str.what}}{...}
    \item{\code{str.when}}{...}
    \item{\code{str.which}}{...}
    \item{\code{str.who}}{...}
    \item{\code{str.will}}{...}
    \item{\code{str.with}}{...}
    \item{\code{str.would}}{...}
    \item{\code{str.your}}{...}
    \item{\code{BookID}}{a factor with levels \code{b1} \code{b2} \code{b3} \code{b4} \code{b5} \code{b6} \code{b7} \code{b8} \code{b9} \code{b10} \code{b11} \code{b12} which are corresponding to the book id}
    \item{\code{Author}}{a factor with levels \code{Austen} \code{London} \code{Milton} \code{Shakespeare} which are corresponding to the author id}
  }
}
\details{
This dataset is used to illustrate text classification.
}
\source{
	Jeffrey S. Simonoff
	\emph{Analyzing Categorical Data}
}
\references{
  ~~ possibly secondary sources and usages ~~
}
\examples{
authorship2 <- authorship

# transform data in order to have percentage
authorship2[, 1:69] <- authorship2[, 1:69] / rowSums(authorship2[, 1:69] )

# create the model 
authorship.som.init  <- som ( formula = ~ .
	, data = authorship2
	, neighborhood = "gaussian"
	, grid = grid ( xdim = 20 , ydim = 20 , type = "hexagonal" ) 
	)

# train the network
authorship.som <- learn( authorship.som.init , number.iter = 1000, max.alpha = 0.5, min.alpha = .001, max.rayon = 5 , step.eval.si = 100)

summary(authorship.som)

plot(authorship.som, "energy")

#--- see the distribution on the map
plot(authorship.som, "effectif", cex.label = 0)

}
\keyword{sysdata}
