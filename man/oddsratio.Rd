\name{oddsratio}
\alias{oddsratio}
\title{Calculate odds ratio and its confidence intervals}
\description{
 Calculate odds ratio and its confidence intervals 
 based on approximation, followed by null-hypothesis
 (odds ratio equals to 1) testing.
}
\usage{ oddsratio(a, b, c, d, conf.level=0.95) }
\arguments{
 \item{a}{The number of individuals who both suffer from exposure and disease.}
 \item{b}{The number of individuals who suffer from disesase but not exposed.}
 \item{c}{The number of individuals who suffer from exposure but are healthy.}
 \item{d}{The number of individuals who neither suffered from exposure nor disease.}
 \item{conf.level}{Probability for confidence intervals.  Default is 0.95.}
}
\value{
 \item{estimate}{Calculated point estimate of odds ratio.}
 \item{conf.int}{A numeric vector of length 2 to give upper/lower limit of confidence intervals.}
 \item{p.value}{The significant probability as the result of null-hypothesis testing.}
}
\author{Minato Nakazawa \email{minato-nakazawa@umin.net} \url{http://minato.sip21c.org/}}
\examples{
 res <- oddsratio(5,10,85,80)
 str(res)
 print(res)
}
