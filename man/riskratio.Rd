\name{riskratio}
\alias{riskratio}
\title{Calculate risk ratio and its confidence intervals}
\description{
 Calculate risk ratio (a kind of relative risk) and its 
 confidence intervals based on approximation, followed by
 null hypothesis (risk ratio equals to 1) testing.
}
\usage{ riskratio(X, Y, m1, m2, conf.level=0.95) }
\arguments{
 \item{X}{The number of disease occurence among exposed cohort.}
 \item{Y}{The number of disease occurence among non-exposed cohort.}
 \item{m1}{The number of individuals in exposed cohort group.}
 \item{m2}{The number of individuals in non-exposed cohort group.}
 \item{conf.level}{Probability for confidence intervals.  Default is 0.95.}
}
\value{
 \item{estimate}{Calculated point estimate of risk ratio.}
 \item{conf.int}{A numeric vector of length 2 to give upper/lower limit of confidence intervals.}
 \item{p.value}{The significant probability of the result of null-hypothesis testing.}
}
\author{Minato Nakazawa \email{minato-nakazawa@umin.net} \url{http://minato.sip21c.org/}}
\examples{
 res <- riskratio(5, 10, 90, 90)
 str(res)
 print(res)
}
