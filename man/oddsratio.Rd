\name{oddsratio}
\alias{oddsratio}
\title{Calculate odds ratio and its confidence intervals}
\description{
 Calculate odds ratio and its confidence intervals 
 based on approximation, followed by null-hypothesis
 (odds ratio equals to 1) testing.
}
\usage{ oddsratio(a, b, c, d, conf.level=0.95, p.calc.by.independence=TRUE) }
\arguments{
 \item{a}{The number of individuals who both suffer from exposure and disease.}
 \item{b}{The number of individuals who suffer from disesase but not exposed.}
 \item{c}{The number of individuals who suffer from exposure but are healthy.}
 \item{d}{The number of individuals who neither suffered from exposure nor disease.}
 \item{conf.level}{Probability for confidence intervals.  Default is 0.95.}
 \item{p.calc.by.independence}{Logical.  If TRUE, calculating p-value by 
  testing the null-hypothesis of independence between exposure and disease.
  Otherwise, calculating p-value by inverse-function of confidence intervals
  calculation (the result becomes the same as the vcd package). Default TRUE.}
}
\value{
 \item{estimate}{Calculated point estimate of odds ratio.}
 \item{conf.int}{A numeric vector of length 2 to give upper/lower limit of confidence intervals.}
 \item{p.value}{The significant probability as the result of null-hypothesis testing.}
}
\author{Minato Nakazawa \email{minato-nakazawa@umin.net} \url{http://minato.sip21c.org/}}
\examples{
 res <- oddsratio(5, 10, 85, 80)
 str(res)
 print(res)
 oddsratio(12, 5, 6, 12)
 oddsratio(12, 5, 6, 12, p.calc.by.independence=FALSE)
}
