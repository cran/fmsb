\name{IRCIPois}
\alias{IRCIPois}
\title{Calculate incidence rate with its confidence
 intervals by exact method}
\description{
Calculate incidence rate with its confidence
 intervals by exact method using Poisson distribution.
}
\usage{ IRCIPois(a, PT, conf.level=0.9) }
\arguments{
 \item{a}{Number of cases}
 \item{PT}{Person-years of observed population at risk}
 \item{conf.level}{Probability for confidence intervals.  Default is 0.9.}
}
\value{
 \item{IR}{Point estimate of incidence rate.}
 \item{IRL}{Lower limit of confidence interval}
 \item{IRU}{Upper limit of confidence interval}
}
\references{
 \url{https://www.statsdirect.com/help/rates/poisson_rate_ci.htm}
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\examples{
IRCIPois(8, 85000)
}
