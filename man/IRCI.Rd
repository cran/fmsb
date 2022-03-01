\name{IRCI}
\alias{IRCI}
\title{Calculate a incidence rate with confidence interval}
\description{
 Calculate a incidence rate with confidence interval.
}
\usage{ IRCI(a, PT, conf.level=0.9) }
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
 Rothman KJ (2012) Epidemiology: An Introduction. 2nd Ed., Oxford University Press, Oxford.
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\examples{
IRCI(8, 85000)
}
