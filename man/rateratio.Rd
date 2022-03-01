\name{rateratio}
\alias{rateratio}
\title{Calculate incidence rate ratio and its confidence intervals}
\description{
 Calculate incidence rate ratio (a kind of relative risk) and its 
 confidence intervals based on approximation, followed by
 null hypothesis (incidence rate ratio equals to 1) testing.
}
\usage{ rateratio(a, b, PT1, PT0, conf.level=0.95) }
\arguments{
 \item{a}{The number of disease occurence among exposed cohort.}
 \item{b}{The number of disease occurence among non-exposed cohort.}
 \item{PT1}{The observed person-time of the exposed cohort.}
 \item{PT0}{The observed person-time of the unexposed cohort.}
 \item{conf.level}{Probability for confidence intervals.  Default is 0.95.}
}
\value{
 \item{estimate}{Calculated point estimate of incidence rate ratio.}
 \item{conf.int}{A numeric vector of length 2 to give upper/lower limit of confidence intervals.}
 \item{p.value}{The significant probability of the result of null-hypothesis testing.}
}
\references{
 Rothman KJ (2012) Epidemiology: An Introduction. 2nd Ed., Oxford University Press, Oxford.
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\examples{
 res <- rateratio(136, 1709, 22050, 127650)
 str(res)
 print(res)
}
