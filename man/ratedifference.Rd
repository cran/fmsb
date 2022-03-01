\name{ratedifference}
\alias{ratedifference}
\title{Calculate incidence rate difference and its confidence intervals}
\description{
 Calculate incidence rate difference (a kind of attributable risk / excess 
 risk) and its confidence intervals based on approximation, followed by 
 null hypothesis (incidence rate difference equals to 0) testing.
}
\usage{ ratedifference(a, b, PT1, PT0, CRC=FALSE, conf.level=0.95) }
\arguments{
 \item{a}{The number of disease occurence among exposed cohort.}
 \item{b}{The number of disease occurence among non-exposed cohort.}
 \item{PT1}{The observed person-time of the exposed cohort.}
 \item{PT0}{The observed person-time of the unexposed cohort.}
 \item{CRC}{Logical.  If TRUE, calculate confidence intervals for each incidence rate.  Default is FALSE.}
 \item{conf.level}{Probability for confidence intervals.  Default is 0.95.}
}
\value{
 \item{estimate}{Calculated point estimate of incidence rate difference.}
 \item{conf.int}{A numeric vector of length 2 to give upper/lower limit of confidence intervals.}
 \item{p.value}{The significant probability of the result of null-hypothesis testing.}
}
\references{
 Rothman KJ (2012) Epidemiology: An Introduction. 2nd Ed., Oxford University Press, Oxford.
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\examples{
 res <- ratedifference(136, 1709, 22050, 127650, CRC=TRUE)
 str(res)
 print(res)
}
