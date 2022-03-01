\name{riskdifference}
\alias{riskdifference}
\title{Calculate risk difference and its confidence intervals}
\description{
 Calculate risk difference (a kind of attributable risk / excess 
 risk) and its confidence intervals based on approximation, followed by 
 null hypothesis (risk difference equals to 0) testing.
}
\usage{ riskdifference(a, b, N1, N0, CRC=FALSE, conf.level=0.95) }
\arguments{
 \item{a}{The number of disease occurence among exposed cohort.}
 \item{b}{The number of disease occurence among non-exposed cohort.}
 \item{N1}{The population at risk of the exposed cohort.}
 \item{N0}{The population at risk of the unexposed cohort.}
 \item{CRC}{Logical.  If TRUE, calculate confidence intervals for each risk.  Default is FALSE.}
 \item{conf.level}{Probability for confidence intervals.  Default is 0.95.}
}
\value{
 \item{estimate}{Calculated point estimate of risk difference.}
 \item{conf.int}{A numeric vector of length 2 to give upper/lower limit of confidence intervals.}
 \item{p.value}{The significant probability of the result of null-hypothesis testing.}
}
\references{
 Rothman KJ (2012) Epidemiology: An Introduction. 2nd Ed., Oxford University Press, Oxford.
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\examples{
 res <- riskdifference(321, 411, 686, 689, CRC=TRUE)
 str(res)
 print(res)
}
