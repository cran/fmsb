\name{riskratio}
\alias{riskratio}
\title{Calculate risk ratio and its confidence intervals}
\description{
 Calculate risk ratio (a kind of relative risk) and its 
 confidence intervals based on approximation, followed by
 null hypothesis (risk ratio equals to 1) testing.
}
\usage{ riskratio(X, Y, m1, m2, conf.level=0.95, p.calc.by.independence=TRUE) }
\arguments{
 \item{X}{The number of disease occurence among exposed cohort.}
 \item{Y}{The number of disease occurence among non-exposed cohort.}
 \item{m1}{The number of individuals in exposed cohort group.}
 \item{m2}{The number of individuals in non-exposed cohort group.}
 \item{conf.level}{Probability for confidence intervals.  Default is 0.95.}
 \item{p.calc.by.independence}{Logical.  If TRUE, calculating p-value by 
  testing the null-hypothesis of independence between exposure and disease.
  Otherwise, calculating p-value by inverse-function of confidence intervals
  calculation (the result becomes the same as the vcd package). Default TRUE.}
}
\value{
 \item{estimate}{Calculated point estimate of risk ratio.}
 \item{conf.int}{A numeric vector of length 2 to give upper/lower limit of confidence intervals.}
 \item{p.value}{The significant probability of the result of null-hypothesis testing.}
}
\references{
 Rothman KJ (2012) Epidemiology: An Introduction. 2nd Ed., Oxford University Press, Oxford.
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\examples{
 res <- riskratio(5, 10, 90, 90)
 str(res)
 print(res)
 riskratio(12, 5, 18, 17)
 riskratio(12, 5, 18, 17, p.calc.by.independence=FALSE)
}
