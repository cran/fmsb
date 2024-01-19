\name{IRDMH}
\alias{IRDMH}
\title{Calculate pooled incidence rate difference and its confidence intervals
 with Mantel-Haenszel's method}
\description{
 Calculate pooled incidence rate difference and its confidence intervals 
 with Mantel-Haenszel's method.
}
\usage{ IRDMH(XTAB, conf.level=0.9) }
\arguments{
 \item{XTAB}{A matrix with 4 columns.
 The first column is the incidence in the exposed cohort.
 The second column is the incidence in the unexposed cohort.
 The third column is the observed person-time of exposed cohort.
 The forth column is the observed person-time of unexposed cohort.
 Rows should be composed of different strata or studies.}
 \item{conf.level}{Probability for confidence intervals.  Default is 0.9.}
}
\value{
 \item{estimate}{Calculated point estimate of pooled incidence rate difference with Manterl-Haenszel's method.}
 \item{conf.int}{A numeric vector of length 2 to give upper/lower limit of confidence intervals.}
 \item{conf.level}{Simply return the value of given conf.level.}
}
\references{
 Rothman KJ (2012) Epidemiology: An Introduction. 2nd Ed., Oxford University Press, Oxford.
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\examples{
# Table 10-5 of Rothman's textbook (Chapter 10).
IRDMH(matrix(c(196, 111, 62119, 15763, 167, 157, 6085, 2780), 2, byrow=TRUE), conf.level=0.9)
}
