\name{RRMH}
\alias{RRMH}
\title{Calculate pooled risk ratio and its confidence intervals
 with Mantel-Haenszel's method}
\description{
 Calculate pooled risk ratio and its confidence intervals 
 with Mantel-Haenszel's method.
}
\usage{ RRMH(XTAB, conf.level=0.9) }
\arguments{
 \item{XTAB}{A matrix with 4 columns.
 The first column is the number of disease occurrence in exposed cohort.
 The second column is the number of disease occurrence in unexposed cohort.
 The third column is the total number of exposed cohort.
 The forth column is the total number of unexposed cohort.
 Rows should be composed of different strata or studies.}
 \item{conf.level}{Probability for confidence intervals.  Default is 0.9.}
}
\value{
 \item{estimate}{Calculated point estimate of pooled risk ratio with Manterl-Haenszel's method.}
 \item{conf.int}{A numeric vector of length 2 to give upper/lower limit of confidence intervals.}
 \item{conf.level}{Simply return the value of given conf.level.}
}
\references{
 Rothman KJ (2012) Epidemiology: An Introduction. 2nd Ed., Oxford University Press, Oxford.
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\examples{
# Table 10-3 of Rothman's textbook (Chapter 10).
RRMH(matrix(c(8, 5, 106, 120, 22, 16, 98, 85), 2, 4, byrow=TRUE), conf.level=0.9)
}
