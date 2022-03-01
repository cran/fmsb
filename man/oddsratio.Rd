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
 \item{a}{A scalar or a matrix.  If matrix, it has to be 2 by 2, which contains
 the number of individuals who both suffer from exposure and disease as [1, 1],
 the number of individuals who suffer from disesase but not exposed as [2, 1],
 the number of individuals who suffer from exposure but are healthy as [1, 2] and
 the number of individuals who neither suffered from exposure nor disease as [2, 2].}
 \item{b}{If a is a scalar, this has to be given as the number of individuals who 
 suffer from disesase but not exposed.  Otherwise, ignored.}
 \item{c}{If a is a scalar, this has to be given as the number of individuals who 
 suffer from exposure but are healthy.  Otherwise, ignored.}
 \item{d}{If a is a scalar, this has to be given as the number of individuals who 
 neither suffered from exposure nor disease.  Otherwise, ignored.}
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
\note{This function can also accept a matrix as argument, as suggested by Dr. Toshiaki Ara
 (\email{toshiaki.ara@gmail.com}).  Thanks for a good suggestion.}
\references{
 Rothman KJ (2012) Epidemiology: An Introduction. 2nd Ed., Oxford University Press, Oxford.
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\examples{
 res <- oddsratio(5, 10, 85, 80)
 str(res)
 print(res)
 oddsratio(12, 5, 6, 12)
 oddsratio(12, 5, 6, 12, p.calc.by.independence=FALSE)
 DH <- sample(c("Disease", "Health"), 100, replace=TRUE)
 EN <- sample(c("Exposed", "Nonexposed"), 100, replace=TRUE)
 x <- table(EN, DH)
 oddsratio(x)
 # same as oddsratio(x[1,1], x[2,1], x[1,2], x[2,2])
}
