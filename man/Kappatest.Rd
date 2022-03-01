\name{Kappa.test}
\alias{Kappa.test}
\title{Calculate Cohen's kappa statistics for agreement}
\description{
 Calculate Cohen's kappa statistics for agreement and its confidence intervals 
 followed by testing null-hypothesis that the extent of agreement is same as random,
 kappa statistic equals zero.
}
\usage{ Kappa.test(x, y=NULL, conf.level=0.95) }
\arguments{
 \item{x}{If y is not given, x must be the square matrix that the rows and columns 
 show the ratings of different rater (or repeated measure) and the values indicate 
 the numbers of data having that combination.  If y is given, x must be the result 
 of ratings by the first rater (or first time measurement).}
 \item{y}{If given, y must be the result of ratings by the second rater (or second 
 time measurement).  As default, it is not given.}
 \item{conf.level}{Probability for confidence intervals for kappa statistics.  Default is 0.95.}
}
\value{
 \item{Result$statistic}{Z score to test null-hypothesis.}
 \item{Result$estimate}{Calculated point estimate of Cohen's kappa statistic.}
 \item{Result$conf.int}{A numeric vector of length 2 to give upper/lower limit of 
  confidence intervals.}
 \item{Result$p.value}{The significant probability as the result of null-hypothesis testing.}
 \item{Judgement}{The judgement for the estimated kappa about the extent of agreement, 
  given by Landis JR, Koch GG (1977) Biometrics, 33: 159-174: If kappa is less than 0, "No 
  agreement", if 0-0.2, "Slignt agreement", if 0.2-0.4, "Fair agreement", if 0.4-0.6, 
  "Moderate agreement", if 0.6-0.8, "Substantial agreement", if 0.8-1.0, 
  "Almost perfect agreement".}
}
\references{
 Landis JR, Koch GG (1977) The measurement of observer agreement for categorical 
  data.  Biometrics, 33: 159-174.
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\seealso{\link[vcd]{Kappa}}
\examples{
 res <- Kappa.test(matrix(c(20, 10, 5, 15), 2, 2))
 str(res)
 print(res)
 Kappa.test(c(1, 1, 3, 1, 1, 2, 1, 2, 1, 1), c(2, 1, 3, 1, 3, 2, 1, 3, 3, 3))
}
