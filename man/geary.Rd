\name{geary.test}
\alias{geary.test}
\title{Geary's test for normality about kurtosis.}
\description{
 Geary's test for normality.  Null hypothesis is that the data obeys to normal distribution.
}
\usage{ geary.test(X) }
\arguments{
 \item{X}{A numeric vector.}
}
\value{
 \item{statistic}{Geary's test statistic G}
 \item{p.value}{The significant probability of the null-hypothesis testing.}
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\examples{
 geary.test(rnorm(100))
 geary.test(20:50)
}
