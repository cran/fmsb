\name{SIQR}
\alias{SIQR}
\title{Calculate semi-interquartile range}
\description{
  Calculate semi-interquartile range, using IQR or fivenum.
}
\usage{ SIQR(X, mode) }
\arguments{
 \item{X}{a numeric vector.}
 \item{mode}{If 1, using fivenum, otherwise using IQR function.  Default is 1.}
}
\value{
 A numeric vector of length 1, giving the semi-interquartile range.
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\examples{
 data <- rnorm(100, 10, 1)
 SIQR(data)
 SIQR(data, 2)
 sd(data)
 idata <- sample(50:80, 100, replace=TRUE)
 SIQR(idata)
 SIQR(idata, 2)
 sd(idata)
}
