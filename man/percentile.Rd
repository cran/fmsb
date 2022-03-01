\name{percentile}
\alias{percentile}
\title{Convert numeric vector into its percentile}
\description{
 Convert numeric vector into its percentile.  For example, 1:5 will become c(0,25,50,75,100).
}
\usage{ percentile(dat) }
\arguments{
 \item{dat}{A numeric vector, which will be converted into percentile value.}
}
\value{
 A integer vector in [0,100].  Minimum value always becomes 0 and maximum always becomes 100.
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\examples{
 percentile(1:5)
 X <- runif(1000, 10, 20)
 percentile(X)
}
