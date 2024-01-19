\name{truemedian}
\alias{truemedian}
\title{Calculate true median for data with ties.}
\description{
  Usually median for data with ties, tied values are treated as
  exactly same.  For example, median of \{3, 3, 4, 4, 4\} will be 4.
  However, the measured values are usually rounded, so that we can
  assume evenly distributed true values for tied values.  For example,
  the previous data can be treated as rounded values of \{2.75, 3.25,
  11/3, 4, 13/3\}.  From this viewpoint, true median of \{3, 3, 4, 4, 4\}
  could be 11/3 (=3.66...).  This function calculates this.
}
\usage{ truemedian(X, h) }
\arguments{
 \item{X}{A numeric vector.  Usually integer.}
 \item{h}{Width of measurement unit.  Default is 1.}
}
\value{
 A numeric vector of length 1, giving "true" median estimate.
}
\references{
 Grimm LG (1993) Statistical Applications for the Behavioral Sciences. John Wiley and Sons.
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\examples{
 median(c(3, 3, 4, 4, 4))
 truemedian(c(3, 3, 4, 4, 4))
}
