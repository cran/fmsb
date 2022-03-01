\name{mhchart}
\alias{mhchart}
\title{Draw maternity history charts.}
\description{To compare the maternity histories among several human populations, this kind of graph is useful, inspired by Wood JW (1994) "Dynamics of Human Reproduction", Aldine de Gruyter, New York.}
\usage{mhchart(LIST, XLIM=c(15,45), COL="black", FILL="white", BWD=1, ...)}
\arguments{
 \item{LIST}{The list of groups with their maternity histories from first birth to the last birth.  Usually the first childbirth age is estimated as median by Kaplan-Meier method, the second childbirth age was given by adding the median of first birth intervals to the first childbirth age by Kaplan-Meier method, and so on.}
 \item{XLIM}{The limit of x axis, which means the range of reproductive ages.  Default is 15 and 45.}
 \item{COL}{The border color.  Default is black.}
 \item{FILL}{The painting color.  Default is white.}
 \item{BWD}{The line width of the boxes.  Default is 1.}
 \item{...}{Other parameters handed to barplot() to draw axes and background.}
}
\value{
 No value is returned.
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\examples{
 Developing <- c(18, 21, 24, 27, 30, 33.5, 37)
 Hutterite <- c(23, 25, 27, 29, 31, 33, 35, 37, 39)
 Gainj <- c(27, 31, 35, 39)
 Japan <- c(29, 34)
 x <- list(
  Developing=Developing,
  Hutterite=Hutterite,
  Gainj=Gainj,
  Japan=Japan)
 mhchart(rev(x), COL="blue", FILL="pink", BWD=2, XLIM=c(15, 45),
  main="Maternity histories for selected populations",
  xlab="Maternal age (years)")
}
