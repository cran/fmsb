\name{pvalueplot}
\alias{pvalueplot}
\title{Drawing p-value function plot by a cross table}
\description{
 Drawing the p-value function plot using a 2 by 2 cross table, which is strongly recommended by Rothman KJ "Epidemiology: An introduction." Oxford Univ. Press.
}
\usage{
 pvalueplot(XTAB, plot.OR, xrange, scale, plot.log)
}
\arguments{
 \item{XTAB}{A 2 by 2 matrix to draw p-value function.}
 \item{plot.OR}{Logical.  If draw odds ratio, it's TRUE, otherwise FALSE. Default FALSE.}
 \item{xrange}{A numeric vector includes 2 elements for minimum and maximum of x axis.  Default is c(0.01, 5).}
 \item{scale}{A numeric value indicates the drawing density.  If 10, p-values are calculated for all possible XTAB counts when the peripheral counts were fixed and 9 in-between values for each possible counts (I'm sorry for such bad explanation).  Default is 1.}
 \item{plot.log}{Logical.  Default FALSE.}
}
\value{
 No value is returned.
}
\author{Minato Nakazawa \email{minato-nakazawa@umin.net} \url{http://minato.sip21c.org/}}
\examples{
pvalueplot(matrix(c(321,411,686-321,689-411),2,2),xrange=c(0.7,0.9))
pvalueplot(matrix(c(4,386,4,1250),2,2),xrange=c(0.1,20),scale=10,plot.log=TRUE)
pvalueplot(matrix(c(468,480,229,205),2,2),plot.OR=TRUE,xrange=c(0.7,1.0))
}
