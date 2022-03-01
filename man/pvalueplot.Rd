\name{pvalueplot}
\alias{pvalueplot}
\title{Drawing p-value function plot by a cross table}
\description{
Drawing the p-value function (a.k.a. nested confidence intervals) plot of risk ratio (RR) or odds ratio (OR) for a given 2 by 2 cross table, which is strongly recommended by Rothman KJ "Epidemiology: An introduction. 2nd Ed." Oxford Univ. Press.
 
Until fmsb-0.4.2, the formula to calculate p-values was not appropriate, so that the curve was not correct.  Through discussion with Professor Rothman, I realized my mistake, then fixed it in fmsb-0.4.3.  The feasible calculation is only possible in the manner of back-calculation from p-values to RR or OR, so that the calculation of p-values is restricted to the given range from 0.0005 to 1.
}
\usage{
 pvalueplot(XTAB, plot.OR, plot.log, xrange, add, ...)
}
\arguments{
 \item{XTAB}{A 2 by 2 matrix to draw p-value function (in another term, nested confidence intervals).  The table should be given as the cross table for the exposure status being column and the health outcome status being row, opposite from usual manner for cross tabulation.  To note, usually the numbers of incidence and the total observed numbers for exposed and nonexposed population as risk data, but in this function, the numbers of incidence and the remaining numbers without disease should be given as rows.}
 \item{plot.OR}{Logical.  If you want to draw the p-value function for the odds ratio, it should be set at TRUE, otherwise the p-value function for the risk ratio is drawn. Default FALSE.}
 \item{xrange}{A numeric vector includes 2 elements for minimum and maximum of x axis.  Default is c(0.01, 5).}
 \item{plot.log}{Logical. If TRUE, the horizontal axis becomes logarythmic scale.   Default FALSE.}
 \item{add}{Logical. If TRUE, the line is overlayed on the existing pvalueplot, otherwise the graph is newly plotted.  Default FALSE.}
 \item{...}{Other options handed down to plot() or lines().  pch, lty or col may be useful.}
}
\value{
 The data.frame containing the set of p-values (ranging from 0.0005 to 1) and corresponding RR or OR is returned.
}
\references{
 Rothman KJ (2012) Epidemiology: An Introduction. 2nd Ed., Oxford University Press, Oxford.
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\examples{
pvalueplot(matrix(c(321, 411, 686-321, 689-411), 2, 2), xrange=c(0.7, 0.9))
pvalueplot(matrix(c(4, 386, 4, 1250), 2, 2), xrange=c(0.1, 20), plot.log=TRUE)
pvalueplot(matrix(c(468, 480, 229, 205), 2, 2), plot.OR=TRUE, xrange=c(0.7, 1.0))
}
