\name{pvpORMH}
\alias{pvpORMH}
\title{Drawing p-value function plot by stratified or pooled cross tables}
\description{
Drawing the p-value function (a.k.a. nested confidence intervals) plot of pooled odds ratios (pORs) for several 2 by 2 crosstables, which are stratified by a confounding variable or pooled for several studies, with Mantel-Haenszel's method.
}
\usage{
 pvpORMH(XTAB, xrange, add, ...)
}
\arguments{
 \item{XTAB}{A matrix with 4 columns.
 The first column is the number of exposed cases.
 The second column is the number of unexposed cases.
 The third column is the number of exposed controls.
 The forth column is the number of unexposed controls.
 Rows should be composed of different strata or studies.}
 \item{xrange}{A numeric vector includes 2 elements for minimum and maximum of x axis.  Default is c(0.6, 1.2).}
 \item{add}{Logical. If TRUE, the line is overlayed on the existing pvalueplot, otherwise the graph is newly plotted.  Default FALSE.}
 \item{...}{Other options handed down to plot() or lines().  pch, lty or col may be useful.}
}
\value{
 A data.frame containing the set of p-values (ranging from 0.0005 to 1) and corresponding pORs are returned.
}
\references{
 Rothman KJ (2012) Epidemiology: An Introduction. 2nd Ed., Oxford University Press, Oxford.
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\examples{
TenStudies <- matrix(
 c(215, 229, 311-215, 306-229,
   38, 33, 59-38, 51-33,
   161, 174, 293-161, 293-174,
   76, 88, 164-76, 163-88,
   103, 105, 129-103, 133-105,
   65, 67, 120-65, 125-67,
   81, 75, 113-81, 110-75,
   48, 63, 160-48, 159-63,
   22, 21, 60-22, 62-21,
   56, 51, 137-56, 140-51
   ), 10, 4, byrow=TRUE)
ElevenStudies <- rbind(TenStudies, c(468, 480, 229, 205))
# Figure 8-4 in Chapter 8 of Rothman's textbook.
pvpORMH(TenStudies)
pvpORMH(ElevenStudies, add=TRUE, lty=2)
segments(1, 0, 1, 1, lwd=2)
}
