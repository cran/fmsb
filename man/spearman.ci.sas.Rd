\name{spearman.ci.sas}
\alias{spearman.ci.sas}
\title{Calculate Spearman's rank correlation with its confidence intervals by SAS method}
\description{
 Calculate Spearman's rank correlation with its 
 confidence intervals by the same method as SAS.
}
\usage{ spearman.ci.sas(x, y, adj.bias=TRUE, conf.level=0.95) }
\arguments{
 \item{x}{A numeric vector.}
 \item{y}{A numeric vector.}
 \item{adj.bias}{Logical.  If TRUE, adjustment for bias is taken.  Default TRUE.}
 \item{conf.level}{Probability for confidence intervals.  Default is 0.95.}
}
\value{
 \item{rho}{Calculated point estimate of Spearman's rank correlation coefficient.}
 \item{rho.ll}{The lower limit of given confidence intervals.}
 \item{rho.ul}{The upper limit of given confidence intervals.}
 \item{adj.bias}{The option for bias adjustment taken.}
}
\author{Minato Nakazawa \email{minato-nakazawa@umin.net} \url{http://minato.sip21c.org/}}
\references{
 \url{http://support.sas.com/documentation/cdl/en/procstat/63104/HTML/default/corr_toc.htm}
}
\examples{
 data(airquality)
 spearman.ci.sas(airquality$Ozone, airquality$Wind)
}
