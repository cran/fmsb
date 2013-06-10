\name{roc}
\alias{roc}
\alias{print.roc}
\alias{plot.roc}
\title{Calculate Receiver Operating Characteristic (ROC) curve}
\description{
 Calculate Receiver Operating Characteristic (ROC) curve's each performance set of [sensitivity, 1-specificity], each distance of the performance from the worst performance [0, 1], and each piece of area under the curve, for each cutoff point, as list.
}
\usage{ roc(values, iscase) }
\arguments{
 \item{values}{A numeric vector of measured values.}
 \item{iscase}{A logical (or 0/1) vector of diagnostics.  Negative result must be given by FALSE or 0.}
}
\value{
 \item{cutoff}{The numeric vector of cutoff points, which are composed of the all unique values among the given measurements and the maximum cutoff is maximum measurement plus 1.  Therefore, the minimum cutoff gives [1, 1] and the maximum cutoff gives [0, 0] as the performance set of [sensitivity, 1-specificity], respectively.}
 \item{sens}{The numeric vector of sensitivities for all cutoff points.}
 \item{falsepos}{The numeric vector of 1-specificities (false positiveness) for all cutoff points.}
 \item{distance}{The numeric vector of distance between actual performance set and the worst performance.}
 \item{aucpiece}{The numeric vector of the pieces of areas under the curve.}
}
\author{Minato Nakazawa \email{minato-nakazawa@umin.net} \url{http://minato.sip21c.org/}}
\seealso{\link[Epi]{ROC}}
\examples{
 scores <- c(15, 20, 19, 28, 26, 17, 13, 22, 23, 24)
 diagno <- c(0, 0, 0, 1, 1, 1, 0, 1, 1, 1)
 res <- roc(scores, diagno)
 print(res)
 plot(res)
}
