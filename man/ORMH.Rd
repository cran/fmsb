\name{ORMH}
\alias{ORMH}
\title{Calculate pooled odds ratio and its confidence intervals
 with Mantel-Haenszel's method}
\description{
 Calculate pooled odds ratio and its confidence intervals 
 with Mantel-Haenszel's method.
}
\usage{ ORMH(TBL, conf.level=0.95) }
\arguments{
 \item{TBL}{A matrix with 4 columns.
 The first column is the number of exposed cases.
 The second column is the number of unexposed cases.
 The third column is the number of exposed controls.
 The forth column is the number of unexposed controls.
 Rows should be composed of different strata or studies.}
 \item{conf.level}{Probability for confidence intervals.  Default is 0.95.}
}
\value{
 \item{estimate}{Calculated point estimate of pooled odds ratio with Manterl-Haenszel's method.}
 \item{conf.int}{A numeric vector of length 2 to give upper/lower limit of confidence intervals.}
 \item{conf.level}{Simply return the value of given conf.level.}
}
\references{
 Rothman KJ (2012) Epidemiology: An Introduction. 2nd Ed., Oxford University Press, Oxford.
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\examples{
# Table 10-6 of Rothman's textbook (Chapter 10).
ORMH(matrix(c(3, 9, 104, 1059, 1, 3, 5, 86), 2, 4, byrow=TRUE), conf.level=0.9)
# Figure 8-4 of Rothman's textbook (Chapter 8)
# https://www.ncbi.nlm.nih.gov/pubmed/7630245
# https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(05)74403-2/fulltext
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
ORMH(TenStudies)
ElevenStudies <- rbind(TenStudies, c(468, 480, 229, 205))
ORMH(ElevenStudies)
}
