\name{CronbachAlpha}
\alias{CronbachAlpha}
\title{Calculate Cronbach's alpha coefficient.}
\description{Calculate Cronbach's alpha coefficient from a matrix or data.frame with more than 2 columns.}
\usage{CronbachAlpha(X)}
\arguments{
 \item{X}{A matrix or data.frame with more than 2 columns.}
}
\value{
 Single numeric value of Cronbach's alpha.
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\references{Bland JM, Altman DG (1997) Statistics notes: Cronbach's alpha. BMJ, 314: 572.}
\examples{
 QUEST <- data.frame(
  Q1=c(1, 5, 2, 3, 4, 2, 3, 4, 3, 2), 
  Q2=c(2, 4, 1, 2, 4, 1, 2, 5, 2, 1), 
  Q3=c(2, 5, 1, 3, 3, 2, 2, 4, 2, 2))
 CronbachAlpha(QUEST)
}
