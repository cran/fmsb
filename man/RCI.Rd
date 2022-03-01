\name{RCI}
\alias{RCI}
\title{Calculate risk and its confidence interval}
\description{
 Calculate risk and its confidence interval by the
 simple asymptotic method.
}
\usage{ RCI(a, N, conf.level=0.9) }
\arguments{
 \item{a}{Number of cases}
 \item{N}{Number of population at risk}
 \item{conf.level}{Probability for confidence intervals.  Default is 0.9.}
}
\value{
 \item{R}{Point estimate of risk.}
 \item{RL}{Lower limit of confidence interval}
 \item{RU}{Upper limit of confidence interval}
}
\references{
 Rothman KJ (2012) Epidemiology: An Introduction. 2nd Ed., Oxford University Press, Oxford.
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\examples{
# By simple asymptotic method
RCI(20, 100)
# By Wilson Score (without continuity correction)
prop.test(20, 100, conf.level=0.9, correct=FALSE)
# By Exact method
binom.test(20, 100, conf.level=0.9)
}
