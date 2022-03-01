\name{Siler}
\alias{Siler}
\alias{fitSiler}
\title{Siler's model mortality for qx and its fitting}
\description{
 Implementing Siler's model mortality function of qx and fitting the model to actual qx of given lifetable.
}
\usage{
 Siler(a1, b1, a2, a3, b3, t)
 fitSiler(initialpar=c(0.01,3,1e-4,1e-5,0.1), data, mode=1, Method="Nelder-Mead", ...)
}
\arguments{
 \item{a1}{The parameter a1 of the Siler model, q(t)=a1*exp(-b1*t)+a2+a3*exp(b3*t).}
 \item{b1}{The parameter b1 of the Siler model, q(t)=a1*exp(-b1*t)+a2+a3*exp(b3*t).}
 \item{a2}{The parameter a2 of the Siler model, q(t)=a1*exp(-b1*t)+a2+a3*exp(b3*t).}
 \item{a3}{The parameter a3 of the Siler model, q(t)=a1*exp(-b1*t)+a2+a3*exp(b3*t).}
 \item{b3}{The parameter b3 of the Siler model, q(t)=a1*exp(-b1*t)+a2+a3*exp(b3*t).}
 \item{t}{Age (vector OK) in years}
 \item{initialpar}{Initial value for the parameters to be estimated.  If not given, c(0.01, 0.0003, 0.07) is used.}
 \item{data}{Actual vector of qx in the lifetable to be used to obtain the best-fit parameters of the Gompertz-Makeham model.}
 \item{mode}{Which of lifetable functions should be used to calculate the RMSE: 1 qx, 2 dx, otherwise lx.  Default is 1.}
 \item{Method}{The method to be used in optim() function.  Default is "Nelder-Mead".}
 \item{...}{Other options to be passed to optim().}}
\value{
 Siler() returns model qx for the same length with t.
 fitSiler() returns the numeric vector of fitted parameters a1, b1, a2, a3 and b3, RMSE for those values, and the flag of convergence.
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\seealso{\link[fmsb]{Jlife}}
\examples{
 res <- fitSiler(,Jlife$qx2005M)
 FLAG <- res[7]
 while (FLAG>0) {
   res <- fitSiler(res[1:5], Jlife$qx2005M)
   FLAG <- res[7]
 }
 print(res)
}
