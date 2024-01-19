\name{GompertzMakeham}
\alias{GompertzMakeham}
\alias{fitGM}
\title{Gompertz-Makeham's model mortality for u(x) and its fitting}
\description{
 Implementing Gompertz-Makeham's model mortality function of the force of mortality u(x) with conversion into qx and fitting the model to actual qx of given lifetable.
}
\usage{
 GompertzMakeham(A, B, C, t)
 fitGM(initialpar=c(0.01, 0.0003, 0.07), data, mode=1, Method="Nelder-Mead", ...)
}
\arguments{
 \item{A}{The parameter A of the Gompertz-Makeham model u(t)=A*exp(B*t)+C.}
 \item{B}{The parameter B of the Gompertz-Makeham model u(t)=A*exp(B*t)+C.}
 \item{C}{The parameter C of the Gompertz-Makeham model u(t)=A*exp(B*t)+C.}
 \item{t}{Age (vector OK) in years}
 \item{initialpar}{Initial value for the parameters to be estimated.  If not given, c(0.01, 0.0003, 0.07) is used.}
 \item{data}{Actual vector of qx in the lifetable to be used to obtain the best-fit parameters of the Gompertz-Makeham model.}
 \item{mode}{Which of lifetable functions should be used to calculate the RMSE, which is to be minimized in optim() function: 1 qx, 2 dx, otherwise lx.  Default is 1.}
 \item{Method}{The method to be used in optim() function.  Default is "Nelder-Mead".}
 \item{...}{Other options to be passed to optim().}
}
\value{
 GompertzMakeham() returns model qx for the same length with t, where u(x) is internally converted into qx.
 fitGM() returns the numeric vector of fitted parameters of A, B and C, RMSE for those values, and the flag of convergence.
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\seealso{\link[fmsb]{Jlife}}
\examples{
 res <- fitGM(,Jlife$qx2005M)
 FLAG <- res[5]
 while (FLAG>0) {
   res <- fitGM(res[1:3], Jlife$qx2005M)
   FLAG <- res[5]
 }
 print(res)
}
