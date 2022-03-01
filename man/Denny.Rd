\name{Denny}
\alias{Denny}
\alias{fitDenny}
\title{Denny's model mortality for lx and its fitting}
\description{
 Implementing Denny's model mortality function of lx and fitting the model to actual lx of given lifetable.
}
\usage{
 Denny(a, b, c, t)
 fitDenny(initialpar=rep(0.1, 3), data, mode=3, Method="Nelder-Mead", ...)
}
\arguments{
 \item{a}{The parameter a of the Denny model, l(t)=1/(1+a*(t/(105-t))^3+b*sqrt(exp(t/(105-t))-1)+c*(1-exp(-2*t))).}
 \item{b}{The parameter b of the Denny model, l(t)=1/(1+a*(t/(105-t))^3+b*sqrt(exp(t/(105-t))-1)+c*(1-exp(-2*t))).}
 \item{c}{The parameter c of the Denny model, l(t)=1/(1+a*(t/(105-t))^3+b*sqrt(exp(t/(105-t))-1)+c*(1-exp(-2*t))).}
 \item{t}{Age (vector OK) in years.  The t must be less than 105, otherwise the value by Denny() become 0.}
 \item{initialpar}{Initial value for the parameters to be estimated.  If not given, rep(0.1, 3) is used.}
 \item{data}{Actual vector of qx in the lifetable to be used to obtain the best-fit parameters of the Denny's model.  If the ages for qx are equal or elder than 105 years old, those will be ignored in fitting.}
 \item{mode}{Which of lifetable functions should be used to calculate the RMSE: 1 qx, 2 dx, otherwise lx.  Default is 3.}
 \item{Method}{The method to be used in optim() function.  Default is "Nelder-Mead".}
 \item{...}{Other options to be passed to optim().}}
\value{
 Denny() returns model lx for the same length with t.
 fitDenny() returns the numeric vector of fitted parameters a, b, and c, RMSE for those values, and the flag of convergence.
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\references{
 Denny C (1997) A model of the probability of survival from birth. Mathematical and Computer Modelling, 26: 69-78.
 \doi{10.1016/S0895-7177(97)00170-2}
}
\seealso{\link[fmsb]{Jlife}}
\examples{
 res <- fitDenny(,qxtolx(Jlife$qx2005M))
 FLAG <- res[5]
 while (FLAG>0) {
   res <- fitDenny(res[1:3], qxtolx(Jlife$qx2005M))
   FLAG <- res[5]
 }
 print(res)
}
