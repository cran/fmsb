\name{Hadwiger}
\alias{Hadwiger}
\alias{fitHad}
\title{Hadwiger's model fertility and its fitting}
\description{
 Implementing Hadwiger's model of age-specific fertility rates and fitting the model to actual ASFR.
}
\usage{
 Hadwiger(a=3.4, b=2.5, c=22.2)
 fitHad(initialpar=c(3.4, 2.5, 22.2), data, Method="Nelder-Mead", ...)
}
\arguments{
 \item{a}{The parameter a of the Hadwiger model, ASFR(x) = a*b/c*(c/x)^1.5*exp(-b^2*(c/x+x/c-2)) for age x from 15 to 54.}
 \item{b}{The parameter b of the Hadwiger model, ASFR(x) = a*b/c*(c/x)^1.5*exp(-b^2*(c/x+x/c-2)) for age x from 15 to 54.}
 \item{c}{The parameter c of the Hadwiger model, ASFR(x) = a*b/c*(c/x)^1.5*exp(-b^2*(c/x+x/c-2)) for age x from 15 to 54.}
 \item{initialpar}{Initial value for the parameters to be estimated.  If not given, c(3.4, 2.5, 22.2) is used.}
 \item{data}{Actual vector of ASFR (which must be given for from ages from 15 to 54 for each age) to be used to obtain the best-fit parameters of the Hadwiger's model.}
 \item{Method}{The method to be used in optim() function.  Default is "Nelder-Mead".}
 \item{...}{Other options to be passed to optim().}}
\value{
 Hadwiger() returns model ASFR for ages from 15 to 54.
 fitHad() returns the numeric vector of fitted parameters a, b and c, RMSE for those values, and the flag of convergence.
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\references{
Chandola T, Coleman DA, Horns RW (1999) Recent European fertility patterns:
fitting curves to 'distorted' distributions. Population Studies, 53(3): 317-329. \doi{10.1080/00324720308089}

}
\seealso{\link[fmsb]{Jfert}}
\examples{
 res <- fitHad(,Jfert$ASFR2000)
 FLAG <- res[5]
 while (FLAG>0) {
   res <- fitHad(res[1:3], Jfert$ASFR2000)
   FLAG <- res[5]
 }
 print(res)
}
