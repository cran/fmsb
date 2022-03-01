\name{CT}
\alias{CT}
\alias{fitCT}
\title{Coale and Trussell's model marital fertility and its fitting}
\description{
 Implementing Coale and Trussell's model of age-specific marital fertility rates and fitting the model to actual ASMFR.
}
\usage{
 CT(M=1, m=0)
 fitCT(initialpar=c(1.0, 1.0), data, Method="Nelder-Mead", ...)
}
\arguments{
 \item{M}{The parameter M of the CT model, the scale (peak height) parameter of fertility}
 \item{m}{The parameter m of the CT model, the strength of downward discordance from natural fertility with aging}
 \item{initialpar}{Initial value for the parameters to be estimated.  If not given, c(1.0, 1.0) is used.}
 \item{data}{Actual vector of ASMFR (which must be given for from age 12 to age 49 for each age) to be used to obtain the best-fit parameters of the CT's model.}
 \item{Method}{The method to be used in optim() function.  Default is "Nelder-Mead".}
 \item{...}{Other options to be passed to optim().}}
\value{
 CT() returns model ASMFR for ages from 12 to 49.
 fitCT() returns the numeric vector of fitted parameters M and m, RMSE for those values, and the flag of convergence.
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\seealso{\link[fmsb]{Jfert}}
\references{
Coale AJ, Trussell TJ (1978) Technical Note: Finding the Two Parameters That Specify a Model Schedule of Marital Fertility.  Population Index, 44(2): 203-213.

}
\examples{
 ASMFR <- c(0, 0, 0, Jfert$ASMFR2000[1:35]) # Jfert's ASMFR should be rearranged to 12:49
 res <- fitCT(,ASMFR)
 FLAG <- res[4]
 while (FLAG>0) {
   res <- fitCT(res[1:2], ASMFR)
   FLAG <- res[4]
 }
 print(res)
}
