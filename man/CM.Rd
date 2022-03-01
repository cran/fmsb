\name{CM}
\alias{CM}
\alias{fitCM}
\title{Coale and McNeil's model nupitiality and its fitting}
\description{
 Implementing Coale and McNeil's model (1972) for the age-specific probability of first marriage and fitting the model to actual data.
}
\usage{
 CM(scale=0.8, a0=18, k=2)
 fitCM(initialpar=c(0.8, 18, 2), data, ages=10:60, mode=1, Method="Nelder-Mead", ...)
}
\arguments{
 \item{scale}{The parameter C of the Coale-McNeil model, the scale parameter for total nupitiality.  Almost same as 1-(probability of never marring for whole life)}
 \item{a0}{The parameter a0 of the Coale-McNeil model, the age of the beginning of first marriage.  The beginning means about 1\% of the population having ever married.}
 \item{k}{The parameter k of the Coale-McNeil model, how fast the population marry after a0.}
 \item{initialpar}{Initial value for the parameters to be estimated.  If not given, c(0.8, 18, 2) is used.}
 \item{data}{Actual vector of the age-specific probability of first marriage when the mode is 1, otherwise the proportion ever married by each age.}
 \item{ages}{The age ranges for the data.  It must be within the range from 10 to 60.  Default is 10:60.  It must have the same lengths as data.}
 \item{mode}{The mode of fitting, 1 means fitting for the probability of first marriage, otherwise fitting for the proportion ever married.  Default is 1.}
 \item{Method}{The method to be used in optim() function.  Default is "Nelder-Mead".}
 \item{...}{Other options to be passed to optim().}}
\value{
 CM() returns model schedule of nupitiality for ages from 10 to 60 as a list, composed of g (the numeric vector for the probability of first marriage happening for each age), G (the numeric vector for the proportion ever married by each age), mu (mean age of first marriage among total population), and sigma (standard deviation of the ages of first marriage).
 fitCM() returns the numeric vector of fitted parameters C, a0 and k, RMSE for those values, and the flag of convergence.
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\seealso{\link[fmsb]{CT}}
\references{
Coale AJ, McNeil DR (1972) The distribution by age of the frequency of first marriage in a female cohort. Journal of the American Statistical Association, 67(340): 743-749.\doi{10.1080/01621459.1972.10481287}

Newell C (1988) Methods and Models in Demography. The Guilford Press, New York.
}
\examples{
# The data of Japanese population census 2010 for the whole country
# The proportion of ever married females for ages from 15 to 60.
# https://www.e-stat.go.jp/SG1/estat/List.do?bid=000001034991&cycode=0
 Ages <- 15:60
 EverMarriedFemales <-  c(0.003081039, 0.003203058, 0.006502558,
 0.014261608, 0.028378604, 0.048903318, 0.07596101, 0.110311095, 
 0.153365573, 0.2090648, 0.273819118, 0.342672073, 0.41259517, 
 0.479789489, 0.536291775, 0.589919881, 0.631937609, 0.663719195,
 0.691411757, 0.71775138, 0.740807817, 0.760155848, 0.775400475,
 0.788445244, 0.799522713, 0.81108241, 0.821591503, 0.830695486,
 0.840776283, 0.846773585, 0.85921777, 0.867991763, 0.876908992,
 0.886388747, 0.894302114, 0.902385961, 0.909329207, 0.914662575,
 0.920327092, 0.925013244, 0.929551158, 0.933150578, 0.935851652,
 0.938421122, 0.940089719, 0.943223398)

 res <- fitCM(initialpar=c(0.8, 18, 2), data=EverMarriedFemales, 
  ages=Ages, mode=2)
 print(res)
 plot(Ages, EverMarriedFemales, 
  main="Proportion ever married by each age\n for 2010 Japanese females")
 fitted <- CM(res[1], res[2], res[3])
 lines(Ages, fitted$G[6:51], col="red")
 NoteForm <- "C=\%3.1f, a0=\%3.1f, k=\%3.1f\n mu=\%3.1f, sd=\%3.1f"
 text(40, 0.2, sprintf(NoteForm, res[1], res[2], res[3], fitted$mu, fitted$sigma))
 # mean age of first marriage happening
 print(sum(Ages*fitted$g[Ages-9]/sum(fitted$g[Ages-9]))) 
}
