\name{lifetable}
\alias{lifetable}
\alias{lifetable2}
\alias{lifetable3}
\alias{clifetable}
\alias{lxtodx}
\alias{dxtolx}
\alias{qxtodx}
\alias{dxtoqx}
\alias{qxtomx}
\alias{mxtoqx}
\alias{qxtolx}
\alias{lxtoqx}
\alias{uxtoqx}
\alias{getax}
\alias{hlifetable}
\alias{Jp2010s}
\title{Lifetable functions with mutual conversions}
\description{
 Lifetable related functions.
}
\usage{
 lifetable(mx, ns=NULL, class=5, mode=1)
 lifetable2(mx, ax=0.5, n=1)
 lifetable3(lx, ax=0.5, n=1)
 clifetable(qx)
 lxtodx(lx)
 dxtolx(dx)
 qxtodx(qx)
 dxtoqx(dx)
 qxtomx(qx, ax=0.5, n=1, mmax=NULL)
 mxtoqx(mx, ax=0.5, n=1)
 qxtolx(qx)
 lxtoqx(lx)
 uxtoqx(ux)
 hlifetable(mx, ax=0.5, n=5, pix=0, Nx=NULL, conf.level=0.95)
 getax(lx, Tx, n=5)
}
\arguments{
 \item{mx}{Lifetable function mx, meaning the age (class)-specific death rates.}
 \item{ns}{If given as a vector with the same length as mx, the duration for each age (class). Default is NULL: same duration with class is automatically used.}
 \item{n}{If given as a vector with the same length as mx or qx, the duration for each age (class).  Default is 1, which means the length for all age-classes being 1 year.}
 \item{class}{Age-class of lifetable() function.  Default is 5.}
 \item{mode}{How to set ax and correction method in conversion from mx to qx.  1 and 11: all ax is 0.5 except the open-ended class [where ax is reciprocal of mx], 2, 4, 12, 14: ax is 0.1 for age 0, 0.4 for age 1-4, 0.5 for the other ages except the open-ended class [where ax is reciprocal of mx], 3, 5, 13, 15: ax is 0.3 for age 0, 0.4 for age 1-4, 0.5 for the other ages except the open-ended class [where ax is reciprocal of mx], 6 and 16: Males value given in Preston SH (2001), pp.48 Table 3.3, 7 and 17: Females value given in Preston SH (2001), pp.48 Table 3.3.  If less than 10, simply calculating qx as n*mx/(1+n*(1-ax)*mx) (Note: In the formula of Preston SH (2001) pp.47, the function is given as n*mx/(1+(n-ax)*mx).  The difference is due to the formulation of ax.  In this function, ax is given for single age, same as Newell C (1988), pp.71) except for the open-ended class where qx=1, otherwise calculating qx by Greville's method.  Default is 1.}
 \item{ax}{Lifetable function ax, fraction of last year lived.  Default is 0.5 (scalar) for all classes.  It can be given as scalar or vector.  Note: This argument can only be specified in lifetable2() or lifetable3(), not in lifetable().}
 \item{qx}{Lifetable function qx, which means the probability of dying between age x and x+1 (for lifetable(), x+class).}
 \item{mmax}{To calculate mx from qx, mx at the maximum open-ended age-class cannot be calculated from qx.  In such situation, mmax gives a value for it.  Default is NULL.}
 \item{lx}{Lifetable function lx, which means number of people left alive at age x from 100,000 newborns.}
 \item{dx}{Lifetable function dx, which means number of people dying between age x and x+1 (for lifetable(), x+class) from 100,000 newborns. Differentials of lx.}
 \item{ux}{The force of mortality.}
 \item{pix}{age-(class-)specific proportions of unhealthy people.}
 \item{Nx}{Population of xth age-class, which is needed to calculate confidence intervals.}
 \item{conf.level}{The level of confidence intervals.  Default is 0.95.}
 \item{Tx}{Lifetable function Tx, which means sum of person-years lived above age x.}
}
\value{
 \item{ages}{Lifetable's exact age x, which is the beginning of each interval.}
 \item{n}{Duration of each interval.  If ns is not given, the value of the class is repeatedly used.}
 \item{mx}{Lifetable function mx, meaning the age (class)-specific death rates.}
 \item{qx}{Lifetable function qx, which means the probability of dying between age x and x+1 (for lifetable(), x+class).}
 \item{ax}{Lifetable function ax, which means the average number of person-years lived in the interval by those dying in the interval.  In lifetable(), it's automatically specified by mode.}
 \item{lx}{Lifetable function lx, which means number of people left alive at age x from 100,000 newborns.}
 \item{dx}{Lifetable function dx, which means number of people dying between age x and x+1 (for lifetable(), x+class) from 100,000 newborns. Differentials of lx.}
 \item{Lx}{Lifetable function Lx, which means person-years lived between age x and x+class.}
 \item{Tx}{Lifetable function Tx, which means person-years lived above age x.}
 \item{ex}{Lifetable function ex, which means expectation of life at age x.}
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\references{
 Preston SH, Heuveline P, Guillot M (2001) Demography: Measuring and Modeling 
 Population Processes.  Blackwell Publishing, Oxford.
 
 Newell C (1988) Methods and Models in Demography.  The Guilford Press, New York.
 
 Sullivan DF (1971) A single index of mortality and morbidity. 
 HSMHA Health Reports, 86: 347-354.
}
\seealso{\link[fmsb]{Jlife}}
\examples{
 lifetable(c(0.0087, 0.00015, 0.00019, 0.00098, 0.0013, 0.0011, 0.0014, 0.0019, 
             0.0029, 0.0048, 0.0071, 0.011, 0.019, 0.028, 0.041, 0.072, 0.11, 
             0.19), class=5, mode=11)
 lifetable2(c(0.008314, 0.000408, 0.000181, 0.000187, 0.000282, 0.000307, 0.000364, 
              0.000566, 0.000884, 0.001445, 0.002485, 0.004210, 0.007219, 0.012054, 
              0.018259, 0.029920, 0.049689, 0.085545, 0.177987), 
              ax = c(0.1, 0.4, rep(0.5, 16), NA), n = c(1, 4, rep(5, 16), NA) )
 lifetable3(lx=c(1.0, 0.8499, 0.8070, 0.7876, 0.7762, 0.7691, 0.7502, 0.7362,
                 0.7130, 0.6826, 0.6525, 0.6223, 0.5898, 0.5535, 0.5106, 0.4585,
                 0.3965, 0.3210, 0.2380, 0.1516, 0.0768, 0.0276, 0.0059, 0.0006, 0),
                 n=c(rep(1, 5), rep(5, 20)), ax=c(0.3, rep(0.5, 24))) # Newell, Table 13.1
 clifetable(Jlife$qx2000F)
}
