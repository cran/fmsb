\name{lifetable}
\alias{lifetable}
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
\title{Lifetable functions with mutual conversions}
\description{
 Lifetable related functions.
}
\usage{
 lifetable(mx, ns=NULL, class=5, mode=1)
 clifetable(qx)
 lxtodx(lx)
 dxtolx(dx)
 qxtodx(qx)
 dxtoqx(dx)
 qxtomx(qx)
 mxtoqx(mx)
 qxtolx(qx)
 lxtoqx(lx)
 uxtoqx(ux)
}
\arguments{
 \item{mx}{Lifetable function mx, meaning the age (class)-specific death rates.}
 \item{ns}{If given as a vector with the same length as mx, the duration for each age (class). Default is NULL: same duration with class is automatically used.}
 \item{class}{Age-class of lifetable() function.  Default is 5.}
 \item{mode}{How to set ax and correction method in conversion from mx to qx.  1 and 11: all ax is 0.5 except the open-ended class [where ax is reciprocal of mx], 2, 4, 12, 14: ax is 0.1 for age 0, 0.4 for age 1-4, 0.5 for the other ages except the open-ended class [where ax is reciprocal of mx], 3, 5, 13, 15: ax is 0.3 for age 0, 0.4 for age 1-4, 0.5 for the other ages except the open-ended class [where ax is reciprocal of mx], 6 and 16: Males value given in Preston SH (2001), pp.48 Table 3.3, 7 and 17: Females value given in Preston SH (2001), pp.48 Table 3.3.  If less than 10, simply calculating qx as n*mx/(1+n*(1-ax)*mx) except for the open-ended class where qx=1, otherwise calculating qx by Greville's method.  Default is 1.}
 \item{qx}{Lifetable function qx, which means the probability of dying between age x and x+1 (for lifetable(), x+class).}
 \item{lx}{Lifetable function lx, which means number of people left alive at age x from 100,000 newborns.}
 \item{dx}{Lifetable function dx, which means number of people dying between age x and x+1 (for lifetable(), x+class) from 100,000 newborns. Differentials of lx.}
 \item{ux}{The force of mortality.}
}
\value{
 \item{ages}{Lifetable's exact age x, which is the beginning of each interval.}
 \item{n}{Duration of each interval.  If ns is not given, the value of the class is repeatedly used.}
 \item{mx}{Lifetable function mx, meaning the age (class)-specific death rates.}
 \item{qx}{Lifetable function qx, which means the probability of dying between age x and x+1 (for lifetable(), x+class).}
 \item{ax}{Lifetable function ax, which means the average number of person-years lived in the interval by those dying in the interval.  Automatically specified by mode.}
 \item{lx}{Lifetable function lx, which means number of people left alive at age x from 100,000 newborns.}
 \item{dx}{Lifetable function dx, which means number of people dying between age x and x+1 (for lifetable(), x+class) from 100,000 newborns. Differentials of lx.}
 \item{Lx}{Lifetable function Lx, which means person-years lived between age x and x+class.}
 \item{Tx}{Lifetable function Tx, which means person-years lived above age x.}
 \item{ex}{Lifetable function ex, which means expectation of life at age x.}
}
\author{Minato Nakazawa \email{minato-nakazawa@umin.net} \url{http://minato.sip21c.org/}}
\references{
 Preston SH, Heuveline P, Guillot M (2001) Demography: Measuring and Modeling Population Processes.  Blackwell Publishing, Oxford.
 
 Newell C (1988) Methods and Models in Demography.  The Guilford Press, New York.
}
\seealso{\link[fmsb]{Jlife}}
\examples{
 lifetable(c(0.0087, 0.00015, 0.00019, 0.00098, 0.0013, 0.0011, 0.0014, 0.0019, 0.0029, 
             0.0048, 0.0071, 0.011, 0.019, 0.028, 0.041, 0.072, 0.11, 0.19), class=5, mode=11)
 clifetable(Jlife$qx2000F)
}
