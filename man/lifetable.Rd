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
 lifetable(mx, class=5, mode=1)
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
 \item{class}{Age-class of lifetable() function.  Default is 5.}
 \item{mode}{Correction method in conversion from mx to qx.  Default is 1, which means simple calculation, otherwise, Greville' correction is applied.}
 \item{qx}{Lifetable function qx, which means the probability of dying between age x and x+1 (for lifetable(), x+class).}
 \item{lx}{Lifetable function lx, which means number of people left alive at age x from 100,000 newborns.}
 \item{dx}{Lifetable function dx, which means number of people dying between age x and x+1 (for lifetable(), x+class) from 100,000 newborns. Differentials of lx.}
 \item{ux}{The force of mortality.}
}
\value{
 \item{mx}{Lifetable function mx, meaning the age (class)-specific death rates.}
 \item{qx}{Lifetable function qx, which means the probability of dying between age x and x+1 (for lifetable(), x+class).}
 \item{lx}{Lifetable function lx, which means number of people left alive at age x from 100,000 newborns.}
 \item{dx}{Lifetable function dx, which means number of people dying between age x and x+1 (for lifetable(), x+class) from 100,000 newborns. Differentials of lx.}
 \item{Lx}{Lifetable function Lx, which means person-years lived between age x and x+class.}
 \item{Tx}{Lifetable function Tx, which means person-years lived above age x.}
 \item{ex}{Lifetable function ex, which means expectation of life at age x.}
}
\author{Minato Nakazawa \email{minato-nakazawa@umin.net} \url{http://phi.med.gunma-u.ac.jp/}}
\seealso{\link[fmsb]{Jlife}}
\examples{
 lifetable(c(0.0087,0.00015,0.00019,0.00098,0.0013,0.0011,0.0014,0.0019,0.0029,0.0048,0.0071,0.011,0.019,0.028,0.041,0.072,0.11,0.19),class=5,mode=2)
 clifetable(Jlife$qx2000F)
}
