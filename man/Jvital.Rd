\name{Jvital}
\docType{data}
\alias{Jvital}
\title{Longitudinal data of several vital statistics in Japan}
\description{
 The data gives longitudinal data of several vital statistics in Japan.
 Included indices are crude birth rates, crude death rates, infant mortality
 rates, and so on.
}
\usage{ Jvital }
\format{ A data frame with 108 observations on 16 variables.

 \tabular{rlll}{
  \code{[, 1]} \tab \code{YEAR} \tab numeric \tab Year\cr
  \code{[, 2]} \tab \code{CBR} \tab numeric \tab Crude birth rates of Japan\cr
  \code{[, 3]} \tab \code{CDR} \tab numeric \tab Crude death rates of Japan\cr
  \code{[, 4]} \tab \code{IMR} \tab numeric \tab Infant mortality rates of Japan\cr
  \code{[, 5]} \tab \code{NMR} \tab numeric \tab Neonatal mortality rates of Japan\cr
  \code{[, 6]} \tab \code{NIR} \tab numeric \tab Natural increase rates of Japan\cr
  \code{[, 7]} \tab \code{SBRPB} \tab numeric \tab Stillbirth rates of Japan\cr
  \code{[, 8]} \tab \code{SARPB} \tab numeric \tab Spontaneous abortion rates of Japan\cr
  \code{[, 9]} \tab \code{ACRPB} \tab numeric \tab Artificial contraception rates of Japan\cr
  \code{[,10]} \tab \code{PNMPB} \tab numeric \tab Perinatal mortalities per birth of Japan\cr
  \code{[,11]} \tab \code{MR} \tab numeric \tab Marriage rates of Japan\cr
  \code{[,12]} \tab \code{DR} \tab numeric \tab Divorce rates of Japan\cr
  \code{[,13]} \tab \code{TFR} \tab numeric \tab Total fertility rates of Japan\cr
  \code{[,14]} \tab \code{ASMRM} \tab numeric \tab Age-standardized mortality rates of males in Japan\cr
  \code{[,15]} \tab \code{ASMRF} \tab numeric \tab Age-standardized mortality rates of females in Japan\cr
  \code{[,16]} \tab \code{PNMPLB} \tab numeric \tab Perinatal mortalities per live births of Japan\cr
 }
}
\details{
   Longitudinal vital statistics in Japan provided as national
   official vital statitistics every year from 1899 to 2009,
   except for 1944-1946.
   
   \itemize{
    \item \code{YEAR}: Calender year.
    \item \code{CBR}: Crude birth rate.  Number of all live birth / mid-year population 1000.
    \item \code{CDR}: Crude death rate.  Number of death / mid-year population 1000.
    \item \code{IMR}: Infant mortality rate. Number of death at age 0 / 1000 live births.
    \item \code{NMR}: Neonatal mortality rate. Number of death within 4 weeks after birth / 1000 live births.
    \item \code{NIR}: Natural increase rate.  \code{CBR}-\code{CDR}.
    \item \code{SBRPB}: Stillbirth rate per birth.  Number of stillbirths / 1000 births.
    \item \code{SARPB}: Spontaneous abortion rate per birth.  Number of spontaneous abortions / 1000 births.
    \item \code{ACRPB}: Artificial contraception (= induced abortion) rate per birth.  Number of induced abortions / 1000 births.
    \item \code{PNMPB}: Perinatal mortality per birth.  [(Number of stillbirths after gestational age 22 weeks) + (Number of early neonatal deaths within a week after birth)] per 1000 births.  The denominator is the sum of the number of stillbirths after gestational age 22 weeks and the number of live births.  This definition was established in 1995, but \code{PNMPB} also includes some values before 1995.
    \item \code{MR}: Marital rate.  The number of marriages / mid-year population 1000.
    \item \code{DR}: Divorce rate.  The number of divorces / mid-year population 1000.
    \item \code{TFR}: Total fertility rate.  The sum of age-specific fertility rates, which is the number of births divided by the number of women's population for each age.
    \item \code{ASMRM}: Age-standardized mortality rate of males, per mid-year population 1000, where the standard population is the model population in 1985.
    \item \code{ASMRF}: Age-standardized mortality rate of females, per mid-year population 1000, where the standard population is the model population in 1985.
    \item \code{PNMPLB}: Perinatal mortality per live birth.  [(Number of stillbirths after gestational age 28 weeks) + (Number of early neonatal deaths within a week after birth)] per 1000 live births.  This definition stood until 1994, but \code{PNMPLB} also includes values after 1995, for comparison.
   }
}
\source{
 \url{http://www.stat.go.jp/english/data/chouki/02.htm}
 \url{http://www.e-stat.go.jp/SG1/estat/List.do?bid=000001025191&amp;cycode=0}
}
\references{
  Ministry of Health, Labor and Welfare of Japan: Vital Statistics.
}
\keyword{datasets}
