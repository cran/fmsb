\name{Jfert}
\docType{data}
\alias{Jfert}
\title{Age-specific fertility rates and age-specific marital fertility rates in Japan.}
\description{
  Age-specific fertility and marital fertility rates for aged 15-54 Japanese wowmen in Japan, from 1950 to 2020, every five years.
}
\usage{ Jfert }
\format{ A data frame with 40 observations on 31 variables.

 \tabular{rlll}{
  \code{[, 1]} \tab \code{Age} \tab integer \tab Ages of women\cr
  \code{[, 2]} \tab \code{ASFR1950} \tab numeric \tab Age-specific fertility rates of Japanese women in 1950.\cr
  \code{[, 3]} \tab \code{ASFR1955} \tab numeric \tab Age-specific fertility rates of Japanese women in 1955.\cr
  \code{[, 4]} \tab \code{ASFR1960} \tab numeric \tab Age-specific fertility rates of Japanese women in 1960.\cr
  \code{[, 5]} \tab \code{ASFR1965} \tab numeric \tab Age-specific fertility rates of Japanese women in 1965.\cr
  \code{[, 6]} \tab \code{ASFR1970} \tab numeric \tab Age-specific fertility rates of Japanese women in 1970.\cr
  \code{[, 7]} \tab \code{ASFR1975} \tab numeric \tab Age-specific fertility rates of Japanese women in 1975.\cr
  \code{[, 8]} \tab \code{ASFR1980} \tab numeric \tab Age-specific fertility rates of Japanese women in 1980.\cr
  \code{[, 9]} \tab \code{ASFR1985} \tab numeric \tab Age-specific fertility rates of Japanese women in 1985.\cr
  \code{[,10]} \tab \code{ASFR1990} \tab numeric \tab Age-specific fertility rates of Japanese women in 1990.\cr
  \code{[,11]} \tab \code{ASFR1995} \tab numeric \tab Age-specific fertility rates of Japanese women in 1995.\cr
  \code{[,12]} \tab \code{ASFR2000} \tab numeric \tab Age-specific fertility rates of Japanese women in 2000.\cr
  \code{[,13]} \tab \code{ASFR2005} \tab numeric \tab Age-specific fertility rates of Japanese women in 2005.\cr
  \code{[,14]} \tab \code{ASFR2010} \tab numeric \tab Age-specific fertility rates of Japanese women in 2010.\cr
  \code{[,15]} \tab \code{ASFR2015} \tab numeric \tab Age-specific fertility rates of Japanese women in 2015.\cr
  \code{[,16]} \tab \code{ASFR2020} \tab numeric \tab Age-specific fertility rates of Japanese women in 2020.\cr
  \code{[,17]} \tab \code{ASMFR1950} \tab numeric \tab Age-specific marital fertility rates of Japanese married women in 1950.\cr
  \code{[,18]} \tab \code{ASMFR1955} \tab numeric \tab Age-specific marital fertility rates of Japanese married women in 1955.\cr
  \code{[,19]} \tab \code{ASMFR1960} \tab numeric \tab Age-specific marital fertility rates of Japanese married women in 1960.\cr
  \code{[,20]} \tab \code{ASMFR1965} \tab numeric \tab Age-specific marital fertility rates of Japanese married women in 1965.\cr
  \code{[,21]} \tab \code{ASMFR1970} \tab numeric \tab Age-specific marital fertility rates of Japanese married women in 1970.\cr
  \code{[,22]} \tab \code{ASMFR1975} \tab numeric \tab Age-specific marital fertility rates of Japanese married women in 1975.\cr
  \code{[,23]} \tab \code{ASMFR1980} \tab numeric \tab Age-specific marital fertility rates of Japanese married women in 1980.\cr
  \code{[,24]} \tab \code{ASMFR1985} \tab numeric \tab Age-specific marital fertility rates of Japanese married women in 1985.\cr
  \code{[,25]} \tab \code{ASMFR1990} \tab numeric \tab Age-specific marital fertility rates of Japanese married women in 1990.\cr
  \code{[,26]} \tab \code{ASMFR1995} \tab numeric \tab Age-specific marital fertility rates of Japanese married women in 1995.\cr
  \code{[,27]} \tab \code{ASMFR2000} \tab numeric \tab Age-specific marital fertility rates of Japanese married women in 2000.\cr
  \code{[,28]} \tab \code{ASMFR2005} \tab numeric \tab Age-specific marital fertility rates of Japanese married women in 2005.\cr
  \code{[,29]} \tab \code{ASMFR2010} \tab numeric \tab Age-specific marital fertility rates of Japanese married women in 2010.\cr
  \code{[,30]} \tab \code{ASMFR2015} \tab numeric \tab Age-specific marital fertility rates of Japanese married women in 2015.\cr
  \code{[,31]} \tab \code{ASMFR2020} \tab numeric \tab Age-specific marital fertility rates of Japanese married women in 2020.\cr
 }
}
\details{
  The calculations were the numbers of live births divided by the numbers of women for ASFR (15-54),
  and the numbers of legitimate live births divided by the numbers of married women for ASMFR (15-54).
  Data sources are all official publication as vital statistics and national population census, so that the data are given with 5 years intervals.
   
   \itemize{
    \item \code{Age}: Ages of women, from 15 to 54.
    \item \code{ASFR1950}-\code{ASFR2020}: Age-specific fertility rates for all women aged 15-54 for 1950-2020, every 5 years.
    \item \code{ASMFR1950}-\code{ASMFR2020}: Age-specific marital fertility rates for married women aged 15-54 for 1950-2020, every 5 years.
   }
}
\source{
 \url{https://warp.da.ndl.go.jp/info:ndljp/pid/1334623/www.stat.go.jp/english/data/chouki/02.htm}
 \url{https://warp.da.ndl.go.jp/info:ndljp/pid/1334623/www.stat.go.jp/data/chouki/zuhyou/02-29-b.xls}
 \url{https://www.ipss.go.jp/syoushika/tohkei/Popular/P_Detail2022.asp?fname=T04-09.htm}
 \url{https://www.e-stat.go.jp/stat-search/file-download?statInfId=000032118572&fileKind=1}
 \url{https://www.e-stat.go.jp/stat-search/file-download?statInfId=000032142474&fileKind=0}
}
\references{
  Ministry of Health, Labor and Welfare of Japan: Vital Statistics. / Ministry of Internal Affairs and Communications, Statistics Bureau: Population Census.
}
\keyword{datasets}
