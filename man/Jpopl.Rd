\name{Jpopl}
\docType{data}
\alias{Jpopl}
\title{Japanese population data (unshrunken version)}
\description{
 The data gives the sex and age specific population for the all census 
 results in Japan.
}
\usage{ Jpopl }
\format{ A data frame with 111 observations on 67 variables.

 \tabular{rlll}{
  \code{[, 1]} \tab \code{Age} \tab Factor w/111 levels \tab Ages (years old, combined for 110+)\cr
  \code{[, 2]} \tab \code{M1888} \tab numeric \tab Age specific population of males in 1888\cr
  \code{[, 3]} \tab \code{F1888} \tab numeric \tab Age specific population of females in 1888\cr
  \code{[, 4]} \tab \code{M1893} \tab numeric \tab Age specific population of males in 1893\cr
  \code{[, 5]} \tab \code{F1893} \tab numeric \tab Age specific population of females in 1893\cr
  \code{[, 6]} \tab \code{M1898} \tab numeric \tab Age specific population of males in 1898\cr
  \code{[, 7]} \tab \code{F1898} \tab numeric \tab Age specific population of females in 1898\cr
  \code{[, 8]} \tab \code{M1903} \tab numeric \tab Age specific population of males in 1903\cr
  \code{[, 9]} \tab \code{F1903} \tab numeric \tab Age specific population of females in 1903\cr
  \code{[,10]} \tab \code{M1908} \tab numeric \tab Age specific population of males in 1908\cr
  \code{[,11]} \tab \code{F1908} \tab numeric \tab Age specific population of females in 1908\cr
  \code{[,12]} \tab \code{M1913} \tab numeric \tab Age specific population of males in 1913\cr
  \code{[,13]} \tab \code{F1913} \tab numeric \tab Age specific population of females in 1913\cr
  \code{[,14]} \tab \code{M1918} \tab numeric \tab Age specific population of males in 1918\cr
  \code{[,15]} \tab \code{F1918} \tab numeric \tab Age specific population of females in 1918\cr
  \code{[,16]} \tab \code{M1920} \tab numeric \tab Age specific population of males in 1920\cr
  \code{[,17]} \tab \code{F1920} \tab numeric \tab Age specific population of females in 1920\cr
  \code{[,18]} \tab \code{M1925} \tab numeric \tab Age specific population of males in 1925\cr
  \code{[,19]} \tab \code{F1925} \tab numeric \tab Age specific population of females in 1925\cr
  \code{[,20]} \tab \code{M1930} \tab numeric \tab Age specific population of males in 1930\cr
  \code{[,21]} \tab \code{F1930} \tab numeric \tab Age specific population of females in 1930\cr
  \code{[,22]} \tab \code{M1935} \tab numeric \tab Age specific population of males in 1935\cr
  \code{[,23]} \tab \code{F1935} \tab numeric \tab Age specific population of females in 1935\cr
  \code{[,24]} \tab \code{M1940} \tab numeric \tab Age specific population of males in 1940\cr
  \code{[,25]} \tab \code{F1940} \tab numeric \tab Age specific population of females in 1940\cr
  \code{[,26]} \tab \code{M1947} \tab numeric \tab Age specific population of males in 1947\cr
  \code{[,27]} \tab \code{F1947} \tab numeric \tab Age specific population of females in 1947\cr
  \code{[,28]} \tab \code{M1950} \tab numeric \tab Age specific population of males in 1950\cr
  \code{[,29]} \tab \code{F1950} \tab numeric \tab Age specific population of females in 1950\cr
  \code{[,30]} \tab \code{M1955} \tab numeric \tab Age specific population of males in 1955\cr
  \code{[,31]} \tab \code{F1955} \tab numeric \tab Age specific population of females in 1955\cr
  \code{[,32]} \tab \code{M1960} \tab numeric \tab Age specific population of males in 1960\cr
  \code{[,33]} \tab \code{F1960} \tab numeric \tab Age specific population of females in 1960\cr
  \code{[,34]} \tab \code{M1965} \tab numeric \tab Age specific population of males in 1965\cr
  \code{[,35]} \tab \code{F1965} \tab numeric \tab Age specific population of females in 1965\cr
  \code{[,36]} \tab \code{M1970} \tab numeric \tab Age specific population of males in 1970\cr
  \code{[,37]} \tab \code{F1970} \tab numeric \tab Age specific population of females in 1970\cr
  \code{[,38]} \tab \code{M1975} \tab numeric \tab Age specific population of males in 1975\cr
  \code{[,39]} \tab \code{F1975} \tab numeric \tab Age specific population of females in 1975\cr
  \code{[,40]} \tab \code{M1980} \tab numeric \tab Age specific population of males in 1980\cr
  \code{[,41]} \tab \code{F1980} \tab numeric \tab Age specific population of females in 1980\cr
  \code{[,42]} \tab \code{M1985} \tab numeric \tab Age specific population of males in 1985\cr
  \code{[,43]} \tab \code{F1985} \tab numeric \tab Age specific population of females in 1985\cr
  \code{[,44]} \tab \code{M1990} \tab numeric \tab Age specific population of males in 1990\cr
  \code{[,45]} \tab \code{F1990} \tab numeric \tab Age specific population of females in 1990\cr
  \code{[,46]} \tab \code{M1995} \tab numeric \tab Age specific population of males in 1995\cr
  \code{[,47]} \tab \code{F1995} \tab numeric \tab Age specific population of females in 1995\cr
  \code{[,48]} \tab \code{M2000} \tab numeric \tab Age specific population of males in 2000\cr
  \code{[,49]} \tab \code{F2000} \tab numeric \tab Age specific population of females in 2000\cr
  \code{[,50]} \tab \code{M2000J} \tab numeric \tab Age specific population of Japanese males in 2000\cr
  \code{[,51]} \tab \code{F2000J} \tab numeric \tab Age specific population of Japanese females in 2000\cr
  \code{[,52]} \tab \code{M2005} \tab numeric \tab Age specific population of males in 2005\cr
  \code{[,53]} \tab \code{F2005} \tab numeric \tab Age specific population of females in 2005\cr
  \code{[,54]} \tab \code{M2005J} \tab numeric \tab Age specific population of Japanese males in 2005\cr
  \code{[,55]} \tab \code{F2005J} \tab numeric \tab Age specific population of Japanese females in 2005\cr
  \code{[,56]} \tab \code{M2010} \tab numeric \tab Age specific population of males in 2010\cr
  \code{[,57]} \tab \code{F2010} \tab numeric \tab Age specific population of females in 2010\cr
  \code{[,58]} \tab \code{M2010J} \tab numeric \tab Age specific population of Japanese males in 2010\cr
  \code{[,59]} \tab \code{F2010J} \tab numeric \tab Age specific population of Japanese females in 2010\cr
  \code{[,60]} \tab \code{M2015} \tab numeric \tab Age specific population of males in 2015\cr
  \code{[,61]} \tab \code{F2015} \tab numeric \tab Age specific population of females in 2015\cr
  \code{[,62]} \tab \code{M2015J} \tab numeric \tab Age specific population of Japanese males in 2015\cr
  \code{[,63]} \tab \code{F2015J} \tab numeric \tab Age specific population of Japanese females in 2015\cr
  \code{[,64]} \tab \code{M2020} \tab numeric \tab Age specific population of males in 2020\cr
  \code{[,65]} \tab \code{F2020} \tab numeric \tab Age specific population of females in 2020\cr
  \code{[,66]} \tab \code{M2020J} \tab numeric \tab Age specific population of Japanese males in 2020\cr
  \code{[,67]} \tab \code{F2020J} \tab numeric \tab Age specific population of Japanese females in 2020\cr
 }
}
\details{
   Japanese population data by sex and age given as national official census record.
   
   \itemize{
    \item \code{Age}: Ages, combined for 110+.
    \item \code{M1888}-\code{M2020}: Age specific number of males' population in Japan for 1888-2020.
    \item \code{F1888}-\code{F2020}: Age specific number of females' population in Japan for 1888-2020.
    \item \code{M2000J}-\code{M2020J}: Age specific number of Japanese males' population in Japan for 2000-2020 by every 5 years.
    \item \code{F2000J}-\code{F2020J}: Age specific number of Japanese females' population in Japan for 2000-2020 by every 5 years.
   }
}
\source{
 \url{https://www.stat.go.jp/english/data/kokusei/index.html}
 \url{https://warp.da.ndl.go.jp/info:ndljp/pid/1334623/www.stat.go.jp/english/data/chouki/02.htm}
 \url{https://www.e-stat.go.jp/stat-search/files/data?fileid=000007809775&rcount=3}
 \url{https://www.e-stat.go.jp/stat-search/file-download?statInfId=000032142404&fileKind=0}
}
\references{
  Statistics Bureau, Ministry of Internal Affairs and Communications: Population Census, 1888-2020.
}
\keyword{datasets}
