\name{Prefe0}
\docType{data}
\alias{Prefe0}
\title{Changes of life expectancy at birth for each prefecture in Japan since 1965}
\description{
 The data gives the estimates of life expectancy at birth (e0) for
 each prefecture in Japan since 1965.
}
\usage{Prefe0}
\format{ A data frame with 47 observations on 26 variables.

 \tabular{rlll}{
  \code{[, 1]} \tab \code{PNAME} \tab factor w/47 levels \tab The name (in roma-ji) for prefectures\cr
  \code{[, 2]} \tab \code{JCODE} \tab numeric \tab Prefecture number defined by Geographical Information 
   Authority of Japan\cr
  \code{[, 3]} \tab \code{E0M.1965} \tab numeric \tab Life expectancy at birth of each prefecture for males in 1965\cr
  \code{[, 4]} \tab \code{E0M.1970} \tab numeric \tab Life expectancy at birth of each prefecture for males in 1970\cr
  \code{[, 5]} \tab \code{E0M.1975} \tab numeric \tab Life expectancy at birth of each prefecture for males in 1975\cr
  \code{[, 6]} \tab \code{E0M.1980} \tab numeric \tab Life expectancy at birth of each prefecture for males in 1980\cr
  \code{[, 7]} \tab \code{E0M.1985} \tab numeric \tab Life expectancy at birth of each prefecture for males in 1985\cr
  \code{[, 8]} \tab \code{E0M.1990} \tab numeric \tab Life expectancy at birth of each prefecture for males in 1990\cr
  \code{[, 9]} \tab \code{E0M.1995} \tab numeric \tab Life expectancy at birth of each prefecture for males in 1995\cr
  \code{[,10]} \tab \code{E0M.2000} \tab numeric \tab Life expectancy at birth of each prefecture for males in 2000\cr
  \code{[,11]} \tab \code{E0M.2005} \tab numeric \tab Life expectancy at birth of each prefecture for males in 2005\cr
  \code{[,12]} \tab \code{E0M.2010} \tab numeric \tab Life expectancy at birth of each prefecture for males in 2010\cr
  \code{[,13]} \tab \code{E0M.2015} \tab numeric \tab Life expectancy at birth of each prefecture for males in 2015\cr
  \code{[,14]} \tab \code{E0M.2020} \tab numeric \tab Life expectancy at birth of each prefecture for males in 2020\cr
  \code{[,15]} \tab \code{E0F.1965} \tab numeric \tab Life expectancy at birth of each prefecture for females in 1965\cr
  \code{[,16]} \tab \code{E0F.1970} \tab numeric \tab Life expectancy at birth of each prefecture for females in 1970\cr
  \code{[,17]} \tab \code{E0F.1975} \tab numeric \tab Life expectancy at birth of each prefecture for females in 1975\cr
  \code{[,18]} \tab \code{E0F.1980} \tab numeric \tab Life expectancy at birth of each prefecture for females in 1980\cr
  \code{[,19]} \tab \code{E0F.1985} \tab numeric \tab Life expectancy at birth of each prefecture for females in 1985\cr
  \code{[,20]} \tab \code{E0F.1990} \tab numeric \tab Life expectancy at birth of each prefecture for females in 1990\cr
  \code{[,21]} \tab \code{E0F.1995} \tab numeric \tab Life expectancy at birth of each prefecture for females in 1995\cr
  \code{[,22]} \tab \code{E0F.2000} \tab numeric \tab Life expectancy at birth of each prefecture for females in 2000\cr
  \code{[,23]} \tab \code{E0F.2005} \tab numeric \tab Life expectancy at birth of each prefecture for females in 2005\cr
  \code{[,24]} \tab \code{E0F.2010} \tab numeric \tab Life expectancy at birth of each prefecture for females in 2010\cr
  \code{[,25]} \tab \code{E0F.2015} \tab numeric \tab Life expectancy at birth of each prefecture for females in 2015\cr
  \code{[,26]} \tab \code{E0F.2020} \tab numeric \tab Life expectancy at birth of each prefecture for females in 2020\cr
 }
}
\details{
   Life expectancy at birth for each prefecture in Japan since 1965.
   
   \itemize{
    \item \code{PNAME}: The name (in roma-ji) for prefectures.
    \item \code{JCODE}: Prefecture number defined by Geographical Information Authority of Japan. 
     From 1 to 47.
    \item \code{E0[M|F].*}: Life expectancy at birth (e0) of each prefecture for males ([M]) or for females ([F]) in year (*).
   }
}
\source{
 \url{https://www.mhlw.go.jp/toukei/saikin/hw/life/tdfk20/dl/tdfk20-08.xls}
}
\references{
  Ministry of Health, Labor and Welfare of Japan: Vital Statistics with Life Expectancy 2020.
  \url{https://minato.sip21c.org/demography/how-to-make-pref-charts.html} (in Japanese), 
  WHO \url{https://www.who.int/data/gho/data/themes/mortality-and-global-health-estimates}
}
\examples{
require(fmsb)
x <- Prefe0
males <- t(x[, 3:14])
colnames(males) <- x$PNAME
females <- t(x[, 15:26])
colnames(females) <- x$PNAME
COL <- ifelse(x$PNAME=="Nagano", "blue", ifelse(x$PNAME=="Okinawa", "red", "lightgrey"))
LWD <- ifelse(x$PNAME=="Nagano", 2, ifelse(x$PNAME=="Okinawa", 2, 1))
LTY <- ifelse(x$PNAME=="Nagano", 1, ifelse(x$PNAME=="Okinawa", 1, 3))
years <- 1965+0:11*5
layout(t(1:2))
matplot(years, males, type="l", col=COL, lwd=LWD, lty=LTY, 
 main="Changes of e0 for males in each prefecture of Japan
 (Blue: Nagano, Red: Okinawa, Grey: Other)")
matplot(years, females, type="l", col=COL, lwd=LWD, lty=LTY, 
 main="Changes of e0 for females in each prefecture of Japan
 (Blue: Nagano, Red: Okinawa, Grey: Other)")
}
\keyword{datasets}
