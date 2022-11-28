\name{Jlife}
\docType{data}
\alias{Jlife}
\title{Completed lifetables in Japan.}
\description{
 The qx column of the completed lifetables in Japan, from "1891-1898" to "2020", mostly every 5 years.
}
\usage{ Jlife }
\format{ A data frame with 117 observations (NAs are filled for the ages with no survivors) on 45 variables.

 \tabular{rlll}{
  \code{[, 1]} \tab \code{Age} \tab integer \tab Ages of women\cr
  \code{[, 2]} \tab \code{qx1895M} \tab numeric \tab qx of completed lifetable functions of Japanese men in 1891-1898.\cr
  \code{[, 3]} \tab \code{qx1895F} \tab numeric \tab qx of completed lifetable functions of Japanese women in 1891-1898.\cr
  \code{[, 4]} \tab \code{qx1901M} \tab numeric \tab qx of completed lifetable functions of Japanese men in 1899-1903.\cr
  \code{[, 5]} \tab \code{qx1901F} \tab numeric \tab qx of completed lifetable functions of Japanese women in 1899-1903.\cr
  \code{[, 6]} \tab \code{qx1911M} \tab numeric \tab qx of completed lifetable functions of Japanese men in 1909-1913.\cr
  \code{[, 7]} \tab \code{qx1911F} \tab numeric \tab qx of completed lifetable functions of Japanese women in 1909-1913.\cr
  \code{[, 8]} \tab \code{qx1923M} \tab numeric \tab qx of completed lifetable functions of Japanese men in 1921-1925.\cr
  \code{[, 9]} \tab \code{qx1923F} \tab numeric \tab qx of completed lifetable functions of Japanese women in 1921-1925.\cr
  \code{[,10]} \tab \code{qx1928M} \tab numeric \tab qx of completed lifetable functions of Japanese men in 1926-1930.\cr
  \code{[,11]} \tab \code{qx1928F} \tab numeric \tab qx of completed lifetable functions of Japanese women in 1926-1930.\cr
  \code{[,12]} \tab \code{qx1935M} \tab numeric \tab qx of completed lifetable functions of Japanese men in 1935-1936.\cr
  \code{[,13]} \tab \code{qx1935F} \tab numeric \tab qx of completed lifetable functions of Japanese women in 1935-1936.\cr
  \code{[,14]} \tab \code{qx1947M} \tab numeric \tab qx of completed lifetable functions of Japanese men in 1947.\cr
  \code{[,15]} \tab \code{qx1947F} \tab numeric \tab qx of completed lifetable functions of Japanese women in 1947.\cr
  \code{[,16]} \tab \code{qx1951M} \tab numeric \tab qx of completed lifetable functions of Japanese men in 1950-1952.\cr
  \code{[,17]} \tab \code{qx1951F} \tab numeric \tab qx of completed lifetable functions of Japanese women in 1950-1952.\cr
  \code{[,18]} \tab \code{qx1955M} \tab numeric \tab qx of completed lifetable functions of Japanese men in 1955.\cr
  \code{[,19]} \tab \code{qx1955F} \tab numeric \tab qx of completed lifetable functions of Japanese women in 1955.\cr
  \code{[,20]} \tab \code{qx1960M} \tab numeric \tab qx of completed lifetable functions of Japanese men in 1960.\cr
  \code{[,21]} \tab \code{qx1960F} \tab numeric \tab qx of completed lifetable functions of Japanese women in 1960.\cr
  \code{[,22]} \tab \code{qx1965M} \tab numeric \tab qx of completed lifetable functions of Japanese men in 1965.\cr
  \code{[,23]} \tab \code{qx1965F} \tab numeric \tab qx of completed lifetable functions of Japanese women in 1965.\cr
  \code{[,24]} \tab \code{qx1970M} \tab numeric \tab qx of completed lifetable functions of Japanese men in 1970.\cr
  \code{[,25]} \tab \code{qx1970F} \tab numeric \tab qx of completed lifetable functions of Japanese women in 1970.\cr
  \code{[,26]} \tab \code{qx1975M} \tab numeric \tab qx of completed lifetable functions of Japanese men in 1975.\cr
  \code{[,27]} \tab \code{qx1975F} \tab numeric \tab qx of completed lifetable functions of Japanese women in 1975.\cr
  \code{[,28]} \tab \code{qx1980M} \tab numeric \tab qx of completed lifetable functions of Japanese men in 1980.\cr
  \code{[,29]} \tab \code{qx1980F} \tab numeric \tab qx of completed lifetable functions of Japanese women in 1980.\cr
  \code{[,30]} \tab \code{qx1985M} \tab numeric \tab qx of completed lifetable functions of Japanese men in 1985.\cr
  \code{[,31]} \tab \code{qx1985F} \tab numeric \tab qx of completed lifetable functions of Japanese women in 1985.\cr
  \code{[,32]} \tab \code{qx1990M} \tab numeric \tab qx of completed lifetable functions of Japanese men in 1990.\cr
  \code{[,33]} \tab \code{qx1990F} \tab numeric \tab qx of completed lifetable functions of Japanese women in 1990.\cr
  \code{[,34]} \tab \code{qx1995M} \tab numeric \tab qx of completed lifetable functions of Japanese men in 1995.\cr
  \code{[,35]} \tab \code{qx1995F} \tab numeric \tab qx of completed lifetable functions of Japanese women in 1995.\cr
  \code{[,36]} \tab \code{qx2000M} \tab numeric \tab qx of completed lifetable functions of Japanese men in 2000.\cr
  \code{[,37]} \tab \code{qx2000F} \tab numeric \tab qx of completed lifetable functions of Japanese women in 2000.\cr
  \code{[,38]} \tab \code{qx2005M} \tab numeric \tab qx of completed lifetable functions of Japanese men in 2005.\cr
  \code{[,39]} \tab \code{qx2005F} \tab numeric \tab qx of completed lifetable functions of Japanese women in 2005.\cr
  \code{[,40]} \tab \code{qx2010M} \tab numeric \tab qx of completed lifetable functions of Japanese men in 2010.\cr
  \code{[,41]} \tab \code{qx2010F} \tab numeric \tab qx of completed lifetable functions of Japanese women in 2010.\cr
  \code{[,42]} \tab \code{qx2015M} \tab numeric \tab qx of completed lifetable functions of Japanese men in 2015.\cr
  \code{[,43]} \tab \code{qx2015F} \tab numeric \tab qx of completed lifetable functions of Japanese women in 2015.\cr
  \code{[,44]} \tab \code{qx2020M} \tab numeric \tab qx of completed lifetable functions of Japanese men in 2020.\cr
  \code{[,45]} \tab \code{qx2020F} \tab numeric \tab qx of completed lifetable functions of Japanese women in 2020.\cr
 }
}
\details{
   qx columns were cited from the completed life tables in Japan for the 1st to 23rd one (7th one was not made, so that it is missing).
   
   \itemize{
    \item \code{Age}: Ages from 0 to 116.
    \item \code{qx1895M}-\code{qx2020M}: qx of 1st to 23rd completed lifetables for Japanese men.
    \item \code{qx1895F}-\code{qx2020F}: qx of 1st to 23rd completed lifetables for Japanese women.
   }
}
\source{
 \url{https://warp.da.ndl.go.jp/info:ndljp/pid/1334623/www.stat.go.jp/english/data/chouki/02.htm}
 \url{https://warp.da.ndl.go.jp/collections/content/info:ndljp/pid/11423429/www.stat.go.jp/data/chouki/zuhyou/02-35.xls}
 \url{https://www.mhlw.go.jp/toukei/saikin/hw/life/20th/index.html}
 \url{https://www.mhlw.go.jp/toukei/saikin/hw/life/21th/index.html}
 \url{https://www.mhlw.go.jp/toukei/saikin/hw/life/22th/index.html}
 \url{https://www.mhlw.go.jp/toukei/saikin/hw/life/23th/index.html}
}
\references{
  Ministry of Health, Labor and Welfare of Japan: Completed lifetables. / Ministry of Internal Affairs and Communications, Statistics Bureau: Historical Statistics of Japan.
}
\keyword{datasets}
