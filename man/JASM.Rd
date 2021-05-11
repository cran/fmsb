\name{JASM}
\docType{data}
\alias{JASM}
\title{Japanese mortality data}
\description{
 The data gives the sex and age-class (by five) specific numbers of death
 in Showa 60 (S60 = 1985), Heisei 2 (H02 = 1990), Heisei 7 (H07 = 1995),
 Heisei 12 (H12 = 2000), Heisei 17 (H17 = 2005), Heisei 22 (H22 = 2010)
 and Heisei 27 (H27 = 2015), and corresponding populations.
}
\usage{ JASM }
\format{ A data frame with 18 observations on 30 variables.

 \tabular{rlll}{
  \code{[, 1]} \tab \code{AGECLASS} \tab Factor w/18 levels \tab Age class (years old)\cr
  \code{[, 2]} \tab \code{S60MODEL} \tab numeric \tab Model population in 1985\cr
  \code{[, 3]} \tab \code{S60M} \tab numeric \tab Number of males' death by age classes in 1985\cr
  \code{[, 4]} \tab \code{H02M} \tab numeric \tab Number of males' death by age classes in 1990\cr
  \code{[, 5]} \tab \code{H07M} \tab numeric \tab Number of males' death by age classes in 1995\cr
  \code{[, 6]} \tab \code{H12M} \tab numeric \tab Number of males' death by age classes in 2000\cr
  \code{[, 7]} \tab \code{H17M} \tab numeric \tab Number of males' death by age classes in 2005\cr
  \code{[, 8]} \tab \code{H22M} \tab numeric \tab Number of males' death by age classes in 2010\cr
  \code{[, 8]} \tab \code{H27M} \tab numeric \tab Number of males' death by age classes in 2015\cr
  \code{[, 9]} \tab \code{S60F} \tab numeric \tab Number of females' death by age classes in 1985\cr
  \code{[,10]} \tab \code{H02F} \tab numeric \tab Number of females' death by age classes in 1990\cr
  \code{[,11]} \tab \code{H07F} \tab numeric \tab Number of females' death by age classes in 1995\cr
  \code{[,12]} \tab \code{H12F} \tab numeric \tab Number of females' death by age classes in 2000\cr
  \code{[,13]} \tab \code{H17F} \tab numeric \tab Number of females' death by age classes in 2005\cr
  \code{[,14]} \tab \code{H22F} \tab numeric \tab Number of females' death by age classes in 2010\cr
  \code{[,14]} \tab \code{H27F} \tab numeric \tab Number of females' death by age classes in 2015\cr
  \code{[,15]} \tab \code{S60MP} \tab numeric \tab Number of males' population by age classes in 1985\cr
  \code{[,16]} \tab \code{H02MP} \tab numeric \tab Number of males' population by age classes in 1990\cr
  \code{[,17]} \tab \code{H07MP} \tab numeric \tab Number of males' population by age classes in 1995\cr
  \code{[,18]} \tab \code{H12MP} \tab numeric \tab Number of males' population by age classes in 2000\cr
  \code{[,19]} \tab \code{H17MP} \tab numeric \tab Number of males' population by age classes in 2005\cr
  \code{[,20]} \tab \code{H22MP} \tab numeric \tab Number of males' population by age classes in 2010\cr
  \code{[,20]} \tab \code{H27MP} \tab numeric \tab Number of males' population by age classes in 2015\cr
  \code{[,21]} \tab \code{S60FP} \tab numeric \tab Number of females' population by age classes in 1985\cr
  \code{[,22]} \tab \code{H02FP} \tab numeric \tab Number of females' population by age classes in 1990\cr
  \code{[,23]} \tab \code{H07FP} \tab numeric \tab Number of females' population by age classes in 1995\cr
  \code{[,24]} \tab \code{H12FP} \tab numeric \tab Number of females' population by age classes in 2000\cr
  \code{[,25]} \tab \code{H17FP} \tab numeric \tab Number of females' population by age classes in 2005\cr
  \code{[,26]} \tab \code{H22FP} \tab numeric \tab Number of females' population by age classes in 2010\cr
  \code{[,26]} \tab \code{H27FP} \tab numeric \tab Number of females' population by age classes in 2015\cr
 }
}
\details{
   Japanese mortality data by sex and age-class (by five) given as national
   official vital statitistics from 1985 to 2015, every 5 years.
   
   \itemize{
    \item \code{AGECLASS}: Labels for age classes.  [0-4] to [85-].
    \item \code{S60MODEL}: Age class specific model population of Japan in 1985.
    \item \code{S60M}-\code{H27M}: Age class specific number of death of males in 1985-2015.
    \item \code{S60F}-\code{H27F}: Age class specific number of death of females in 1985-2015.
    \item \code{S60MP}-\code{H27MP}: Age class specific number of males' population in 1985-2015.
    \item \code{S60FP}-\code{H27FP}: Age class specific number of females' population in 1985-2015.
   }
}
\source{
 \url{https://www.stat.go.jp/english/data/nenkan/66nenkan/index.html}
}
\references{
  Ministry of Health, Labor and Welfare of Japan: Vital Statistics.
}
\keyword{datasets}
