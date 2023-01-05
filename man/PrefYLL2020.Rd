\name{PrefYLL2020}
\docType{data}
\alias{PrefYLL2020}
\title{Years of Life Lost by several causes in Japan 2020 for each prefecture}
\description{
 The data gives years of life lost by several causes in Japan 2020 for each prefecture.
 There are several definitions of YLL.  For example, WHO's Global Burden of Disease
 defines the YLL as  the number of deaths multiplied by the standard life expectancy 
 at the age at which death occurs, for a given cause, age and sex (WHO).
 However, Japanese Ministry of Health, Labor and Welfare gives the expected increase 
 of the life expectancy at birth if the mortality due to each cause of death is 
 removed from the age-specific mortality as the measure of YLL, and thus this
 dataset implements such data derived from the report of regional life tables
 in Japan (Ministry of Health, Labor and Welfare, 2020).
 Until 2015, deaths caused by tuberculosis were analyzed, but in 2020, deaths
 caused by COVID-19 are calculated instead.
}
\usage{PrefYLL2020}
\format{ A data frame with 47 observations on 28 variables.

 \tabular{rlll}{
  \code{[, 1]} \tab \code{PNAME} \tab factor w/47 levels \tab The name (in roma-ji) for prefectures\cr
  \code{[, 2]} \tab \code{JCODE} \tab numeric \tab Prefecture number defined by Geographical Information 
   Authority of Japan\cr
  \code{[, 3]} \tab \code{CancerM} \tab numeric \tab Years of Life Lost (YLL) of males by cancer\cr
  \code{[, 4]} \tab \code{CardioM} \tab numeric \tab Years of Life Lost (YLL) of males by heart diseases
   except for hypertention\cr
  \code{[, 5]} \tab \code{CerebroM} \tab numeric \tab Years of Life Lost (YLL) of males by 
   cerebrovascular disease\cr
  \code{[, 6]} \tab \code{Top3M} \tab numeric \tab Years of Life Lost (YLL) of males by cancer, heart 
   disease or cerebrovascular disease\cr
  \code{[, 7]} \tab \code{PneumoniaM} \tab numeric \tab Years of Life Lost (YLL) of males by pneumonia\cr
  \code{[, 8]} \tab \code{AccidentM} \tab numeric \tab Years of Life Lost (YLL) of males by accident\cr
  \code{[, 9]} \tab \code{TrafficM} \tab numeric \tab Years of Life Lost (YLL) of males by 
   traffic accidents\cr
  \code{[,10]} \tab \code{SuicideM} \tab numeric \tab Years of Life Lost (YLL) of males by suicide\cr
  \code{[,11]} \tab \code{KidneyM} \tab numeric \tab Years of Life Lost (YLL) of males by kidney failure\cr
  \code{[,12]} \tab \code{LiverM} \tab numeric \tab Years of Life Lost (YLL) of males by liver disease\cr
  \code{[,13]} \tab \code{DiabetesM} \tab numeric \tab Years of Life Lost (YLL) of males by diabetes\cr
  \code{[,14]} \tab \code{HypertensM} \tab numeric \tab Years of Life Lost (YLL) of males by hypertension\cr
  \code{[,15]} \tab \code{Covid19M} \tab numeric \tab Years of Life Lost (YLL) of males by tuberculosis\cr
  \code{[,16]} \tab \code{CancerF} \tab numeric \tab Years of Life Lost (YLL) of females by cancer\cr
  \code{[,17]} \tab \code{CardioF} \tab numeric \tab Years of Life Lost (YLL) of females by heart 
   diseases except for hypertention\cr
  \code{[,18]} \tab \code{CerebroF} \tab numeric \tab Years of Life Lost (YLL) of females by 
   cerebrovascular disease\cr
  \code{[,19]} \tab \code{Top3F} \tab numeric \tab Years of Life Lost (YLL) of females by cancer, heart 
   disease or cerebrovascular disease\cr
  \code{[,20]} \tab \code{PneumoniaF} \tab numeric \tab Years of Life Lost (YLL) of females by pneumonia\cr
  \code{[,21]} \tab \code{AccidentF} \tab numeric \tab Years of Life Lost (YLL) of females by accident\cr
  \code{[,22]} \tab \code{TrafficF} \tab numeric \tab Years of Life Lost (YLL) of females by 
   traffic accidents\cr
  \code{[,23]} \tab \code{SuicideF} \tab numeric \tab Years of Life Lost (YLL) of females by suicide\cr
  \code{[,24]} \tab \code{KidneyF} \tab numeric \tab Years of Life Lost (YLL) of females by kidney failure\cr
  \code{[,25]} \tab \code{LiverF} \tab numeric \tab Years of Life Lost (YLL) of females by liver disease\cr
  \code{[,26]} \tab \code{DiabetesF} \tab numeric \tab Years of Life Lost (YLL) of females by diabetes\cr
  \code{[,27]} \tab \code{HypertensF} \tab numeric \tab Years of Life Lost (YLL) of females by hypertension\cr
  \code{[,28]} \tab \code{Covid19F} \tab numeric \tab Years of Life Lost (YLL) of females by tuberculosis\cr
 }
}
\details{
   Years of Life Lost by several causes in Japan 2020 for each prefecture.
   
   \itemize{
    \item \code{PNAME}: The name (in roma-ji) for prefectures.
    \item \code{JCODE}: Prefecture number defined by Geographical Information Authority of Japan. 
     From 1 to 47.
    \item \code{Cancer[M|F]}: YLL by cancer for males ([M]) or for females ([F]).
    \item \code{Cardio[M|F]}: YLL by heart disease for males ([M]) or for females ([F]).
    \item \code{Cerebro[M|F]}: YLL by cerebrovascular disease for males ([M]) or for females ([F]).
    \item \code{Top3[M|F]}: YLL by above 3 major diseases for males ([M]) or for females ([F]).
    \item \code{Peumonia[M|F]}: YLL by pneumonia for males ([M]) or for females ([F]).
    \item \code{Accident[M|F]}: YLL by accidents for males ([M]) or for females ([F]).
    \item \code{Traffic[M|F]}: YLL by traffic accidents (it's also included in \code{Accident[M|F]}
     for males ([M]) or for females ([F]).
    \item \code{Suicide[M|F]}: YLL by suicide for males ([M]) or for females ([F]).
    \item \code{Kidney[M|F]}: YLL by kidney failure for males ([M]) or for females ([F]).
    \item \code{Liver[M|F]}: YLL by liver disease for males ([M]) or for females ([F]).
    \item \code{Diabetes[M|F]}: YLL by diabates for males ([M]) or for females ([F]).
    \item \code{Hypertension[M|F]}: YLL by hypertension for males ([M]) or for females ([F]).
    \item \code{Covid19[M|F]}: YLL by COVID-19 for males ([M]) or for females ([F]).
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
x <- PrefYLL2020
COL <- ifelse(x$PNAME=="Nagano", "blue", ifelse(x$PNAME=="Okinawa", "red", 
              ifelse(x$PNAME=="Shiga", "green", "lightgrey")))
LWD <- ifelse(x$PNAME=="Nagano", 2, ifelse(x$PNAME=="Okinawa", 2, ifelse(x$PNAME=="Shiga", 2, 1)))
LTY <- ifelse(x$PNAME=="Nagano", 1, ifelse(x$PNAME=="Okinawa", 1, ifelse(x$PNAME=="Shiga", 1, 3)))
VX <- c("Cancer","Heart\n Disease","Cerebrovascular\n Disease","Top 3 causes","Pneumonia",
 "Accident","(Traffic\n Accident)","Suicide","Kidney\n Failure","Liver\n Disease",
 "Diabetes","Hypertension","Tuberculosis")
males <- x[,3:15]
females <- x[,16:28]
layout(t(1:2))
radarchart(males, maxmin=FALSE, pcol=COL, axistype=2, pty=32, plty=LTY, plwd=LWD, vlabels=VX,
 title="YLLs in males (2020)\n (Blue: Nagano, Green: Shiga,\n Red: Okinawa, Gray: Others)")
radarchart(females, maxmin=FALSE, pcol=COL, axistype=2, pty=32, plty=LTY, plwd=LWD, vlabels=VX,
 title="YLL in females (2020)\n (Blue: Nagano, Green: Shiga,\n Red: Okinawa, Gray: Others)")}
\keyword{datasets}
