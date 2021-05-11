\name{Jvital2013byPref}
\docType{data}
\alias{Jvital2013byPref}
\title{Cross sectional data of several vital statistics in Japan 2013 for each prefecture}
\description{
 The data gives cross sectional data of several vital statistics in Japan 2013
 for each prefecture.
 Included indices are crude birth rates, crude death rates, infant mortality
 rates, and so on.
}
\usage{ Jvital2013byPref }
\format{ A data frame with 47 observations on 34 variables.

 \tabular{rlll}{
  \code{[, 1]} \tab \code{PNAME} \tab factor w/47 levels \tab The name (in roma-ji) for prefectures\cr
  \code{[, 2]} \tab \code{JCODE} \tab numeric \tab Prefecture number defined by Geographical Information Authority of Japan\cr
  \code{[, 3]} \tab \code{CBR} \tab numeric \tab Crude birth rates\cr
  \code{[, 4]} \tab \code{CDR} \tab numeric \tab Crude death rates\cr
  \code{[, 5]} \tab \code{IMR} \tab numeric \tab Infant mortality rates\cr
  \code{[, 6]} \tab \code{NMR} \tab numeric \tab Neonatal mortality rates\cr
  \code{[, 7]} \tab \code{NIR} \tab numeric \tab Natural increase rates\cr
  \code{[, 8]} \tab \code{SBRPB} \tab numeric \tab Stillbirth rates\cr
  \code{[, 9]} \tab \code{SARPB} \tab numeric \tab Spontaneous abortion rates\cr
  \code{[,10]} \tab \code{ACRPB} \tab numeric \tab Artificial contraception rates\cr
  \code{[,11]} \tab \code{PNMPB} \tab numeric \tab Perinatal mortalities per birth\cr
  \code{[,12]} \tab \code{SBRA22W} \tab numeric \tab Stillbirth rate after gestational
   age of 22 weeks per birth\cr
  \code{[,13]} \tab \code{ENMR} \tab numeric \tab Early neonatal mortality rate per live birth\cr
  \code{[,14]} \tab \code{MR} \tab numeric \tab Marriage rates\cr
  \code{[,15]} \tab \code{DR} \tab numeric \tab Divorce rates\cr
  \code{[,16]} \tab \code{TFR} \tab numeric \tab Total fertility rates\cr
  \code{[,17]} \tab \code{CSM.ALL} \tab numeric \tab Cause-specific mortality for all causes\cr
  \code{[,18]} \tab \code{CSM.CANCER} \tab numeric \tab Cause-specific mortality for cancer\cr
  \code{[,19]} \tab \code{CSM.HD} \tab numeric \tab Cause-specific mortality for heart disease
   except for hypertention\cr
  \code{[,20]} \tab \code{CSM.PNEUM} \tab numeric \tab Cause-specific mortality for pneumonia\cr
  \code{[,21]} \tab \code{CSM.CEVD} \tab numeric \tab Cause-specific mortality for 
   cerebrovascular disease\cr
  \code{[,22]} \tab \code{CSM.SEN} \tab numeric \tab Cause-specific mortality for senescence\cr
  \code{[,23]} \tab \code{CSM.ACC} \tab numeric \tab Cause-specific mortality for accidents\cr
  \code{[,24]} \tab \code{CSM.SUI} \tab numeric \tab Cause-specific mortality for suicide\cr
  \code{[,25]} \tab \code{CSM.KF} \tab numeric \tab Cause-specific mortality for kidney failure\cr
  \code{[,26]} \tab \code{CSM.COPD} \tab numeric \tab Cause-specific mortality for chronic
   obstructive pulmonary disease\cr
  \code{[,27]} \tab \code{CSM.AA} \tab numeric \tab Cause-specific mortality for aneuysm of aorta\cr
  \code{[,28]} \tab \code{CSM.LIVD} \tab numeric \tab Cause-specific mortality for liver disease\cr
  \code{[,29]} \tab \code{CSM.DIAB} \tab numeric \tab Cause-specific mortality for diabetes\cr
  \code{[,30]} \tab \code{CSM.SEP} \tab numeric \tab Cause-specific mortality for sepsis\cr
  \code{[,31]} \tab \code{CSM.MNP} \tab numeric \tab Cause-specific mortality for miscellaneous 
   neoplasms\cr
  \code{[,32]} \tab \code{CSM.DEM} \tab numeric \tab Cause-specific mortality for dementia\cr
  \code{[,33]} \tab \code{CSM.TB} \tab numeric \tab Cause-specific mortality for tuberculosis\cr
  \code{[,34]} \tab \code{CSM.TA} \tab numeric \tab Cause-specific mortality for traffic accidents\cr
 }
}
\details{
   Official vital statistics in Japan in 2013 for each prefecture.
   
   \itemize{
    \item \code{PNAME}: The name (in roma-ji) for prefectures.
    \item \code{JCODE}: Prefecture number defined by Geographical Information Authority of Japan. From 1 to 47.
    \item \code{CBR}: Crude birth rate.  Number of all live birth / mid-year population 1000.
    \item \code{CDR}: Crude death rate.  Number of death / mid-year population 1000.
    \item \code{IMR}: Infant mortality rate. Number of death at age 0 / 1000 live births.
    \item \code{NMR}: Neonatal mortality rate. Number of death within 4 weeks after birth / 1000 live births.
    \item \code{NIR}: Natural increase rate.  \code{CBR}-\code{CDR}.
    \item \code{SBRPB}: Stillbirth rate per birth.  Number of stillbirths / 1000 births.
    \item \code{SARPB}: Spontaneous abortion rate per birth.  Number of spontaneous abortions / 1000 births.
    \item \code{ACRPB}: Artificial contraception (= induced abortion) rate per birth.  Number of induced abortions / 1000 births.
    \item \code{PNMPB}: Perinatal mortality per birth.  [(Number of stillbirths after gestational age 22 weeks) + (Number of early neonatal deaths within a week after birth)] per 1000 births.  The denominator is the sum of the number of stillbirths after gestational age 22 weeks and the number of live births.  This definition was established in 1995, but \code{PNMPB} also includes some values before 1995.
    \item \code{SBRA22W}: Stillbirth rate after gestational age of 22 weeks per 1000 births.
    \item \code{ENMR}: Early neonatal mortality rate per 1000 live births.
    \item \code{MR}: Marital rate.  The number of marriages / mid-year population 1000.
    \item \code{DR}: Divorce rate.  The number of divorces / mid-year population 1000.
    \item \code{TFR}: Total fertility rate.  The sum of age-specific fertility rates, which is the number of births divided by the number of women's population for each age.
    \item \code{CSM.ALL}: Cause-specific mortality for all causes.  Similar to \code{CDR}, but the denominator is mid-year population 100000 instead of 1000.
    \item \code{CSM.CANCER}: Cause-specific mortality for cancer.  The number of deaths caused by cancer / mid-year population 100000.
    \item \code{CSM.HD}: Cause-specific mortality for heart disease except for hypertention / mid-year population 100000.
    \item \code{CSM.PNEUM}: Cause-specific mortality for pneumonia / mid-year population 100000.
    \item \code{CSM.CEVD}: Cause-specific mortality for cerebrovascular disease / mid-year population 100000.
    \item \code{CSM.SEN}: Cause-specific mortality for senescence / mid-year population 100000.
    \item \code{CSM.ACC}: Cause-specific mortality for accidents / mid-year population 100000.
    \item \code{CSM.SUI}: Cause-specific mortality for suicide / mid-year population 100000.
    \item \code{CSM.KF}: Cause-specific mortality for kidney failure / mid-year population 100000.
    \item \code{CSM.COPD}: Cause-specific mortality for chronic obstructive pulmonary disease / mid-year population 100000.
    \item \code{CSM.AA}: Cause-specific mortality for aneuysm of aorta / mid-year population 100000.
    \item \code{CSM.LIVD}: Cause-specific mortality for liver disease / mid-year population 100000.
    \item \code{CSM.DIAB}: Cause-specific mortality for diabetes / mid-year population 100000.
    \item \code{CSM.SEP}: Cause-specific mortality for sepsis / mid-year population 100000.
    \item \code{CSM.MNP}: Cause-specific mortality for miscellaneous neoplasms / mid-year population 100000.
    \item \code{CSM.DEM}: Cause-specific mortality for dementia / mid-year population 100000.
    \item \code{CSM.TB}: Cause-specific mortality for tuberculosis / mid-year population 100000.
    \item \code{CSM.TA}: Cause-specific mortality for traffic accidents / mid-year population 100000.
   }
}
\source{
 \url{https://www.mhlw.go.jp/toukei/saikin/hw/jinkou/kakutei13/xls/hyo.xls}
 \url{https://www.mhlw.go.jp/toukei/saikin/hw/jinkou/kakutei13/xls/sankou.xls}
}
\references{
  Ministry of Health, Labor and Welfare of Japan: Vital Statistics 2013.
}
\keyword{datasets}
