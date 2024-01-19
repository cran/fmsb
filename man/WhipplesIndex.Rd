\name{WhipplesIndex}
\alias{WhipplesIndex}
\title{Whipple's Index}
\description{
 Whipple's Index for age-heaping
}
\usage{
 WhipplesIndex(X)
}
\arguments{
 \item{X}{ The integer vector to give age-specific population from age 0 to more than 63 for each age. }
}
\value{
 \item{WI}{ The Whipple's Index. }
 \item{JUDGE}{ Based on Whipple's Index, accuracy of age-reporting is judged. }
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\references{
 Preston SH, Heuveline P, Guillot M (2001) Demography: Measuring and Modeling 
 Population Processes.  Blackwell Publishing, Oxford.
 
 Newell C (1988) Methods and Models in Demography.  The Guilford Press, New York.
 
 Rowland DT (2003) Demographic methods and concepts. Oxford Univ. Press, Oxford.

 Ministry of Home Affairs, India (2011) 2011 Census C-13. \url{https://censusindia.gov.in/nada/index.php/catalog/22542/download/25673/PC01_C13_00.xls}
}
\examples{
 WhipplesIndex(Jpop$M2000)
# India <- read.delim("https://minato.sip21c.org/ldaR/India2011census.txt")
# CRAN requires the example code can work without internet connection
# since 2024, and thus I modified the code to include the data here.
 India <- data.frame(
  Males = c(10633298, 11381468, 11952853, 12331431, 12333024, 13725480, 
  13394700, 12903364, 14061937, 12214985, 16089436, 12962604, 14637892, 
  12563775, 13165128, 13739746, 13027935, 11349449, 15020851, 10844415, 
  14892165, 10532278, 12392976, 9674189, 10093085, 14311524, 10315030, 
  8552032, 10719926, 7445696, 15628996, 7157502, 8801105, 6108879, 
  6964192, 15036666, 8067568, 5784879, 8090401, 5939867, 15173411, 
  6172297, 6856826, 4468914, 4873938, 12685175, 5735540, 4043122, 
  5568554, 4105723, 11379329, 4323584, 4068700, 2808043, 3263610, 
  7769352, 3666804, 2339391, 3072508, 2607957, 8677046, 3095448, 
  2892015, 1977207, 2060033, 6275854, 2278670, 1353711, 1640034, 
  1396057, 5393714, 1584873, 1176727, 708381, 787804, 2278704, 832251, 
  438394, 506957, 434297, 1725200, 491522, 306378, 192946, 210994, 
  580527, 215850, 112348, 112374, 99007, 360237, 118606, 75430, 46220, 
  51972, 124950, 57894, 35238, 48393, 28284, 289325),
 Females = c(9677936, 10373729, 11103415, 11642610, 11377014, 12328750, 
 12259545, 11923276, 12906436, 11209653, 14462671, 11778342, 13239415, 
 11716908, 12093041, 12159708, 11564358, 9868018, 12937296, 10014673, 
 13990570, 9446694, 11135249, 9479866, 9787150, 13456554, 9761967, 
 8157318, 11407090, 7286828, 14770033, 6665743, 8812439, 6655662, 
 7030400, 13385965, 7760149, 5907352, 9381357, 5786480, 13355581, 
 5395597, 6523816, 4865438, 4752294, 11187786, 5257138, 3908175, 
 6081038, 3746076, 10083093, 3562382, 3666464, 2782747, 3131302, 
 7838194, 3405033, 2259635, 3646426, 2540755, 9133643, 2931365, 
 2853128, 2016898, 2026924, 6746498, 2233276, 1251371, 1908339, 
 1371173, 5592566, 1499310, 1074202, 658155, 733110, 2493642, 834882, 
 396654, 561458, 455264, 2059738, 536294, 297415, 187239, 212503, 
 684271, 232048, 109063, 123266, 114413, 472835, 139691, 78131, 
 48410, 55002, 147584, 62374, 36175, 56118, 36287, 316453),
 Age = 0:100)
 WhipplesIndex(India$Males)
# To check age-heaping graphically,
# you can install and load pyramid package from cran
# and pyramid(India, Cstep=5) may be useful.
}
