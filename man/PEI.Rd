\name{PEI}
\alias{PEI}
\title{Population Expansion Index}
\description{
 Population Expansion Index (Bulge Index) for movement.
}
\usage{
PEI(X, CLS, MODE)
}
\arguments{
 \item{X}{The vector to give age-specific population from age 0.}
 \item{CLS}{The width of age-class in X. default is 1.}
 \item{MODE}{If the MODE is 1, the ages of 20 to 39 years old are assumed as 
 "easily movable ages" comparing with "relatively unmovable" 10 to 19 and 
 40 to 49, otherwise the ages of 15 to 34 years old are assumed as the former
 and 5 to 14 and 35 to 44 are assumed as the latter, as Dr. Toshio Kuroda 
 suggested in his book. Then PEI (originally named as bulge index, but I
 prefer to use PEI instead) is calculated as the ratio of the population of
 "easily movable ages" to the population of "relatively unmovavle ages" times
 100.  If PEI is larger than 100, net migration may be positive and vise versa.
 Default 1.}
}
\value{
 The value of PEI is returned.
}
\references{
 Kuroda T (1976) Japan's Changing Population Structure (in Japanese). Kokon-Shoin, Tokyo.

 Kuroda T (1971) A study on population composition: Special reference to Japan. 
 (in Japanese, with abstract in English)  Journal of Population Problems
 (Jinko-Mondai-Kenkyu), No. 119: 1-12.
 \url{https://www.ipss.go.jp/syoushika/bunken/data/pdf/j119.pdf}
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\examples{
# Prefectural population estimates in 2018 (unit=1000 persons)
# total of males and females, by 5 year age-class
# (Data source) Download Excel file and extracted
# \url{https://www.e-stat.go.jp/stat-search/file-download?statInfId=000031807147&fileKind=0}
PPT2018 <- data.frame(
 Hokkaido = c(175, 195, 207, 229, 235, 232, 266, 304, 367, 381, 344, 341, 
  354, 452, 368, 310, 252, 274),
 Aomori = c(41, 45, 51, 58, 48, 48, 59, 69, 82, 86, 83, 88, 94, 112, 89, 
  75, 67, 69),
 Iwate = c(41, 47, 52, 57, 47, 50, 59, 69, 81, 82, 78, 84, 91, 106, 82, 
  74, 67, 75),
 Miyagi = c(85, 93, 98, 109, 123, 119, 130, 144, 164, 163, 144, 146, 155, 
  182, 139, 116, 98, 108),
 Akita = c(28, 33, 37, 40, 30, 33, 42, 51, 61, 62, 59, 69, 78, 93, 72, 
  63, 61, 69),
 Yamagata = c(38, 42, 47, 51, 40, 43, 53, 61, 70, 69, 66, 73, 80, 95, 71, 
  62, 58, 73),
 Fukushima = c(67, 70, 79, 89, 73, 80, 94, 105, 122, 124, 117, 129, 139, 
  160, 119, 102, 89, 106),
 Ibaraki = c(106, 117, 127, 140, 132, 130, 152, 172, 204, 216, 185, 176, 
  190, 233, 196, 162, 116, 125),
 Tochigi = c(73, 80, 87, 92, 83, 90, 107, 122, 142, 146, 125, 121, 131, 
  157, 127, 101, 76, 84),
 Gumma = c(70, 79, 87, 97, 90, 88, 99, 113, 139, 148, 127, 117, 124, 153, 
  135, 110, 82, 94),
 Saitama = c(279, 299, 312, 343, 405, 381, 407, 458, 552, 608, 508, 428, 
  415, 524, 488, 411, 275, 235),
 Chiba = c(233, 251, 264, 289, 323, 310, 344, 387, 465, 513, 428, 367, 
  361, 460, 430, 358, 248, 225),
 Tokyo = c(539, 516, 495, 554, 867, 911, 961, 1013, 1109, 1167, 1005, 
  810, 687, 797, 750, 647, 500, 494),
 Kanagawa = c(351, 374, 385, 423, 518, 490, 529, 592, 706, 788, 679, 
  551, 486, 599, 558, 479, 343, 326),
 Niigata = c(78, 88, 94, 103, 92, 94, 111, 127, 152, 154, 141, 141, 
  155, 190, 151, 129, 111, 135),
 Toyama = c(37, 40, 45, 50, 44, 44, 50, 58, 75, 79, 65, 62, 65, 84, 
  80, 63, 49, 60),
 Ishikawa = c(44, 48, 51, 57, 57, 53, 58, 65, 82, 86, 71, 68, 70, 86, 
  80, 63, 47, 58),
 Fukui = c(30, 33, 36, 39, 33, 34, 39, 43, 52, 54, 48, 49, 50, 62, 50, 
  43, 36, 44),
 Yamanashi = c(29, 32, 36, 41, 38, 35, 39, 44, 53, 60, 56, 53, 54, 64, 
  54, 47, 37, 46),
 Nagano = c(76, 85, 94, 101, 79, 83, 98, 114, 143, 150, 133, 127, 130, 
  159, 142, 122, 98, 129),
 Gifu = c(75, 86, 92, 101, 92, 86, 98, 112, 139, 148, 128, 122, 123, 
  155, 139, 116, 90, 96),
 Shizuoka = c(137, 155, 164, 175, 150, 162, 193, 214, 256, 276, 241, 
  225, 230, 284, 251, 212, 161, 173),
 Aichi = c(319, 339, 344, 374, 420, 414, 451, 484, 567, 610, 507, 434, 
  397, 492, 461, 387, 276, 259),
 Mie = c(67, 75, 80, 87, 82, 81, 92, 102, 125, 135, 118, 111, 109, 135, 
  122, 103, 80, 88),
 Shiga = c(61, 67, 69, 74, 74, 71, 79, 88, 105, 108, 89, 83, 81, 100, 
  85, 69, 52, 58),
 Kyoto = c(94, 102, 107, 122, 159, 136, 137, 149, 183, 197, 167, 148, 
  142, 189, 180, 150, 112, 118),
 Osaka = c(334, 352, 370, 418, 486, 461, 489, 529, 645, 727, 608, 507, 
  465, 614, 592, 519, 367, 329),
 Hyogo = c(212, 230, 242, 267, 267, 247, 279, 317, 392, 430, 367, 
  334, 324, 410, 377, 313, 234, 244),
 Nara = c(48, 54, 59, 67, 66, 58, 63, 71, 89, 100, 88, 81, 82, 
  107, 99, 85, 60, 63),
 Wakayama = c(33, 37, 39, 45, 36, 37, 44, 48, 61, 67, 60, 60, 61, 77, 
  69, 59, 47, 54),
 Tottori = c(22, 24, 25, 27, 22, 23, 28, 32, 37, 37, 33, 35, 39, 46, 
  38, 30, 27, 36),
 Shimane = c(26, 28, 29, 32, 25, 27, 32, 37, 43, 43, 38, 42, 46, 58, 
  48, 39, 37, 49),
 Okayama = c(75, 81, 84, 93, 97, 90, 99, 107, 130, 136, 112, 108, 114, 
  142, 133, 106, 86, 104),
 Hiroshima = c(113, 125, 126, 135, 136, 133, 150, 165, 201, 212, 175, 
  163, 167, 211, 194, 155, 119, 137),
 Yamaguchi = c(49, 55, 58, 63, 57, 55, 64, 72, 89, 94, 80, 80, 90, 
  119, 104, 88, 71, 83),
 Tokushima = c(26, 28, 29, 33, 29, 30, 36, 41, 48, 50, 45, 46, 52, 
  64, 53, 43, 37, 46),
 Kagawa = c(37, 40, 42, 46, 40, 40, 48, 54, 68, 69, 58, 57, 61, 79, 
  70, 54, 45, 56),
 Ehime = c(49, 55, 58, 63, 52, 54, 65, 74, 91, 93, 83, 84, 90, 115, 
  98, 79, 66, 81),
 Kochi = c(24, 26, 29, 32, 26, 26, 32, 37, 47, 48, 41, 44, 48, 61, 
  55, 44, 37, 49),
 Fukuoka = c(218, 229, 226, 241, 275, 251, 285, 321, 363, 365, 311, 
  300, 315, 394, 320, 262, 205, 228),
 Saga = c(34, 38, 39, 42, 35, 35, 42, 47, 53, 52, 48, 52, 57, 67, 
  51, 43, 37, 46),
 Nagasaki = c(53, 57, 60, 63, 51, 54, 64, 72, 83, 87, 82, 89, 98, 117, 
  90, 77, 66, 79),
 Kumamoto = c(74, 80, 81, 83, 75, 77, 91, 100, 112, 109, 104, 111, 
  122, 142, 111, 94, 85, 105),
 Oita = c(44, 48, 49, 54, 47, 46, 56, 64, 75, 75, 66, 70, 78, 97, 81, 
  67, 57, 69),
 Miyazaki = c(45, 50, 50, 51, 41, 42, 52, 61, 70, 67, 62, 69, 78, 
  92, 72, 60, 54, 64),
 Kagoshima = c(68, 74, 74, 74, 62, 65, 80, 91, 100, 97, 96, 106, 
  123, 135, 101, 88, 80, 102),
 Okinawa = c(82, 84, 81, 81, 72, 76, 88, 93, 104, 102, 89, 90, 93, 
  97, 60, 56, 48, 52)
)
# Calculate PEI for all prefectures
# for (i in 1:47) {
#  print(PEI(PPT2018[, i], CLS=5))
# }
print(apply(PPT2018, 2, PEI, CLS=5))
}
