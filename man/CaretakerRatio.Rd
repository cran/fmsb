\name{CaretakerRatio}
\alias{CaretakerRatio}
\title{Caretaker Ratio}
\description{
 Caretaker ratio.  Defined as the ratio of the aged population who may need care 
 to caretaking females population.
}
\usage{
 CaretakerRatio(PM, PF)
}
\arguments{
 \item{PM}{ The integer vector to give age-specific population from age 0 to more than 80 for males.}
 \item{PF}{ The integer vector to give age-specific population from age 0 to more than 80 for females.}
}
\value{
 \item{CR}{ Caretaker Ratio.}
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\references{
 Preston SH, Heuveline P, Guillot M (2001) Demography: Measuring and Modeling 
 Population Processes.  Blackwell Publishing, Oxford.
 
 Newell C (1988) Methods and Models in Demography.  The Guilford Press, New York.
 
 Rowland DT (2003) Demographic methods and concepts. Oxford Univ. Press, Oxford.
}
\examples{
 # Caretaker Ratio in Japan in 2015.  The value 81.72 is much higher than
 # 46 observed in UK in 1990.
 CaretakerRatio(PM=Jpop$M2015, PF=Jpop$F2015)
}
