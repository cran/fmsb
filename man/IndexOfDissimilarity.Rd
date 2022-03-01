\name{IndexOfDissimilarity}
\alias{IndexOfDissimilarity}
\title{Index of dissimilarity (ID)}
\description{
 Index of dissimilarity between the 2 age-distributions.
}
\usage{
 IndexOfDissimilarity(X, Y)
}
\arguments{
 \item{X}{ A vector of age-specific standard populations (or percentage) for each age. }
 \item{Y}{ A vector of age-specific target populations (or percentage) for each age. }
}
\value{
 \item{ID}{ Index of dissimilarity, which is a half of sum of absolute differences 
 of percentages for each age, where NA is automatically treated as 0.}
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\references{
 Preston SH, Heuveline P, Guillot M (2001) Demography: Measuring and Modeling 
 Population Processes.  Blackwell Publishing, Oxford.
 
 Newell C (1988) Methods and Models in Demography.  The Guilford Press, New York.
 
 Rowland DT (2003) Demographic methods and concepts. Oxford Univ. Press, Oxford.
}
\examples{
 # Dissimilarities of Japanese population structure were increasing
 # from 1960-1980 (0.132) to 1980-2000 (0.156).
 IndexOfDissimilarity(Jpopl$M1980+Jpopl$F1980, Jpopl$M2000+Jpopl$F2000)
 IndexOfDissimilarity(Jpopl$M1980+Jpopl$F1980, Jpopl$M1960+Jpopl$F1960)
}
