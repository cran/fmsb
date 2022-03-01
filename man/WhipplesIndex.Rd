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

 Ministry of Home Affairs, India (2011) 2011 Census C-13. \url{https://censusindia.gov.in/2011census/C-series/C-13.html}
}
\examples{
 WhipplesIndex(Jpop$M2000)
 India <- read.delim("https://minato.sip21c.org/ldaR/India2011census.txt")
 WhipplesIndex(India$Males)
}
