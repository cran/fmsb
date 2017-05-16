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
\author{Minato Nakazawa \email{minato-nakazawa@umin.net} \url{http://minato.sip21c.org/}}
\references{
 Preston SH, Heuveline P, Guillot M (2001) Demography: Measuring and Modeling 
 Population Processes.  Blackwell Publishing, Oxford.
 
 Newell C (1988) Methods and Models in Demography.  The Guilford Press, New York.
 
 Rowland DT (2003) Demographic methods and concepts. Oxford Univ. Press, Oxford.
}
\examples{
 WhipplesIndex(Jpop$M2000)
}
