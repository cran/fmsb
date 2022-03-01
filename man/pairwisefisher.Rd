\name{pairwise.fisher.test}
\alias{pairwise.fisher.test}
\title{Exact version of pairwise.prop.test}
\description{
 By conducting repeatedly Fisher's exact tests instead of chi-square tests, this function can test the null-hypothesis of no difference in any pair of proportions for more than 2 groups, with adjustment of type I error for multiple comparison.
}
\usage{
 pairwise.fisher.test(x, n, p.adjust.method, ...)
}
\arguments{
 \item{x}{A integer vector of event occurences}
 \item{n}{A integer vector of trials}
 \item{p.adjust.method}{A name in p.adjust.methods to specify the method to adjust type I error for multiple comparisons.  Default is "holm".}
 \item{...}{Miscellaneous arguments to be given for fisher.test().}
}
\value{
 An object of adjusted p-values for all possible comparisons of pairs with class pairwise.htest.
}
\seealso{pairwise.prop.test, p.adjust.methods}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}.  The code of this function was provided by Dr. Shigenobu AOKI (Gunma Univ.).}
\examples{
 pairwise.fisher.test(c(2, 4, 5), c(10, 14, 17), p.adjust.method="bonferroni")
 smoker <- c(2, 1, 7)
 total <- c(11, 14, 10)
 names(total) <- c("A", "B", "C")
 pairwise.fisher.test(smoker, total)
}
