\name{gstem}
\alias{gstem}
\title{Draw the stem and leaf plot in a graphic device.}
\description{
 Capture the output of stem() function and plot them into graphic devices.
 However, the result of setting scale parameter as 2 may be controversial.
}
\usage{ gstem(X, scale) }
\arguments{
 \item{X}{A numeric vector.}
 \item{scale}{Parameter to control plot length of graph.  Default is 1.}
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\examples{
 x <- rnorm(100, 10, 1)
 stem(x)
 stem(x, 2)
 layout(t(1:2))
 gstem(x)
 gstem(x, 2)
}
