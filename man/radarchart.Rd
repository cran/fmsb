\name{radarchart}
\alias{radarchart}
\alias{radarchartcirc}
\title{Drawing radar chart (a.k.a. spider plot)}
\description{
 Drawing the radar chart with several lines from a data frame, 
 which must be composed of more than 3 variables as axes and the 
 rows indicate cases as series.
 The radatchart() uses the polygons as radar grid, radarchartcirc()
 uses circles as radar grid.
}
\usage{
 radarchart(df, axistype, seg, pty, pcol, plty, plwd, pdensity, pangle, pfcol, 
 cglty, cglwd, cglcol, axislabcol, title, maxmin, na.itp, centerzero, 
 vlabels, vlcex, caxislabels, calcex, paxislabels, palcex, ...)

 radarchartcirc(df, axistype, seg, pty, pcol, plty, plwd, pdensity, pangle, pfcol, 
 cglty, cglwd, cglcol, axislabcol, title, maxmin, na.itp, centerzero, 
 vlabels, vlcex, caxislabels, calcex, paxislabels, palcex, ...)
}
\arguments{
 \item{df}{The data frame to be used to draw radarchart.  If maxmin is TRUE, this must include maximum values as row 1 and minimum values as row 2 for each variables, and actual data should be given as row 3 and lower rows.  The number of columns (variables) must be more than 2.}
 \item{axistype}{The type of axes, specified by any of 0:5.  0 means no axis label.  1 means center axis label only.  2 means around-the-chart label only.  3 means both center and around-the-chart (peripheral) labels.  4 is *.** format of 1, 5 is *.** format of 3.  Default is 0.}
 \item{seg}{The number of segments for each axis (default 4).}
 \item{pty}{A vector to specify point symbol: Default 16 (closed circle), if you don't plot data points, it should be 32.  This is repeatedly used for data series.}
 \item{pcol}{A vector of color codes for plot data: Default 1:8, which are repeatedly used.}
 \item{plty}{A vector of line types for plot data: Default 1:6, which are repeatedly used.}
 \item{plwd}{A vector of line widths for plot data: Default 1, which is repeatedly used.}
 \item{pdensity}{A vector of filling density of polygons: Default NULL, which is repeatedly used.}
 \item{pangle}{A vector of the angles of lines used as filling polygons: Default 45, which is repeatedly used.}
 \item{pfcol}{A vector of color codes for filling polygons: Default NA, which is repeatedly usd.}
 \item{cglty}{Line type for radar grids: Default 3, dotted line.  For radatchartcirc(), default 1, solid line.}
 \item{cglwd}{Line width for radar grids: Default 1, which means thinnest line.}
 \item{cglcol}{Line color for radar grids: Default "navy"}
 \item{axislabcol}{Color of axis label and numbers: Default "blue"}
 \item{title}{if any, title should be typed.}
 \item{maxmin}{Logical.  If true, data frame includes possible maximum values as row 1 and possible minimum values as row 2.  If false, the maximum and minimum values for each axis will be calculated as actual maximum and minimum of the data.  Default TRUE.}
 \item{na.itp}{Logical. If true, items with NA values are interpolated from nearest neighbor items and connect them.  If false, items with NA are treated as the origin (but not pointed, only connected with lines).  Default FALSE.}
 \item{centerzero}{Logical. If true, this function draws charts with scaling originated from (0,0).  If false, charts originated from (1/segments).  Default FALSE.}
 \item{vlabels}{Character vector for the names for variables.  If NULL, the names of the variables as colnames(df) are used.  Default NULL.}
 \item{vlcex}{Font size magnification for vlabels.  If NULL, the font size is fixed at text()'s default.  Default NULL.}
 \item{caxislabels}{Character vector for center axis labels, overwriting values specified in axistype option.  If NULL, the values specified by axistype option are used. Default is NULL.}
 \item{calcex}{Font size magnification for caxislabels.  If NULL, the font size is fixed at text()'s default.  Default NULL.}
 \item{paxislabels}{Character vector for around-the-chart (peripheral) labels, overwriting values specified in axistype option.  If NULL, the values specified by axistype option are used. Default is NULL.}
 \item{palcex}{Font size magnification for paxislabels.  If NULL, the font size is fixed at text()'s default.  Default NULL.}
 \item{...}{Miscellaneous arguments to be given for plot.default().}
}
\value{
 No value is returned.
}
\author{Minato Nakazawa \email{minatonakazawa@gmail.com} \url{https://minato.sip21c.org/}}
\examples{
# Data must be given as the data frame, where the first cases show maximum.
maxmin <- data.frame(
 total=c(5, 1),
 phys=c(15, 3),
 psycho=c(3, 0),
 social=c(5, 1),
 env=c(5, 1))
# data for radarchart function version 1 series, minimum value must be omitted from above.
RNGkind("Mersenne-Twister")
set.seed(123)
dat <- data.frame(
 total=runif(3, 1, 5),
 phys=rnorm(3, 10, 2),
 psycho=c(0.5, NA, 3),
 social=runif(3, 1, 5),
 env=c(5, 2.5, 4))
dat <- rbind(maxmin, dat)
VARNAMES <- c("Total\nQOL", "Physical\naspects", "Phychological\naspects", 
 "Social\naspects", "Environmental\naspects")
op <- par(mar=c(1, 2, 2, 1), mfrow=c(2, 3))
radarchart(dat, axistype=1, seg=5, plty=1, vlabels=VARNAMES, 
 title="(axis=1, 5 segments, with specified vlabels)", vlcex=0.5)
radarchart(dat, axistype=2, pcol=topo.colors(3), plty=1, pdensity=c(5, 10, 30), 
 pangle=c(10, 45, 120), pfcol=topo.colors(3), 
 title="(topo.colors, fill with hatching, axis=2)")
radarchart(dat, axistype=2, pcol=topo.colors(3), plty=1, 
 pfcol=adjustcolor(topo.colors(3), 0.3), 
 title="(topo.colors, fill with transparency, axis=2)")
radarchart(dat, axistype=3, pty=32, plty=1, axislabcol="grey", na.itp=FALSE,
 title="(no points, axis=3, na.itp=FALSE)")
radarchartcirc(dat, axistype=3, pty=32, plty=1, axislabcol="grey", na.itp=FALSE,
 title="(no points, axis=3, na.itp=FALSE, circular radar grid)")
radarchart(dat, axistype=1, plwd=1:5, pcol=1, centerzero=TRUE, 
 seg=4, caxislabels=c("worst", "", "", "", "best"),
 title="(use lty and lwd but b/w, axis=1,\n centerzero=TRUE, with centerlabels)")
par(op)
}
