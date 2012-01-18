\name{neuroblastoma}
\alias{neuroblastoma}
\docType{data}
\title{
Neuroblastoma copy number profiles and breakpoint annotations
}
\description{
Tumors from patients at the Institut Curie were assayed using aCGH
and SNP microarrays. Their normalized copy number profiles are available
as \code{neuroblastoma$profiles} and the breakpoint annotations are available
as \code{neuroblastoma$annotations}.
}
\usage{data(neuroblastoma)}
\format{
  A named list of 2 data frames:
  \describe{
    \item{\code{profiles}}{A data.frame with
      one row for each probe, and
      these variables:
    \describe{
      \item{\code{profile.id}}{
	factor: id of copy number profile.}
      \item{\code{chromosome}}{
	factor: chromosome on which the probe was mapped.}
      \item{\code{position}}{
	integer: probe was mapped to this position
	in base pairs.}
      \item{\code{logratio}}{
	numeric: normalized logratio of the probe,
	which should be
	proportional to copy number.}
    }}
    \item{\code{annotations}}{a data.frame with
      one row for each annotated region,
      and these variables:
    \describe{
      \item{\code{profile.id}}{
	factor: id of copy number profile.}
      \item{\code{chromosome}}{
	factor: chromosome of this annotation.}
      \item{\code{min}}{
	integer: lower limit of this region in base pairs.}
      \item{\code{max}}{
	integer: upper limit of this region in base pairs.}
      \item{\code{annotation}}{
	factor: annotation of this region,
	either \code{"normal"} for no breakpoints or
	\code{"breakpoint"}
	for at least one breakpoint.}
    }}
  }
}
\source{
  Gudrun Schleiermacher and Isabelle Janoueix-Lerosey,
  Institut Curie.
}
\examples{
data(neuroblastoma)
pids <- levels(neuroblastoma$profiles$profile.id)[1:2]
pro <- subset(neuroblastoma$profiles,profile.id\%in\%pids)
ann <- subset(neuroblastoma$annotations,profile.id\%in\%pids)
library(ggplot2)
ggplot()+
  geom_tallrect(aes(xmin=min/1e6,xmax=max/1e6,
                    fill=annotation),data=ann)+
  geom_point(aes(position/1e6,logratio),pch=1,data=pro)+
  scale_fill_manual(values=c(breakpoint="red",
                             normal="yellow"))+
  scale_x_continuous(breaks=c(100,200))+
  facet_grid(profile.id~chromosome,
             scales="free_x",space="free")+
  theme_bw()+
  opts(axis.title.x=theme_text(vjust = 0),
       panel.margin=unit(0,"lines"))
}
\keyword{datasets}
