\name{matchMotif}
\alias{matchMotif,MotifList-method}
\alias{matchMotif}
\title{matchMotif}
\description{
A very general search tool, returning all matrices whose metadata, in ANY
column, is matched by the query string.
}
\usage{
\S4method{matchMotif}{MotifList}(object, motifs, genomeName, regions, pval.cutoff, fimoDataFrameStyle=FALSE)
}
\arguments{
  \item{object}{a \code{MotifList} object.}
  \item{motifs}{a \code{MotifList} typically the result of a MotifDb query}
  \item{genomeName}{a \code{character} string vector, e.g. "hg38" or "mm10"}
  \item{regions}{a \code{GRange} specifies regions in which to try to match the motifs}
  \item{pval.cutoff}{a \code{numeric} value, for example, 1e-5}
  \item{fimoDataFrameStyle}{a \code{logical} value, default FALSE}

}

\value{
A list of the matrices
}
\author{Paul Shannon}

\examples{
   motifs <- query(MotifDb, "ZNF263", c("jaspar2018", "swissregulon"))
   gr.region <- GRanges(seqnames="chr1", IRanges(start=47229520, end=47229560))
   gr.match <- matchMotif(MotifDb, motifs, "hg38", gr.region, 1e-5)
   tbl.match <- matchMotif(MotifDb, motifs, "hg38", gr.region, 1e-5, fimoDataFrameStyle=TRUE)
}
\seealso{
  MotifDb,
  query,
  subset,
  export,
  flyFactorSurvey,
  hPDI,
  jaspar,
  ScerTF,
  uniprobe
}

\keyword{utilities}
