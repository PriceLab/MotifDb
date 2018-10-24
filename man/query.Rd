\name{query}
\alias{query,MotifList-method}
\alias{query}
\title{query}
\description{
A very general search tool, returning all matrices whose metadata, in ANY
column, is matched by the query string.
}
\usage{
\S4method{query}{MotifList}(object, andStrings, orStrings, notStrings, ignore.case=TRUE)
}
\arguments{
  \item{object}{a \code{MotifList} object.}
  \item{andStrings}{a \code{character} string vector, length one or
    more, every element of which must be found in the metadata}
  \item{orStrings}{a \code{character} string vector, length one or
    more, any element of which must be found in the metadata}
  \item{notStrings}{a \code{character} string vector, length one or
    more, none of which may appear in the metadata}
  \item{ignore.case}{a \code{logical} value, default TRUE}

}

\value{
A list of the matrices
}
\author{Paul Shannon}

\examples{
  matrices.human <- query(MotifDb, 'hsapiens')
  matrices.sox4 <- query(MotifDb, 'sox4')
  uniprobe.sox.matrices <- query(MotifDb, c('uniprobe', 'sox'))
    # two approaches to selectinve extraction of TFEB matrices
  tfeb.human.1 <- query(MotifDb, andStrings=c("TFEB", "hsapiens"), notStrings=c("hpdi", "jolma", "cisbp"))
  tfeb.human.2 <- query(MotifDb, andStrings=c("TFEB", "hsapiens"), orStrings=c("hocomoco", "jaspar", "swissregulon"),
                        notStrings="2016")

}
\seealso{
  MotifDb,
  subset,
  export,
  flyFactorSurvey,
  hPDI,
  jaspar,
  ScerTF,
  uniprobe
}

\keyword{utilities}
