\name{query}
\alias{query,MotifList-method}
\alias{query}
\title{query}
\description{
A very general search tool, returning all matrices whose metadata, in ANY
column, is matched by the query string.  
}
\usage{
\S4method{query}{MotifList}(object, queryString, ignore.case=TRUE)
}
\arguments{
  \item{object}{a \code{MotifList} object.}
  \item{queryString}{a \code{character} string}
  \item{ignore.case}{a \code{logical} value, default TRUE}
  
}

\value{
A list of the matrices
}
\author{Paul Shannon}

\examples{
  mdb <- MotifDb
  matrices.human = query (mdb, 'hsapiens')
  matrices.sox4 = query (mdb, 'sox4')
  uniprobe.sox.matrices = query (query (mdb, 'uniprobe'), '^sox')
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
