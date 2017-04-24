\name{export}
\alias{export}
\alias{export,MotifList,connection,character-method}
\alias{export,MotifList,character,character-method}
\alias{export,MotifList,missing,character-method}
\title{export}
\description{
Writes all matrices in the supplied list, in the specified format, to
the specified connection.
}
\usage{
\S4method{export}{MotifList,connection,character}(object, con, format, ...)
\S4method{export}{MotifList,character,character}(object, con, format, ...)
\S4method{export}{MotifList,missing,character}(object, con, format, ...)
}
\arguments{
  \item{object}{a \code{MotifList} object. }
  \item{con}{either a \code{file connection} or a \code{filename}
    or missing, implying stdout.}
  \item{format}{a \code{character} string, currently only 'meme' and
    'transfac', which both produce the same result}
  \item{...}{ignore this}
}

\value{
The matrices list is written to the specified connection in the
specified format.
}
\author{Paul Shannon}

\examples{
  library (MotifDb)
    # identify all the SOX genes
  sox.indices = grep ('^sox', values (MotifDb)$geneSymbol, ignore.case=TRUE)
  matrices = MotifDb [sox.indices]
  export (matrices, con='SoxGenes-meme.txt', format='meme')
}

\seealso{
  MotifDb,
  query,
  subset,
  flyFactorSurvey,
  hPDI,
  jaspar,
  ScerTF,
  uniprobe
}

\keyword{utilities}
