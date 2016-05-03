\name{subset}
\alias{subset,MotifList-method}
\title{subset}
\description{
  An analog of the base package subset method, this version will return
  all the matrices whose metadata match the (possibly intricate) logical
  expression in the "subset" argument.

  Note: just as with the base subset method, this method is unreliable
  except when used interactively.  Batch, script or other programmatic
  use of this function is to be avoided.
}
\usage{
\S4method{subset}{MotifList}(x, subset, select, drop=FALSE, ...)
}
\arguments{
  \item{x}{a \code{MotifList} object. }
  \item{subset}{a \code{logical} expression whose terms are
    predicates on the column names of the metadata table}
  \item{select, drop, ...}{these are ignored, appearing here only in
    fidelity to the generic definition of the method.}
  
}

\value{
A list of the matrices whose metadata satisfies the supplied subset
}
\author{Paul Shannon}

\examples{
  mdb <- MotifDb
  if (interactive ()) {
    matrices <- subset (mdb, dataSource=='UniPROBE')
    egr1.matrices <- subset (mdb, geneSymbol=='Egr1')
    jaspar.egr1.matrices <- subset (mdb, geneSymbol=='Egr1' &
                                         dataSource == 'JASPAR_CORE')
      # one of the mouse egr1 matrices has a geneSymbol 'Zif268', but
      # has the proper entrez geneId.  
    all.egr1.matrices <- subset (mdb, geneId=='13653')
    }
}


\seealso{
  MotifDb,
  query,
  export,
  flyFactorSurvey,
  hPDI,
  jaspar,
  ScerTF,
  uniprobe
}


\keyword{utilities}
