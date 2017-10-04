\name{motifToGene}
\alias{motifToGene,MotifList-method}
\alias{motifToGene}
\title{motifToGene}
\description{
Using either of our two sources ("MotifDb" or "TFClass") this method retrieves the
the transcription factor (its gene symbol) for each of the supplied motifs.
 Slightly different
information is returned in each case but the columns "geneSymbol",
"motif", "pubmedID", "source" are returned by both.
The TFClass source is described here:
\url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4383905/}.
The MotifDb source is in fact the (typically) 1:1 gene/motif mapping
provided by each of the data sources upon which MotifDb is built.
}
\usage{
\S4method{motifToGene}{MotifList}(object, motifs, source)
}
\arguments{
  \item{object}{a \code{MotifList} object.}
  \item{motifs}{a \code{character} string}
  \item{source}{a \code{character} string, either 'MotifDb' or "TFclass' (case insensitive)}
}

\value{
  A data.frame with these columns: geneSymbol, motif, pubmedID, source.
  The MotifDb source also include dataSource and organism.
}
\author{Paul Shannon}

\examples{
  motifs <- c("MA0592.2", "ELF1.SwissRegulon", "UP00022")
  motifToGene(MotifDb, motifs, source="TFClass")
  motifToGene(MotifDb, motifs, source="MotifDb")
}
\seealso{
  MotifDb,
  geneToMotif,
  associateTranscriptionFactors,
  subset,
  query
}

\keyword{utilities}
