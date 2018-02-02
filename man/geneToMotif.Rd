\name{geneToMotif}
\alias{geneToMotif,MotifList-method}
\alias{geneToMotif}
\title{geneToMotif}
\description{
Using either of our two sources ("MotifDb" or "TFClass") retrieve the
names of the transcription factor binding motifs associated with the
gene symbol for each transcription factor.  Slightly different
information is returned in each case but the columns "geneSymbol",
"motif", "pubmedID", "source" are returned by both sources.
The TFClass source is described here:
\url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4383905/}.
The MotifDb source is in fact the usually 1:1 gene/motif mapping
provided by each of the data sources upon which MotifDb is built.
}
\usage{
\S4method{geneToMotif}{MotifList}(object, geneSymbols, source, ignore.case)
}
\arguments{
  \item{object}{a \code{MotifList} object.}
  \item{geneSymbols}{a \code{character} string}
  \item{source}{a \code{character} string, either 'MotifDb' or 'TFclass' (case insensitive)}
  \item{ignore.case}{a \code{logical} variable, default FALSE, guiding gene name matching}
}

\value{
  A data.frame with these columns: geneSymbol, motif, pubmedID, source.
  The MotifDb source alos include dataSource and organism.
}
\author{Paul Shannon}

\examples{
  genes <- c("ATF5", "FOS")
  geneToMotif(MotifDb, genes, source="TFClass")
  geneToMotif(MotifDb, genes, source="MotifDb")
}
\seealso{
  MotifDb,
  motifToGene,
  associateTranscriptionFactors,
  subset,
  query
}

\keyword{utilities}
