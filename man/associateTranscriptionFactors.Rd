\name{associateTranscriptionFactors}
\alias{associateTranscriptionFactors,MotifList-method}
\alias{associateTranscriptionFactors}
\title{associateTranscriptionFactors}
\description{
In the analysis of, or exploration of gene regulatory networks, one
often creates a data.frame of possible genomic regulatory sites, genomic
locations where a TF binding motif matches some DNA sequence.  A common
next step is to associate each of these motifs with its related
transcription factor/s.     We provide two sources for those
relationships.  When you specify the "MotifDb" source, we return
the motif/TF relationships provided by each of the constituent
public MotifDb sources.  When you specify the "TFClass" source,
transcription factor family memberships (described in
\url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4383905/}) are -
sometimes expansively - provided for each motif you supply.

This method uses, and therefore expects, different columns of the incoming data.frame
to be used with each method.  The MotifDb source uses the "motifName" column of
the incoming data.frame.  The TFClass source expects a "shortName"
column in the incoming database.

A new column, "geneSymbol", is added to the incoming data.frame.  This
new column identifies the transcription factor associated with the motif
for each row in the data.frame.

}
\usage{
\S4method{associateTranscriptionFactors}{MotifList}(object, tbl.withMotifs, source, expand.rows, motifColumnName="motifName")
}
\arguments{
  \item{object}{a \code{MotifList} object.}
  \item{tbl.withMotifs}{a \code{data.frame}}
  \item{source}{a \code{character} string, either "MotifDb" or "TFClass" (case insensitive)}
  \item{expand.rows}{a \code{logical} value, recommended especially for the TFClass source,
    in which sometimes many TFs are mapped to the same motif}
  \item{motifColumnName}{a \code{character} string identifying the
    column in tbl.withMotifs which contains the motifs to be associated
    with transcription factors}
}

\value{
  A data.frame with one column ("geneSymbol") and possibly multiple rows added
}
\author{Paul Shannon}

\examples{
  tbl.tfClassExample <- data.frame(motifName=c("MA0006.1", "MA0042.2", "MA0043.2"),
                                   chrom=c("chr1", "chr1", "chr1"),
                                   start=c(1000005, 1000085, 1000105),
                                   start=c(1000013, 1000092, 1000123),
                                   score=c(0.85, 0.92, 0.98),
                                   stringsAsFactors=FALSE)
      # here we illustrate how to add a column with the required name:
  tbl.tfClassExample$shortMotif <- tbl.tfClassExample$motifName
  tbl.out <- associateTranscriptionFactors(MotifDb, tbl.tfClassExample, source="TFClass",
                                           expand.rows=TRUE)
  dim(tbl.out)   #    MANY tfs mapped, mostly FOX family genes
  tbl.motifDbExample <- data.frame(motifName=c("Mmusculus-jaspar2016-Ahr::Arnt-MA0006.1",
                                               "Hsapiens-jaspar2016-FOXI1-MA0042.2",
                                               "Hsapiens-jaspar2016-HLF-MA0043.2"),
                                    chrom=c("chr1", "chr1", "chr1"),
                                    start=c(1000005, 1000085, 1000105),
                                    start=c(1000013, 1000092, 1000123),
                                    score=c(0.85, 0.92, 0.98),
                                    stringsAsFactors=FALSE)

   tbl.out <- associateTranscriptionFactors(MotifDb, tbl.motifDbExample, source="MotifDb",
                                            expand.rows=TRUE)
   dim(tbl.out)   # one new column ("geneSymbol"), no new rows
}
\seealso{
  MotifDb,
  geneToMotif,
  motifToGene,
  subset,
  query
}

\keyword{utilities}
