\name{FlyFactorSurvey} 
\docType{data}
\alias{FlyFactorSurvey}
\title{The FlyFactorSurvey Collection}

\description{
  FlyFactorSurvey is a collection of 614 position count matrices
  from the University of Massachusetts Medical School, Department of
  Molecular Medicine.  The matrices were  obtained primarily via the
  Bacterial one-hybrid system.

The FlyFactorSurvey provides position count matrices, which MotifDb
converts into position frequency matrices, accompanied by an inferred
sequenceCount integer, one for each matrix, which can be used to
reconstitute the original count matrix if you wish.

The matrices have been determined by three different techniques:
\itemize{
  \item bacterial 1-hybrid, SANGER sequencing: 387 matrices
  \item bacterial 1-hybrid SOLEXA sequencing: 156 matrices
  \item DNASE I footprinting: 65 matrices
  }
}

\examples{
  hits <- which (values (MotifDb)$dataSource=='FlyFactorSurvey')
  matrices.ffs <- MotifDb [hits]
    # examine one matrix
  print (matrices.ffs [[100]])
  info <- values (matrices.ffs) [100,]
    # examine the metadata, skipping missing columns
  print (info [, which (!is.na (info))])

}

\seealso{
  MotifDb,
  query,
  subset,
  export,
  hPDI,
  jaspar,
  ScerTF,
  uniprobe
}
  

\references{
  Zhu, L.J. et al. Nucleic Acids Research, Jan 2011.  FlyFactorSurvey: a database of Drosophila transcription factor
  binding specificities determined using the bacterial one-hybrid
  system. (\url{http://www.ncbi.nlm.nih.gov/pubmed/21097781})

 \url{http://en.wikipedia.org/wiki/Bacterial_one-hybrid_system}
}



\keyword{datasets}
