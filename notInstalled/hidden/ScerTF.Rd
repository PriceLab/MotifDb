\name{ScerTF} 
\docType{data}
\alias{ScerTF}
\title{The ScerTF Collection}

\description{
  ScerTF, from the Stormo lab (\url{http://stormo.wustl.edu/}) is curated
  from multiple sources of transcription factor binding profiles. Quoting from
  \url{http://stormo.wustl.edu/ScerTF/}
  
  \cite{ScerTF catalogs over 1,200 position weight matrices (PWMs) for 196
different yeast transcription factors. We've curated 11 literature
sources, benchmarked the published position-specific scoring matrices
against in-vivo TF occupancy data and TF deletion experiments, and
combined the most accurate models to produce a single collection of the
best performing weight matrices for Saccharomyces cerevisiae.}

The ScerTF project provides these matrices as normalized position frequency matrices,
which we use without further transformation, and in which each column
(each position) sums to one.  No sequenceCount information is available

}
	    

\examples{
  hits <- which (values (MotifDb)$dataSource=='ScerTF')
  matrices.scertf <- MotifDb [hits]
    # examine one matrix
  print (matrices.scertf [[100]])
  info <- values (matrices.scertf) [100,]
    # examine the metadata, skipping missing columns
  print (info [, which (!is.na (info))])
  
}

\seealso{
  MotifDb,
  query,
  subset,
  export,
  flyFactorSurvey,
  hPDI,
  jaspar,
  uniprobe
}

\keyword{datasets}
