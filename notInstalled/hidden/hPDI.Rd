\name{hPDI} 
\docType{data}
\alias{hPDI}
\title{hPDI: Human Protein-DNA Interactome (hPDI)}

\description{
As described on the hPDI website \url{http://bioinfo.wilmer.jhu.edu/PDI}: The hPDI database holds experimental protein-DNA interaction data for humans identified by protein microarray assays. The current release of hPDI contains 17,718 protein-DNA interactions for 1013 human DNA-binding proteins. These DNA-binding proteins include 493 human transcription factors (TFs) and 520 unconventional DNA binding proteins (uDBPs). Here TFs are defined as proteins either are annotated as "regulation of transcription" defined by GO database or contain known DNA binding domains defined by Pfam database and, uDBPs are defined proteins that have not been previously found to binding DNA specifically.
\url{http://jaspar.genereg.net/}. 
}
	    

\examples{
  hits <- which (values (MotifDb)$dataSource=='hPDI')
  matrices.hpdi <- MotifDb [hits]
    # examine one matrix
  print (matrices.hpdi [[100]])
  info <- values (matrices.hpdi) [100,]
    # examine the metadata, skipping missing columns
  print (info [, which (!is.na (info))])


}

\seealso{
  MotifDb,
  query,
  subset,
  export,
  flyFactorSurvey,
  jaspar,
  ScerTF,
  uniprobe
}

\keyword{datasets}
