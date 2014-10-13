\name{uniprobe} 
\docType{data}
\alias{uniprobe}
\title{The UniPROBE Collection}

\description{
  UniPROBE is a collection of 380 position frequency matrices derived
  from universal PBM (protein binding microarray) experiments.  UniPROBE
  is hosted by the Bulyk lab, which is part of the
  Division of Genetics at Brigham and Women's Hospital, Boston. 
  
UniPROBE provides position frequency matrices, which are used here
without modification.  Sequence count data, and the actual sequences of
the DNA plated on to the microarray, are available, sometimes
accompanined by the protein binding ('insertion') sequence, but these
are not yet part of MotifDb.

The 380 matrices are distributed like this:

\tabular{rlll}{

    \tab \bold{ORGANISM} \tab \bold{COUNT} \cr
    \tab Mmusculus \tab 281 \cr
    \tab Scerevisiae \tab 89 \cr
    \tab Scerevisiae;Hsapiens;Mmusculus \tab 5 \cr
    \tab Apicomplexa \tab 3 \cr
    \tab Celegans \tab 1 \cr
    \tab Hsapiens \tab 1 \cr
  }

  Two UniPROBE matrices appear to have rounding errors.  Both are from
  the PNAS08 paper.  They are \emph{almost} normalized, but not quite.
  We present them as we found them.
  
\preformatted{
Protein: Cgd2_3490	Seed k-mer: CGCATGCA	Enrichment Score: 0.499452
A:	0.31	0.27	0.28	0	0.01	0	0.99	0	0	0	0.65	0.16	0.11	0.13
C:	0.31	0.29	0.43	0.64	0	0.99	0	0	0	0.98	0	0.27	0.37	0.48
G:	0.14	0.2	0.11	0.01	0.99	0	0	0	0.99	0	0.34	0.18	0.31	0.17
T:	0.24	0.24	0.19	0.35	0	0	0	0.99	0	0.01	0	0.38	0.21	0.22

Protein: PF14_0633	Seed k-mer: TGCATGCA	Enrichment Score: 0.499939
A:	0.25	0.29	0.38	0	0.09	0.01	0.99	0.01	0	0	0.74	0.32	0.14	0.29
C:	0.21	0.24	0.2	0.26	0	0.98	0	0	0	0.91	0	0.13	0.17	0.26
G:	0.26	0.13	0.12	0	0.91	0	0	0	0.98	0	0.26	0.34	0.44	0.17
T:	0.28	0.34	0.29	0.74	0	0	0.01	0.99	0.01	0.09	0	0.2	0.25	0.28

}
  
}


\examples{
  hits <- which (values (MotifDb)$dataSource=='UniPROBE')
  matrices.uniprobe <- MotifDb [hits]
    # examine one matrix
  print (matrices.uniprobe [[100]])
  info <- values (matrices.uniprobe) [100,]
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
  ScerTF
}
  

\references{
  Newburger DE and Bulyk ML. UniPROBE:  an online database of protein binding microarray data on protein-DNA
  interactions. Nucl. Acids Res. 2009 37: D77-D82;
  doi:10.1093/nar/gkn660. Epub October 8, 2008).

  \url{http://en.wikipedia.org/wiki/Protein-binding_microarray}
}

\keyword{datasets}
