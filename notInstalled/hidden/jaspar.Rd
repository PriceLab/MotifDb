\name{jaspar} 
\docType{data}
\alias{jaspar}
\title{The JASPAR CORE transcription factor binding profiles}

\description{
459 transcription factor binding site matrices downdloaded the open-source JASPAR project.   See 
\url{http://jaspar.genereg.net/}.  These matrices are unevenly distributed across twenty organims, as 
shown below. In normal use, the MotifDb constructor loads all the profiles, including JASPAR_CORE.

JASPAR representes these profiles as \emph{position count matrices}.
We represent them as \emph{position frequency matrices} in which the
JASPAR-reported counts are normalized, after which, at every position
in the profile, the sum of frequencies of the four nucleotides (A, C,
G, T) is 1.0.  In the metadata (see example below) a 'count' is
stored, representing the number of sequences which contributed to the
profile.  This counts makes it possible to reproduce the original
position count matrix.

\tabular{rlll}{
      \tab \bold{ORGANISM}  \tab \bold{COUNT}\cr
      \tab Scerevisiae  \tab 177\cr
      \tab Dmelanogaster  \tab 125\cr
      \tab Hsapiens  \tab 66\cr
      \tab Mmusculus  \tab 47\cr
      \tab Rnorvegicus  \tab 8\cr
      \tab Zmays  \tab 6\cr
      \tab Athaliana  \tab 5\cr
      \tab Celegans  \tab 5\cr
      \tab Amajus  \tab 3\cr
      \tab Psativum  \tab 3\cr
      \tab Gallus  \tab 2\cr
      \tab Hroretzi  \tab 1\cr
      \tab Hvulgare  \tab 1\cr
      \tab Nsylvestris  \tab 1\cr
      \tab Ocuniculus  \tab 1\cr
      \tab Phybrida  \tab 1\cr
      \tab Rrattus  \tab 1\cr
      \tab Taestivam  \tab 1\cr
      \tab Vertebrata  \tab 1\cr
      \tab Xlaevis  \tab 1\cr
    }
}
	    

\examples{
  hits <- which (values (MotifDb)$dataSource=='JASPAR_CORE')
  matrices.jaspar <- MotifDb [hits]
    # examine one matrix
  print (matrices.jaspar [[100]])
  info <- values (matrices.jaspar) [100,]
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
  ScerTF,
  uniprobe
}

\keyword{datasets}
