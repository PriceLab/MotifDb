\name{MotifDb-package}
\alias{MotifDb-package}
\docType{package}
\title{
  MotifDb:  Protein-DNA binding Motifs Collected from Public Sources
}
\description{
Sequence motifs summarize protein-DNA binding sites using 4-row by
n-column position weight matrices (PWMs) with one row for
each of the four possible nucelotides,
and one column for each position in the sequence.  Matrix values are
traditionally either raw counts, frequencies which sum to 1.0 for each
column, or weights -- often log-likelihood or log-odds scores.  Here we provide
more than two thousand such matrices, each accompanied by metadata,
from five sources and twenty-two organisms (but
primiarly fly, human, yeast and mouse).
}
\details{
The matrices are obtained from ive public sources:
\tabular{ll}{
FlyFactorSurvey: \tab 614\cr
hPDI: \tab 437\cr
JASPAR_CORE: \tab 459\cr
ScerTF: \tab 196\cr
UniPROBE: \tab 380\cr
}

Representing primarily four organsisms:
\tabular{ll}{
Dmelanogaster: \tab 739\cr
Hsapiens: \tab 505\cr
Scerevisiae: \tab 464\cr
Mmusculus: \tab 329\cr
Rnorvegicus: \tab 8\cr
Celegans: \tab 7\cr
Zmays: \tab 6\cr
Athaliana: \tab 5\cr
Psativum: \tab 3\cr
Amajus: \tab 3\cr
Pfalciparum: \tab 2\cr
Gallus: \tab 2\cr
Xlaevis: \tab 1\cr
Vertebrata: \tab 1\cr
Taestivam: \tab 1\cr
Rrattus: \tab 1\cr
Phybrida: \tab 1\cr
Ocuniculus: \tab 1\cr
Nsylvestris: \tab 1\cr
Hvulgare: \tab 1\cr
Hroretzi: \tab 1\cr
Cparvum: \tab 1\cr
}

All the matrices are stored as position frequency matrices, in which
each columm (each position) sums to 1.0.  When the number of sequences
which contributed to the motif are known, that number will be found in
the matrix's metadata.  With this information, one can transform the
matrices into either PCM (position count matrices), or PWM (position
weight matrices), also known as PSSM (position-specific-scoring
matrices). The latter transformation requires that a model of the
background distribution be known, or assumed.


Index:
\preformatted{
MotifDb                                   MotifDb
MotifList                                 MotifList
MotifList-class                           MotifList objects
                                          export
FlyFactorSurvey                           The FlyFactorSurvey Collection
hPDI                                      hPDI: Human Protein-DNA Interactome (hPDI)
jaspar                                    The JASPAR CORE transcription factor binding profiles
ScerTF                                    The ScerTF Collection
uniprobe                                  The UniPROBE Collection
export,MotifList,connection,character-method
length,MotifList-method                   length
subset,MotifList-method                   subset
matrices                                  The named list which holds the position frequency matrices fundamental to this
                                          class
tbl.md                                    The named list which holds the position frequency tbl.md fundamental to this class
}

Further information is available in the following vignettes:
\tabular{ll}{
\code{MotifDb} \tab MotifDb (source)\cr
}
}
\author{
Paul Shannon
Maintainer: <pshannon@fhcrc.org>
}

\examples{
  # load all matrices from all sources, with metadata
  # get a summary 
MotifDb
  # learn the name of the first matrix
names (MotifDb [1])
  # examine the matrix
MotifDb [[1]]
  # see the metadata
values (MotifDb) [1,]
  # transform the frequency matrix into a count matrix
MotifDb [[1]] * values (MotifDb) [1,]$sequenceCount
}
