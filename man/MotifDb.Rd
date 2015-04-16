\name{MotifDb}
\alias{MotifDb}
\alias{MotifDb-package}
\title{MotifDb: An Annotated Collection of DNA-binding sequence motifs}
\description{
  Approximately 2000 position frequency matrices collected from public
  sources, with ample accompanying metadata, and search and export
  capabilities provided.

}
\details{
 MotifDb is an R object of class MotifList, whose entries are numeric
 matrices, accompanied by a 'parallel' metadata structure, a DataFrame, in
 which each row provides information about the corresponding matrix. 
This object is automatically created and fully populated by data from
five public sources (see below) when the package is loaded into your R environment via the
\code{library} call. 
The matrices are obtained from six public sources:
\tabular{ll}{
FlyFactorSurvey: \tab 614\cr
hPDI: \tab 437\cr
JASPAR_CORE: \tab 459\cr
jolma2013: \tab 843\cr
ScerTF: \tab 196\cr
stamlab: \tab 683\cr
UniPROBE: \tab 380\cr
}

Representing primarily four organsisms:
\tabular{ll}{
Dmelanogaster: \tab 739\cr
Hsapiens: \tab 1898\cr
Scerevisiae: \tab 464\cr
Mmusculus: \tab 462\cr
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


The names of the matrices are the same as rownames of the metadata
DataFrame, and have been chosen to balance the needs of concision and
full description, including the organism in which the motif was
discovered, the data source, and the name of the motif in the data
source from which it was obtained.  For example:
"Hsapiens-JASPAR_CORE-SP1-MA0079.2" and "Scerevisiae-ScerTF-GSM1-badis".

Subsets of the Matrices may be obtainted in several ways:

\itemize{
  \item By integer index, eg, \code{MotifDb [[1]]}
  \item By query, eg, \code{as.list (query (MotifDb, 'FBgn0000014'))}
  \item (Interactively only) by subset \code{as.list (subset (MotifDb,
    geneSymbol=='Abda' & !is.na (pubmedID)))}
}

The matrices are stored in a \code{SimpleList} which has semantics very
similar to the familiar list of R base.  To examine a matrix, however,
you must sidestep the MotifDb \code{show} method.  These three commands
display quite different results:

\preformatted{
> MotifDb [1]
MotifDb object of length 1
| Created from downloaded public sources: 2012-Jul6
| 1 position frequency matrices from 1 source:
|    FlyFactorSurvey:    1
| 1 organism/s
|      Dmelanogaster:    1
Dmelanogaster-FlyFactorSurvey-ab_SANGER_10_FBgn0259750 

> MotifDb [[1]]
    1    2    3    4 5 6 7 8 9   10   11   12   13   14   15   16   17   18   19   20   21
A 0.0 0.50 0.20 0.35 0 0 1 0 0 0.55 0.35 0.05 0.20 0.45 0.20 0.10 0.40 0.40 0.25 0.50 0.30
C 0.3 0.15 0.25 0.00 1 1 0 0 0 0.10 0.65 0.70 0.45 0.25 0.10 0.25 0.25 0.10 0.10 0.25 0.25
G 0.4 0.05 0.50 0.65 0 0 0 1 1 0.00 0.00 0.05 0.05 0.15 0.05 0.20 0.05 0.15 0.55 0.15 0.45
T 0.3 0.30 0.05 0.00 0 0 0 0 0 0.35 0.00 0.20 0.30 0.15 0.65 0.45 0.30 0.35 0.10 0.10 0.00

> as.list (MotifDb [1])
$`Dmelanogaster-FlyFactorSurvey-ab_SANGER_10_FBgn0259750`
    1    2    3    4 5 6 7 8 9   10   11   12   13   14   15   16   17   18   19   20   21
A 0.0 0.50 0.20 0.35 0 0 1 0 0 0.55 0.35 0.05 0.20 0.45 0.20 0.10 0.40 0.40 0.25 0.50 0.30
C 0.3 0.15 0.25 0.00 1 1 0 0 0 0.10 0.65 0.70 0.45 0.25 0.10 0.25 0.25 0.10 0.10 0.25 0.25
G 0.4 0.05 0.50 0.65 0 0 0 1 1 0.00 0.00 0.05 0.05 0.15 0.05 0.20 0.05 0.15 0.55 0.15 0.45
T 0.3 0.30 0.05 0.00 0 0 0 0 0 0.35 0.00 0.20 0.30 0.15 0.65 0.45 0.30 0.35 0.10 0.10 0.00

}

There are fifteen kinds of metadata -- though not all matrices have a
full complement: not all of the public sources are complete in this
regard. The information falls into these categories, using the
\emph{Dmelanogaster-FlyFactorSurvey-ab_SANGER_10_FBgn0259750} entry as
an example (see below for the associated position frequency matrix):

\enumerate{
   \item providerName: "ab_SANGER_10_FBgn0259750"
   \item providerId: "FBgn0259750"
   \item dataSource: "FlyFactorSurvey"
   \item geneSymbol: "Ab"
   \item geneId: "FBgn0259750"
   \item geneIdType: "FLYBASE"
   \item proteinId: "E1JHF4"
   \item proteinIdType: "UNIPROT"
   \item organism: "Dmelanogaster"
   \item sequenceCount: 20
   \item bindingSequence: NA
   \item bindingDomain: NA
   \item tfFamily: NA
   \item experimentType: "bacterial 1-hybrid, SANGER sequencing"
   \item pubmedID: NA
}

} %% details

\seealso{
  query,
  subset,
  export,
  flyFactorSurvey,
  hPDI,
  jaspar,
  ScerTF,
  uniprobe
}

\examples{
     # are there any matrices for Sox4?  we find two
   mdb.sox4 <- MotifDb [grep ('sox4', values (MotifDb)$geneSymbol, ignore.case=TRUE)]
     # the same two matrices can be obtained this way also
   if (interactive ()) 
     mdb.sox4 <- subset (MotifDb, tolower(geneSymbol)=='sox4')
     # and like this
   mdb.sox4 <- query (MotifDb, 'sox4')  # matches against all fields in the metadata
     # implicitly invoke the 'show' method
   mdb.sox4
     # get their full names
   names (mdb.sox4)
     # examine their metadata
   values (mdb.sox4)
     # examine the matrices with names include
   as.list (mdb.sox4)
     # export the matrices in meme format 
   destination.file = tempfile ()
   export (mdb.sox4, destination.file, 'meme')

}


\references{

\itemize{

\item Neph S, Stergachis AB, Reynolds A, Sandstrom R, Borenstein E,
Stamatoyannopoulos JA. Circuitry and dynamics of human transcription factor regulatory networks.
Cell. 2012 Sep 14;150(6):1274-86.
  
\item Portales-Casamar E, Thongjuea S, Kwon AT, Arenillas D, Zhao X, Valen E, Yusuf D, Lenhard B, Wasserman WW, Sandelin A. JASPAR 2010: the greatly expanded open-access database of transcription factor binding profiles. Nucleic Acids Res. 2010 Jan;38(Database issue):D105-10. Epub 2009 Nov 11.

\item Robasky K, Bulyk ML. UniPROBE, update 2011: expanded content and search tools in the online database of protein-binding microarray data on protein-DNA interactions. Nucleic Acids Res. 2011 Jan;39(Database issue):D124-8. Epub 2010 Oct 30.

\item Spivak AT, Stormo GD. ScerTF: a comprehensive database of benchmarked position weight matrices for Saccharomyces species. Nucleic Acids Res. 2012 Jan;40(Database issue):D162-8. Epub 2011 Dec 2.

\item Xie Z, Hu S, Blackshaw S, Zhu H, Qian J. hPDI: a database of experimental human protein-DNA interactions. Bioinformatics. 2010 Jan 15;26(2):287-9. Epub 2009 Nov 9.

\item Zhu LJ, et al. 2011. FlyFactorSurvey: a database of Drosophila transcription factor binding specificities determined using the bacterial one-hybrid system. Nucleic Acids Res. 2011 Jan;39(Database issue):D111-7. Epub 2010 Nov 19.

\item Jolma A, et al. 2013. DNA-binding specificities of human transcription factors. Cell  2013 Jan 17.

} %% itemize

} %% references


\keyword{datasets}

