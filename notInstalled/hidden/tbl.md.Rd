\name{tbl.md} 
\docType{data}
\alias{tbl.md}
\title{The named list which holds the position frequency tbl.md
  fundamental to this class}

\description{
  A DataFrame with these columns
  \itemize{
    \item name:  3 parts: source-organism-gene, with a source-specific ID added to the gene, occasionally augmented by a uniquifying 4-letter code.
    \item nativeName: the name by which the matrix is described in its source database
    \item geneSymbol: always provided by the source, this is the sometimes ambiguous gene sybmol of the transcription factor, or binding molecule
    \item sequenceCount: when known, this describes the number of aligned sequences which contributed to the motif
    \item organism:  of the binding molecule
    \item bindingMolecule: a standard identifier
    \item bindingMoleculeIdType: which standard?  uniprot, entrez geneID, SGD orf, flybase FBGn, ....
    \item bindingDomainType:  when known, of the transcription factor's DNA-binding domain.
    \item dataSource:  from which the collection of matrices and accompanying metadata were obtained
    \item experimentType: sometimes known:  bacterial one-hybrid, protein-binding microarray, ChIP-chip, ...
    \item pubmedID: of the journal article describing the experiment from which the profile comes
    \item tfFamily: sometimes known, the class of this transcription factor
    }

}
	    

\examples{
  library (MotifDb)
  print (values (MotifDb) [1,])
}

\keyword{datasets}
