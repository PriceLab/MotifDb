\name{MotifList-class}
\docType{class}

\alias{class:MotifList}
\alias{MotifList-class}
\alias{show,MotifList-method}

%%\alias{export,MotifList-method}
%%\alias{show,MotifList-method}
%%\alias{subset,MotifList-method}
%%\alias{query,MotifList-method}

\title{MotifList}

\description{
  A direct subclass of SimpleList, having no extra slots, in which
  listData is a list of position frequency matrices (PFMs), and the
  elementMetadata slot is a DataFrame with fifteen columns describing
  each matrix.  Upon loading the MotifDb class, one MotifList object 
  is instantiated and filled with matrices and their metadata.  
  There should be no need for users to explicitly create objects of this
  class.  When you load the MotifDb package, a fully-populated instance
  of this class is created, with > 2000 matrices with metadata
}


\section{Methods}{
  \describe{
    \item{}{
      \code{subset(x)}: extract matrices by metadata.
     }
   \item{}{
      \code{export(x)}: write matrices
    }
   \item{}{
      \code{show(x)}: describe matrices compactly
    }
   \item{}{
      \code{query(x)}: find matrices
    }
   }
 }

\author{Paul Shannon}

\examples{
    # Examine the number of matrices contributed by each source.
  print (table (values (MotifDb)$dataSource))

}
\keyword{methods}
\keyword{classes}
