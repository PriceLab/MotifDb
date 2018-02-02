setGeneric('query', signature='object', function(object, queryString, ignore.case=TRUE)
              standardGeneric ('query'))
setGeneric('motifToGene', signature='object', function(object, motifs, source) standardGeneric('motifToGene'))
setGeneric('geneToMotif', signature='object', function(object, geneSymbols, source, ignore.case=FALSE) standardGeneric('geneToMotif'))
setGeneric('associateTranscriptionFactors', signature='object',
           function(object, tbl.withMotifs,  source, expand.rows) standardGeneric('associateTranscriptionFactors'))
#------------------------------------------------------------------------------------------------------------------------
setClass ('MotifList',
          contains='SimpleList',
          representation (elementMetadata='DataFrame',
                          manuallyCuratedGeneMotifAssociationTable="data.frame"),
          prototype(elementType='matrix')
          )

setValidity('MotifList', function(object) {
    msg = NULL
      ## what makes for a valid MotifList?
    if (length (object) ==  0)
       return (TRUE)
    if (is.null(names(object)))
        msg = c(msg, '"names()" must not be NULL')
    else if (any(duplicated(names(object))))
        msg = c(msg, 'all "names()" must be unique')
    if (!all(sapply(object, is.matrix)))
        msg = c(msg, 'all matrices must be of class "matrix"')
    if (!all(sapply(object, nrow) == 4))
        msg = c(msg, 'all matrices must have 4 rows')
      # all columns of all matrices should be normalized, summing to one.
      # in fact, 2/2086 matrices
      #   Cparvum-UniPROBE-Cgd2_3490.UP00395
      #   Pfalciparum-UniPROBE-PF14_0633.UP00394
      # fail to pass this test.  rounding to one decimal place allows these matrices,
      # also, to pass.   round (0.98, digits=1) --> 1.0
      # see the UniPROBE manpage for the full story on these two matrices
    ok = sapply(object, function(elt) all (round (colSums(elt), digits=1) == 1))
    if (!all(ok))
       msg = c(msg, 'all elements must have colSums equal to 1')

    if (is.null(msg)) TRUE else msg
})
#-------------------------------------------------------------------------------
MotifList = function (matrices=list(), tbl.metadata=data.frame ())
{
  if (nrow (tbl.metadata) == 0)
    df = DataFrame ()
  else
    df = DataFrame (tbl.metadata, row.names = rownames (tbl.metadata))

  tbl.tfClass.filename <- system.file(package="MotifDb", "extdata", "tfClass.tsv")
  stopifnot(file.exists(tbl.tfClass.filename))
  tbl.tfClass <- read.table(tbl.tfClass.filename, header=TRUE, as.is=TRUE, sep="\t")

  object = new ('MotifList', listData=matrices, elementMetadata = df,
                manuallyCuratedGeneMotifAssociationTable=tbl.tfClass)

  if (length (matrices) == 0)
    names (object) = character ()
  else
    names (object) = rownames (tbl.metadata)

  object

} # ctor
#-------------------------------------------------------------------------------
setMethod ('subset', signature = 'MotifList',

  function (x, subset, select, drop = FALSE, ...) {

    if (missing (subset))
      i = TRUE
    else {
      i = eval(substitute (subset), elementMetadata (x), parent.frame ())
      i = try(as.logical(i), silent=TRUE)
      if (inherits(i, 'try-error'))
        stop('"subset" must be coercible to logical')
      i = i & !is.na(i)
      } # else
    return (x [i])
  })

#-------------------------------------------------------------------------------
# transpose 4-row matrices to 4-column, so that there are as many rows as
# there are nucleotides in the motif sequence. meme requires normalized matrices
# exactly as we store them
transformMatrixToMemeRepresentation = function (m)
{
  return (t (m))

} # transformMatrixToMemeRepresentation
#-------------------------------------------------------------------------------
# http://stuff.mit.edu/afs/athena/software/meme_v3.5.4/etc/meme-explanation.html
# The motif itself is a position-specific probability matrix giving, for each
# position in the pattern, the observed  frequency ('probability') of each
# possible letter. The probability matrix is printed 'sideways'--columns
# correspond  to the letters in the alphabet (in the same order as shown in
# the simplified motif) and rows corresponding to the  positions of the motif,
# position one first. The motif is preceded by a line starting with
# 'letter-probability matrix:' and containing the length of the alphabet,
# width of the motif, number of occurrences of the motif, and the E-value
# of the motif.
matrixToMemeText = function (matrices)
{
  matrix.count = length (matrices)

    # incoming matrices have nucleotide rows, position columns.  meme
    # format, however, requires position-rows, and nucleotide-columns
    # calculate the number of lines of text by counting columns
  total.transposed.matrix.rows = sum (as.integer (sapply (matrices, ncol)))
  predicted.line.count = 12 + (3 * length (matrices)) +
                           total.transposed.matrix.rows
  #s = vector ('character', predicted.line.count)
  s = character (predicted.line.count)

  s [1] = 'MEME version 4'
  s [2] = ''
  s [3] = 'ALPHABET= ACGT'
  s [4] = ''
  s [5] = 'strands: + -'
  s [6] = ''
  s [7] = 'Background letter frequencies'
  s [8] = 'A 0.250 C 0.250 G 0.250 T 0.250 '
  s [9] = ''

  index = 10
  for (name in names (matrices)) {
       # transpose the frequency matrix version of the incoming matrix,
       # hence 'tfMat'
    tfMat = transformMatrixToMemeRepresentation (matrices [[name]])
       # meme output may be used by tomtom, which uses matrix names as
       # part of image filenames. removed all file-system-unfriendly
       # characters here
    fixed.name = gsub ('\\/', '_', name)
    s [index] = sprintf ('MOTIF %s', fixed.name)
    index = index + 1
    new.line =
       sprintf ('letter-probability matrix: alength= 4 w= %d nsites= %d E=8.1e-020',
          nrow (tfMat), 45, 8.1e-020)
    s [index] =  new.line
    index = index + 1
    for (r in 1:nrow (tfMat)) {
      new.row = sprintf (' %12.10f  %12.10f  %12.10f  %12.10f', tfMat [r,1],
                          tfMat [r,2], tfMat [r,3], tfMat [r,4])
      s [index] = new.row
      index = index + 1
      }
    s [index] = ''
    index = index + 1
    } # for name

  invisible (s)

} # matrixToMemeText
#-------------------------------------------------------------------------------
# connection is a character string, create a file by that name, open the file.
# dispatch to export which dispatches on con='connection'

setMethod ('export',  signature=c(object='MotifList', con='character',
                                  format='character'),

  function (object, con, format, ...) {
      ## do minimum work unique to this method, then dispatch to avoid
      ## code duplication
    con = file (con, 'w')
    on.exit(close(con))
    export(object, con, format, ...)
    })

#-------------------------------------------------------------------------------
# write to connection with specified format
# format includes TRANSFAC, meme (also good for tomtom), and tsv

setMethod ('export',  signature=c(object='MotifList', con='connection',
                                 format='character'),

  function (object, con, format, ...) {

    fmt = match.arg (tolower (format), c ('meme', 'transfac','jaspar'))
    ## match.arg fails if !fmt %in% c('meme', 'transfac'), so no need
    ## for test
    ## let the user manage opened cons
    if (!isOpen(con)) {
        open(con)
        on.exit(close(con))
    }
    if (fmt == 'meme') {
        text = matrixToMemeText (object)
    } else if (fmt == 'jaspar') {
        text = matrixToJasparText (object)
    }
    cat (text, sep='\n', file=con)
  })


#-------------------------------------------------------------------------------
# write to connection, using default format,  ??? for matrix list, tsv for
# metadata
setMethod ('export',  signature=c(object='MotifList',  con='missing',
                                 format='character'),

  function (object, con, format,  ...) {

    fmt = match.arg (tolower (format), c ('meme','jaspar')) # , 'transfac'
    if (fmt == 'meme') {
        text = paste (matrixToMemeText (object), collapse='\n')
      cat (text)
      invisible (text)
    } else if (fmt == 'jaspar') {
        text = paste (matrixToJasparText (object), collapse='\n')
      cat (text)
      invisible (text)
    }

  })

#-------------------------------------------------------------------------------
setMethod('show', 'MotifList',

    function(object) {
      msg = sprintf ('MotifDb object of length %d', length (object))
      cat (msg, '\n', sep='')
      if (length (object) == 0)
        return ()

      cat ('| Created from downloaded public sources: 2013-Aug-30', '\n', sep='')

      tbl.dataSource = as.data.frame (table (mcols (object)$dataSource))
      tbl.org = as.data.frame (table (mcols (object)$organism))
      tbl.org = head (tbl.org [order (tbl.org$Freq, decreasing=TRUE),])
      totalMatrixCount = length (object)
      totalOrganismCount = length (unique (mcols (object)$organism))
      dataSourceCount = nrow (tbl.dataSource)

      source.singular.or.plural = 'sources'
      if (dataSourceCount == 1)
        source.singular.or.plural = 'source'

      msg = sprintf ('| %d position frequency matrices from %d %s:',
                     totalMatrixCount, dataSourceCount, source.singular.or.plural)
      cat (msg, '\n', sep='')
      for (r in 1:nrow (tbl.dataSource)) {
         dataSource = tbl.dataSource$Var1 [r]
         matrixCount = tbl.dataSource$Freq [r]
         msg = sprintf ('| %18s: %4d', dataSource, matrixCount)
         cat (msg, '\n', sep='')
         } # for r
      msg = sprintf ('| %d organism/s', totalOrganismCount)
      cat (msg, '\n', sep='')
      for (r in 1:nrow (tbl.org)) {
         organism = tbl.org$Var1 [r]
         orgCount = tbl.org$Freq [r]
         msg = sprintf ('| %18s: %4d', organism, orgCount)
         cat (msg, '\n', sep='')
         } # for r
      otherOrgCount = totalMatrixCount - sum (tbl.org$Freq)
      if (otherOrgCount > 0) {
        msg = sprintf ('| %18s: %4d', 'other', otherOrgCount)
        cat (msg, '\n', sep='')
        }
    if (!is.null (names (object))) {
      all.names = names (object)
      count = length (all.names)
      if (count <= 10)
        cat (paste (all.names, '\n'), sep='')
      else {
        cat (paste (all.names [1:5], '\n'), sep='')
        cat ('...', '\n', sep='')
        end = length (all.names)
        start = end - 4
        cat (paste (all.names [start:end], '\n'), sep='')
        }
      }
    })
#-------------------------------------------------------------------------------
setMethod ('query', 'MotifList',

   function (object, queryString, ignore.case=TRUE) {
       indices = unique (as.integer (unlist (sapply (colnames (mcols (object)),
                    function (colname)
                       grep (queryString, mcols (object)[, colname],
                             ignore.case=ignore.case)))))
        object [indices]
      })
#-------------------------------------------------------------------------------
# Addition on 2017/06/15 from Matt Richards

# This will not exactly match JASPAR because units are PFM and JASPAR uses PCM
# General JASPAR Format:

# > "Motif Name"\t"Transcription Factor"
# A [ PCMS ]
# C [ PCMS ]
# G [ PCMS ]
# T [ PCMS ]
#
# ...

# Note: the PCMs are space-delimited

matrixToJasparText <- function (matrices)
{
  matrix.count <- length (matrices)

  # Incoming matrices have nucleotide rows, position columns.
  # This is the correct orientation for JASPAR; however, we need to also
  # add brackets and letters to them

  # Calculate the number of lines of text by counting matrices and assuming
  # 6 lines per matrix

  predicted.line.count <- 6*matrix.count

  #s = vector ('character', predicted.line.count)
  s <- character (predicted.line.count)

  index <- 1

  for (name in names (matrices)) {

      # Print the name with an arrow, follwed by the motif
      s[index] <- sprintf('>%s',name)
      index <- index + 1

      # For each line of the matrix, print the correct letter and the
      # matrix row surrounded by brackets
      motif.matrix <- matrices[name][[1]]


      for (r in 1:nrow(motif.matrix)) {
          s[index] <- sprintf("%s [ %s ]",
                              rownames(motif.matrix)[r],
                              paste(motif.matrix[r,],collapse=" "))
          index <- index + 1
      }

      s[index] <- ""
      index <- index + 1

  } # for name

  # Remove the last blank line
  s <- s[-length(s)]

  invisible (s)

} # matrixToJasparText
#-------------------------------------------------------------------------------
# returns a data.frame with motif, geneSymbol, source, pubmedID columns
setMethod ('motifToGene', 'MotifList',

   function (object, motifs, source) {
     source <- tolower(source)
     stopifnot(source %in% c("motifdb", "tfclass"))
     tbl <- data.frame()
     if(source %in% c("motifdb")){
        providerId <- NULL   # avoid R CMD check note
        tbl <- as.data.frame(subset(mcols(object), providerId %in% motifs))
        if(nrow(tbl) == 0)
           return(data.frame())
        tbl <- unique(tbl [, c("geneSymbol", "providerId", "dataSource", "organism", "pubmedID")])
        colnames(tbl) <- c("geneSymbol", "motif", "dataSource", "organism", "pubmedID")
        tbl <- tbl[, c("motif", "geneSymbol", "dataSource", "organism", "pubmedID")]
        if(nrow(tbl) > 0)
           tbl$source <- "MotifDb"
        }
     if(source %in% c("tfclass")){
        motif <- NULL
        tbl <- subset(object@manuallyCuratedGeneMotifAssociationTable, motif %in% motifs)
        if(nrow(tbl) == 0)
           return(data.frame())
        tbl <- unique(tbl[, c("motif", "tf.gene", "pubmedID")])
        tbl <- tbl[order(tbl$motif),]
        rownames(tbl) <- NULL
        colnames(tbl) <- c("motif", "geneSymbol", "pubmedID")
        if(nrow(tbl) > 0)
           tbl$source <- "TFClass"
        }
     tbl
     })

#-------------------------------------------------------------------------------
# returns a data.frame with motif, geneSymbol, source, pubmedID columns
setMethod ('geneToMotif', 'MotifList',

   function (object, geneSymbols, source, ignore.case=FALSE) {
     source <- tolower(source)
     stopifnot(source %in% c("motifdb", "tfclass"))
     extract.mdb <- function(gene){
        geneSymbol <- NULL # workaround the R CMD check "no visible binding for global variable"
        if(ignore.case)
           tbl <- as.data.frame(subset(mcols(object), tolower(geneSymbol) == tolower(gene)))
        else
           tbl <- as.data.frame(subset(mcols(object), geneSymbol == gene))

        tbl <- unique(tbl [, c("geneSymbol", "providerId", "dataSource", "organism", "pubmedID")])
        colnames(tbl) <- c("geneSymbol", "motif", "dataSource", "organism", "pubmedID")
        tbl
        }
     if(source %in% c("motifdb")){
        tbls <- lapply(geneSymbols, extract.mdb)
        result <- do.call(rbind, tbls)
        if(nrow(result) > 0)
           result$source <- "MotifDb"
        }
     if(source %in% c("tfclass")){
        if(ignore.case)
           tbl <- subset(object@manuallyCuratedGeneMotifAssociationTable, tolower(tf.gene) %in% tolower(geneSymbols))
        else
           tbl <- subset(object@manuallyCuratedGeneMotifAssociationTable, tf.gene %in% geneSymbols)
        tf.gene <- NULL; motif <- NULL  # workaround R CMD CHECK "no visible binding ..." bogus error
        tbl <- unique(tbl[, c("motif", "tf.gene", "pubmedID")])
        tbl <- tbl[order(tbl$tf.gene),]
        rownames(tbl) <- NULL
        colnames(tbl) <- c("motif", "geneSymbol", "pubmedID")
        result <- tbl[, c("geneSymbol", "motif", "pubmedID")]
        if(nrow(result) > 0)
           result$source <- "TFClass"
        }
     result
     })

#-------------------------------------------------------------------------------
setMethod('associateTranscriptionFactors', 'MotifList',

   function(object, tbl.withMotifs, source, expand.rows){
     source <- tolower(source)
     stopifnot(source %in% c("motifdb", "tfclass"))
     tbl.out <- data.frame()
     if(source %in% c("motifdb")){
           # lookup up in the object metadata, expect one TF geneSymbol per matrix name
        pfm.ids <- tbl.withMotifs[, "motifName"]
        matched.rows <- match(pfm.ids, names(as.list(object)))
        #if(length(matched.rows) == nrow(tbl.withMotifs)) {
        tbl.new <- mcols(object)[matched.rows, c("geneSymbol", "pubmedID")]
        tbl.new$geneSymbol[nchar(tbl.new$geneSymbol)==0] <- NA
        tbl.new$pubmedID[nchar(tbl.new$pubmedID)==0] <- NA
        tbl.out <- as.data.frame(cbind(tbl.withMotifs, tbl.new))
        } # direct
     if(source %in% c("tfclass")){
        if(! "shortMotif" %in% colnames(tbl.withMotifs)){
           stop("MotifDb::assoicateTranscriptionFactors needs a 'shortMotif' column with the TFClass source")
           }
        tbl.tfClass <- read.table(system.file(package="MotifDb", "extdata", "tfClass.tsv"), sep="\t", as.is=TRUE, header=TRUE)
        motif.ids <- tbl.withMotifs[, "shortMotif"]
        geneSymbols <- lapply(motif.ids, function(id)
                                 paste(tbl.tfClass$tf.gene[grep(id, tbl.tfClass$motif, fixed=TRUE)], collapse=";"))
        geneSymbols <- unlist(geneSymbols)
        pubmedIds   <- lapply(motif.ids, function(id)
                                 unique(tbl.tfClass$pubmedID[grep(id, tbl.tfClass$motif, fixed=TRUE)]))
        pubmedIds   <- as.character(pubmedIds)
        pubmedIds   <- gsub("integer(0)", "", pubmedIds, fixed=TRUE)
        tbl.new     <- data.frame(geneSymbol=geneSymbols, pubmedID=pubmedIds, stringsAsFactors=FALSE)
        tbl.new$geneSymbol[nchar(tbl.new$geneSymbol)==0] <- NA
        tbl.new$pubmedID[nchar(tbl.new$pubmedID)==0] <- NA
        tbl.out <- as.data.frame(cbind(tbl.withMotifs, tbl.new))

        if(expand.rows){
           rows.with.na <- which(is.na(tbl.out$geneSymbol))
           rows.with.geneSymbol <- setdiff(1:nrow(tbl.out), rows.with.na)
           tbl.asIs <- tbl.out[rows.with.na,]
           tbl.toExpand <- tbl.out[rows.with.geneSymbol,]
           geneSymbols.split <- strsplit(tbl.toExpand$geneSymbol, ";")
           counts <- unlist(lapply(geneSymbols.split, length))
           geneSymbols.split.vec <- unlist(geneSymbols.split)
           tbl.expanded <- splitstackshape::expandRows(tbl.toExpand, counts, count.is.col=FALSE, drop=FALSE)
           stopifnot(length(geneSymbols.split.vec) == nrow(tbl.expanded))
           tbl.expanded$geneSymbol <- geneSymbols.split.vec
           tbl.out <- rbind(tbl.expanded, tbl.asIs)
           }
        } # indirect
     tbl.out
     })

#-------------------------------------------------------------------------------
