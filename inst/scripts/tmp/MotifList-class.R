setClass ('MotifList',
          contains='SimpleList',
          representation (elementMetadata="DataFrame"),
          prototype(elementType="matrix")
          )

setValidity("MotifList", function(object) {
    msg <- NULL

    ## what makes for a valid MotifList?
    if (is.null(names(object)))
        msg <- c(msg, "'names()' must not be NULL")
    else if (any(duplicated(names(object))))
        msg <- c(msg, "all 'names()' must be unique")
    ## all(sapply(object, class) == "matrix") works for 0-length object
    if (!all(sapply(object, class) == "matrix"))
        msg <- c(msg, "all elements must be of class 'matrix'")
    if (!all(sapply(object, nrow) == 4))
        msg <- c(msg, "all elements must have 4 rows")
    ok <- sapply(object, function(elt) all(colSums(elt) == 1))
    if (!all(ok))
        msg <- c(msg, "all elements must have colSums equal to 1")

    if (is.null(msg)) TRUE else msg
})

setGeneric ('add', function (obj, matrices, tbl.metadata) standardGeneric ('add'),
            signature='obj')
#------------------------------------------------------------------------------------------------------------------------
MotifList = function (matrices=list(),
                      providerNames=character(),              # badis.ABF2, Cell08/Rhox11_2205.1_pwm.txt
                      providerIds=character(),              # badis.ABF2, Cell08/Rhox11_2205.1_pwm.txt
                      geneSymbols=character(),              # ABF2, Rhox11
                      geneIds=character(),
                      geneIdTypes=character(),
                      proteinIds=character(),
                      proteinIdTypes=character(),
                      sequenceCounts=integer(),             # often NA
                      organisms=character(),                # Scerevisiae, Mmusculus
                      bindingDomains=character(),       # NA, Homeo
                      bindingSequences=character(),
                      dataSources=character(),              # ScerTF, UniPROBE
                      experimentTypes=character(),          # NA, protein-binding microarray
                      pubmedIDs=character(),                # 19111667, 1858359
                      tfFamilies=character())
{  
  ## these are redundant with DataFrame validity method
  elementMetadata = DataFrame (#name=names,
                               providerName=providerNames,
                               providerId=providerIds,
                               geneSymbol=geneSymbols,
                               geneId=geneIds,
                               geneIdType=geneIdTypes,
                               proteinId=proteinIds,
                               proteinIdType=proteinIdTypes,
                               sequenceCount=sequenceCounts,
                               organism=organisms,
                               bindingDomain=bindingDomains,
                               bindingSequence=bindingSequences,
                               dataSource=dataSources,
                               experimentType=experimentTypes,
                               pubmedID=pubmedIDs,
                               tfFamily=tfFamilies)
  new ('MotifList', listData = matrices, elementMetadata=elementMetadata)
  
} # ctor
#------------------------------------------------------------------------------------------------------------------------
## quiet vs. verbose
MotifDb = function (loadAllSources=TRUE, quiet=TRUE)
{
  mdb <- MotifList ()
  if (loadAllSources) {
    data.path = system.file(package="MotifDb", "data")
    data.files = dir (data.path, full.names=TRUE)
    sources.to.exclude = "ScerTF.RData"   # pending minor revisions
    data.files = data.files[!basename(data.files) %in% sources.to.exclude]
      for (data.file in data.files) {
        ## see globalVariables(tbl.md, matrices)
        tbl.md = NA; matrices = NA;  # define these to keep 'check' happy.  they are loaded by 'load'
        load (data.file)
        mdb <- append(mdb, do.call("MotifList", c(matrices, tbl.md)))
        if (!quiet)
          message(sprintf ('added %s (%d) matrices, length now: %d',
                           basename(data.file), length (matrices),
                           length (mdb)))
      } # for data.file
  
    if (!quiet)
      print (table (mdb$dataSource))
    } # if loadAllSources

  mdb
} # MotifDb
#------------------------------------------------------------------------------------------------------------------------
setMethod ('subset', signature = 'MotifList',

  function (x, subset, select, drop = FALSE, ...) {

    if (missing (subset)) 
      i <- TRUE
    else {
      i <- eval(substitute (subset), elementMetadata (x), parent.frame ())
      i <- try(as.logical(i), silent=TRUE)
      if (inherits(i, "try-error"))
        stop("'subset' must be coercible to logical")
      i <- i & !is.na(i)
      } # else
    x [i]
  })

#------------------------------------------------------------------------------------------------------------------------
## 'add' is too general a name; implement as do.call + append, but
## then not much added value so omit method?
setMethod ('add', signature = 'MotifList',

  function (obj, matrices, tbl.metadata) { 

     if (!is.list(matrices))            # allow solo matrix
         matrices <- list(matrices)
     ## 1:1 mapping between tbl.metadata (data.frame?) column names
     ## and MotifDb elementMetadata
     new.list <- do.call("MotifList", c(matrices, tbl.metadata))
     append(obj, new.list)
    })

#------------------------------------------------------------------------------------------------------------------------
setMethod("$", signature="MotifList",   # is this a good idea?

  function(x, name) {
    elementMetadata (x)[[name]]
    })
#------------------------------------------------------------------------------------------------------------------------
# transpose 4-row matrices to 4-column, so that there are as many rows as there are nucleotides in the motif sequence
# normalize the values, so that at each position the frequency of the 4 bases sums to 1
#
# for example, this count matrix for a sequence 10 bases long
#
#    1  2  3  4  5  6  7  8  9 10
# A  0  3 79 40 66 48 65 11 65  0
# C 94 75  4  3  1  2  5  2  3  3
# G  1  0  3  4  1  0  5  3 28 88
# T  2 19 11 50 29 47 22 81  1  6
#
# becomes this
#            1          2          3          4          5          6          7          8          9         10
# A 0.00000000 0.03092784 0.81443299 0.41237113 0.68041237 0.49484536 0.67010309 0.11340206 0.67010309 0.00000000
# C 0.96907216 0.77319588 0.04123711 0.03092784 0.01030928 0.02061856 0.05154639 0.02061856 0.03092784 0.03092784
# G 0.01030928 0.00000000 0.03092784 0.04123711 0.01030928 0.00000000 0.05154639 0.03092784 0.28865979 0.90721649
# T 0.02061856 0.19587629 0.11340206 0.51546392 0.29896907 0.48453608 0.22680412 0.83505155 0.01030928 0.06185567
#
# and then this
# 
#             A          C          G          T
# 1  0.00000000 0.96907216 0.01030928 0.02061856
# 2  0.03092784 0.77319588 0.00000000 0.19587629
# 3  0.81443299 0.04123711 0.03092784 0.11340206
# 4  0.41237113 0.03092784 0.04123711 0.51546392
# 5  0.68041237 0.01030928 0.01030928 0.29896907
# 6  0.49484536 0.02061856 0.00000000 0.48453608
# 7  0.67010309 0.05154639 0.05154639 0.22680412
# 8  0.11340206 0.02061856 0.03092784 0.83505155
# 9  0.67010309 0.03092784 0.28865979 0.01030928
# 10 0.00000000 0.03092784 0.90721649 0.06185567
transformMatrixToMemeRepresentation = function (m)
{
  stopifnot (class (m) == 'matrix')
  stopifnot (rownames (m) == c ('A', 'C', 'G', 'T'))
  col.sums = round (colSums (m))
  stopifnot (all (col.sums == 1))
  return (t (m))

} # transformMatrixToMemeRepresentation
#------------------------------------------------------------------------------------------------------------------------
matrixToMemeText = function (matrices)
{
  matrix.count = length (matrices)

  ## incoming matrices have nucleotide rows, position columns.  meme
  ## format, however, requires position-rows, and nucleotide-columns
  ## calculate the number of lines of text by counting columns
  #printf ("length of pre-allocated vector: %d", length (s))

  header = c(
    'MEME version 4',
    '',
    'ALPHABET= ACGT',
    '',
    'strands: + -',
    '',
    'Background letter frequencies',
    'A 0.250 C 0.250 G 0.250 T 0.250 ',
    '')

  tfStrings <- Map(function(name, tfMat) {
    tfMat = transformMatrixToMemeRepresentation (tfMat)
    ## meme output may be used by tomtom, which uses matrix names as
    ## part of image filenames.  removed all file-system-unfriendly
    ## characters here
    fixed.name = gsub ('\\/', '_', name)
    id <- sprintf ('MOTIF %s', fixed.name)
    desc <-
        sprintf ('letter-probability matrix: alength= 4 w= %d nsites= 45 E=8.1e-020',
                 nrow (tfMat))
    mat <- sprintf (' %12.10f  %12.10f  %12.10f  %12.10f',
                    tfMat[,1], tfMat[,2], tfMat[,3], tfMat[,4])
    c(id, desc, mat, "")
  }, names(matricies), matricies)

  invisible(c(header, unlist(tfStrings)))
} # matrixToMemeText
#------------------------------------------------------------------------------------------------------------------------
# write to connection with specified format
# if connection is a character string, create a file by that name, open the file.
# if it is an already-open connection, use that directly
# format includes TRANSFAC, meme (also good for tomtom), and tsv

setMethod ('export',  signature=c(object='MotifList', con='character', format='character'),

  function (object, con, format, ...) {
      ## do minimum work unique to this method, then dispatch to avoid
      ## code duplication
    con = file (con, 'w')
    on.exit(close(con))
    export(object, con, format, ...)
    })

#------------------------------------------------------------------------------------------------------------------------
# write to connection with specified format
# if connection is a character string, create a file by that name, open the file.
# if it is an already-open connection, use that directly
# format includes TRANSFAC, meme (also good for tomtom), and tsv

setMethod ('export',  signature=c(object='MotifList', con='connection', format='character'),

  function (object, con, format, ...) {

    fmt = match.arg (tolower (format), c ('meme', 'transfac'))
    ## match.arg fails if !fmt %in% c('meme', 'transfac'), so no need
    ## for test
    ## let the user manage opened cons
    if (!isOpen(con)) {
        open(con)
        on.exit(close(con))
    }
    cat (text, sep='\n', file=con)
    })

#------------------------------------------------------------------------------------------------------------------------
# write to connection, using default format,  ??? for matrix list, tsv for metadata
## already handled by rtracklayer?
