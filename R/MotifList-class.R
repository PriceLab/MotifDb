setGeneric('query', signature='object', function(object, queryString, ignore.case=TRUE)
            standardGeneric ('query'))
#-------------------------------------------------------------------------------
setClass ('MotifList',
          contains='SimpleList',
          representation (elementMetadata='DataFrame'),
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
  
  object = new ('MotifList', listData=matrices, elementMetadata = df)

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

    fmt = match.arg (tolower (format), c ('meme', 'transfac'))
    ## match.arg fails if !fmt %in% c('meme', 'transfac'), so no need
    ## for test
    ## let the user manage opened cons
    if (!isOpen(con)) {
        open(con)
        on.exit(close(con))
    }
    if (fmt == 'meme') 
      text = matrixToMemeText (object)
    cat (text, sep='\n', file=con)
    })

#-------------------------------------------------------------------------------
# write to connection, using default format,  ??? for matrix list, tsv for
# metadata
setMethod ('export',  signature=c(object='MotifList',  con='missing',
                                 format='character'),

  function (object, con, format,  ...) {

    fmt = match.arg (tolower (format), c ('meme')) # , 'transfac'
    if (fmt == 'meme') {
      text = paste (matrixToMemeText (object), collapse='\n')
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

      tbl.dataSource = as.data.frame (table (values (object)$dataSource))
      tbl.org = as.data.frame (table (values (object)$organism))
      tbl.org = head (tbl.org [order (tbl.org$Freq, decreasing=TRUE),])
      totalMatrixCount = length (object)
      totalOrganismCount = length (unique (values (object)$organism))
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
       indices = unique (as.integer (unlist (sapply (colnames (values (object)), 
                    function (colname) 
                       grep (queryString, values (object)[, colname], 
                             ignore.case=ignore.case)))))
        object [indices]
      })
#-------------------------------------------------------------------------------
