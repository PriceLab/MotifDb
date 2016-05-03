# MotifDb/inst/scripes/import/HOCOMOCO/import.R
#------------------------------------------------------------------------------------------------------------------------
options (stringsAsFactors=FALSE)
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
run = function (dataDir)
{
  dataDir <- file.path(dataDir) # Comment:Tanja: Taken out after dataDir: "HOCOMOCO"
  rawMatrixList <- readRawMatrices (dataDir)
  matrices <- extractMatrices (rawMatrixList)
  
  #TODO change tsv
  #tbl.md <- createMetadataTable (dataDir, matrices,
                                 #raw.metadata.filename="md-raw.tsv")
  matrices <- normalizeMatrices (matrices)
  #matrices <- renameMatrices (matrices, tbl.md)
  
  serializedFile <- file.path(dataDir, "HOCOMOCO.RData")
  printf("writing %s to %s", "HOCOMOCO.RData", dataDir)
  
  save (matrices, file=serializedFile) # TODO reinsert "tbl.md," without the quotation marks before "file =..."
  printf("saved %d matrices to %s", length(matrices), serializedFile)
  printf("next step:  copy %s to <packageRoot>/MotifDb/inst/extdata, rebuild package", serializedFile)
  
} # run

#------------------------------------------------------------------------------------------------------------------------
readRawMatrices = function (dataDir)
{
  # our convention is that there is a shared "dataDir" visible to
  # the importer, and that within that directory there is one
  # subdirectory for each data source.
  # for this example importer, that directory will be <dataDir>/test
  # within which we will look for one small file "sample.pcm"
  
  
  filename <- file.path(dataDir, "HOCOMOCO", "hoco.pcm")
  printf("checking for readable matrix file:")
  printf("     %s", filename)
  stopifnot(file.exists(filename))
  
  all.lines = scan (filename, what=character(0), sep='\n', quiet=TRUE)
  title.lines = grep ('^>', all.lines)
  title.line.count <<- length (title.lines)
  max = title.line.count - 1
  
  pwms = list ()
  #loops through all motifs in the matrix file, one motif at a time
  for (i in 1:max) {
    start.line = title.lines [i]
    end.line = title.lines [i+1] - 1
    new.pwm = parsePwm (all.lines [start.line:end.line])
    pwms = c (pwms, list (new.pwm))
  } # for i
  
  
  invisible (pwms)
  
} # readRawMatrices
#------------------------------------------------------------------------------------------------------------------------
extractMatrices = function (pwm.list)
{
  matrices = sapply (pwm.list, function (element) element$matrix)
  matrix.names <- sapply (pwm.list, function (element) element$title)
  matrix.names <- sub("^> ", "", matrix.names)
  names (matrices) <- matrix.names
  
  matrices
  
} # extractMatrices
#------------------------------------------------------------------------------------------------------------------------
renameMatrices = function (matrices, tbl.md)
{
  stopifnot (length (matrices) == nrow (tbl.md))
  names (matrices) = rownames (tbl.md)
  invisible (matrices)
  
} # renameMatrices
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
normalizeMatrices = function (matrices)
{
  mtx.normalized = sapply (matrices,
                           function (mtx) apply (mtx, 2, function (colvector) colvector / sum (colvector)))
  
  invisible (mtx.normalized)
  
} # normalizeMatrices
#------------------------------------------------------------------------------------------------------------------------
parsePwm = function (text)
{
  lines = strsplit (text, '\t')
  #browser()
  stopifnot(length(lines)==5) # title line, one line for each base
  title = lines [[1]][1]
  line.count = length(lines)
  
  # expect 4 rows, and a number of columns we can discern from
  # the incoming text.
  secondLineParsed <- strsplit(lines[[2]], " ")[[1]]
  class(secondLineParsed) <- "numeric"
  cols <- length(secondLineParsed)
  result <- matrix (nrow=4, ncol=cols,
                    dimnames=list(c('A','C','G','T'),
                                  as.character(1:cols)))
  # loop over the four lines (for each base respectively)
  row = 1

  for(i in 2:line.count){
    linesParsed <- strsplit(lines[[i]], " ")[[1]]
    class(linesParsed) <- "numeric"
    result [row,] = as.numeric (linesParsed)
    row = row + 1
  } # for i
  
  #return (list (title=title, consensus.sequence=consensus.sequence, matrix=result))
  return (list (title=title, matrix=result))
  
} # parsePwm
#----------------------------------------------------------------------------------------------------

