# MotifDb/inst/scripts/import/SwissRegulon/import.R
#------------------------------------------------------------------------------------------------------------------------
options (stringsAsFactors=FALSE)
printf <- function(...) print(noquote(sprintf(...)))
library(RCurl)
#------------------------------------------------------------------------------------------------------------------------
run = function (dataDir)
{
  dataDir <- file.path(dataDir, "SwissRegulon")
  rawMatrixList <- readRawMatrices (dataDir)
  matrices <- extractMatrices (rawMatrixList)
 
  tbl.md <- createMetadataTable(matrices) #(dataDir, matrices,

  matrices <- normalizeMatrices (matrices)
  matrices <- renameMatrices (matrices, tbl.md)
  
  serializedFile <- file.path(dataDir, "SwissRegulon.RData")
  printf("writing %s to %s", "SwissRegulon.RData", dataDir)
  
  save (matrices, tbl.md, file=serializedFile)
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
  
  filename <- file.path(dataDir, "swiss_regulon_hg19_f5.txt")
  printf("checking for readable human matrix file:")
  printf("     %s", filename)
  stopifnot(file.exists(filename))
  
    all.lines = scan (filename, what=character(0), sep='\n', quiet=TRUE)
    # Remove lines with "//"
    all.lines <- grep("^[^//]",all.lines,value = TRUE)
    # Grab title lines
  title.lines <- grep ('^NA', all.lines)
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
    
    #Add the final matrix
    start.line <- title.lines[title.line.count]
    end.line <- length(all.lines)
    new.pwm = parsePwm (all.lines [start.line:end.line])
    pwms = c (pwms, list (new.pwm))
    
  invisible (pwms)
  
} # readRawMatrices
#------------------------------------------------------------------------------------------------------------------------
extractMatrices = function (pwm.list)
{
  matrices = sapply (pwm.list, function (element) element$matrix)
  matrix.names <- sapply (pwm.list, function (element) element$title)
  # Add the "SwissRegulon" suffix to all motifs
  matrix.names <- paste(matrix.names,"SwissRegulon",sep=".")
  #browser()
  names (matrices) <- matrix.names
  
  matrices
  
} # extractMatrices
#------------------------------------------------------------------------------------------------------------------------
createMetadataTable = function (matrices)
{
    #browser()
    tbl.md = data.frame ()
    matrix.ids = names(matrices)
  geturlname <- function(name){
    h = getCurlHandle()
    z <- getURL(paste0("www.uniprot.org/uniprot/?query=",name),
                followlocation=TRUE, curl=h)
    getCurlInfo(h)$effective.url # catch the url redirect
  }
    # Assume we have either Human or Mouse
  for (matrix.id in matrix.ids) {
      my.matrix <- matrices[[matrix.id]]
      # Split up the ID pieces (symbol, organism, database, version)
      #    short.matrix.name <- sub("\\..*$", "", matrix.id)
      id.pieces <- unlist(strsplit(matrix.id, "\\."))
      # Piece 1 has TF; Piece 2 has Origin; Piece 3 has program/tool/resource
      tf <- id.pieces[1]

      # These are all human
      organism <- "Hsapiens"
    
      dataSource <- "SwissRegulon"
      
    # split.matrix.name <- unlist(strsplit(short.matrix.name, "_"))[1]
    # shorter.matrix.name <- split.matrix.name
    #if (grepl(split.matrix.name, "+")){
    #  shorter.matrix.name <- unlist(strsplit(split.matrix.name, "+"))[1]
    #}
    
    #uri <- paste0("www.uniprot.org/uniprot/?query=",idStr)
    protID <- rep(NA,1)
      
    new.row = list (providerName=matrix.id,
                    providerId=matrix.id, 
                    dataSource= dataSource,
                    geneSymbol= tf, #md$symbol
                    geneId= NA,
                    geneIdType= NA,
                    proteinId=protID,
                    proteinIdType="UNIPROT",
                    organism= organism,
                    sequenceCount=max(colSums(my.matrix)),
                    bindingSequence=NA_character_,
                    bindingDomain=NA,
                    tfFamily=NA, #family
                    experimentType="low- and high-throughput methods",
                    pubmedID="19377474")
      printf("matrix.id: %s", matrix.id)      
      tbl.md = rbind (tbl.md, data.frame (new.row, stringsAsFactors=FALSE))      
      full.name = sprintf ('%s-%s-%s', organism, dataSource, matrix.id)      
      rownames (tbl.md) [nrow (tbl.md)] = full.name
      
  } # for matrix.id
 
  invisible (tbl.md)
  
} # createMetadataTable
#------------------------------------------------------------------------------------------------------------------------
renameMatrices = function (matrices, tbl.md)
{
  stopifnot (length (matrices) == nrow (tbl.md))
  names (matrices) = rownames (tbl.md)
  invisible (matrices)
  
} # renameMatrices
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
  lines <- strsplit (text, '\\s+')
  # Grab the TF title
  title <- lines[[1]][2]
  line.count = length(lines)
  
  # browser()
  result = matrix (nrow=line.count-2, ncol=4, dimnames=list(1:(line.count-2), c ('A','C','G','T')))
  row = 1
  # Loop through lines, ignoring the 1st (title) and 2nd (labels)
  for (line in lines [3:line.count]) {
      # Grab only the last 4 elements, stripping the data record number
    result [row,] = as.numeric (line[2:5])
    row = row + 1
  } # for line

  result <- t(result)
  
  return (list (title=title, matrix=result))
  
} # parsePwm
#----------------------------------------------------------------------------------------------------
