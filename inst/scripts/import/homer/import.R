# MotifDb/inst/scripts/import/homer/import.R
#------------------------------------------------------------------------------------------------------------------------
options (stringsAsFactors=FALSE)
printf <- function(...) print(noquote(sprintf(...)))
library(RCurl)
#------------------------------------------------------------------------------------------------------------------------
run = function (dataDir)
{
  dataDir <- file.path(dataDir, "homer")
  rawMatrixList <- readRawMatrices (dataDir)
  matrices <- extractMatrices (rawMatrixList)
 browser() 
  tbl.md <- createMetadataTable(matrices) #(dataDir, matrices,
  #raw.metadata.filename="md-raw.tsv")
#  browser()
  matrices <- normalizeMatrices (matrices)
  matrices <- renameMatrices (matrices, tbl.md)
  
  serializedFile <- file.path(dataDir, "homer.RData")
  printf("writing %s to %s", "homer.RData", dataDir)
  
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
  
  filename <- file.path(dataDir, "custom.motifs.mod.txt")
  printf("checking for readable human matrix file:")
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
    
  #browser()
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
createMetadataTable = function (matrices)
{
  #  browser()
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
      id.pieces <- unlist(strsplit(matrix.id, "\\/"))
      # Piece 1 has TF; Piece 2 has Origin; Piece 3 has program/tool/resource
      tf <- sub("\\(.*\\)","", id.pieces[[1]])
     
      organism <- NA
    
      dataSource <- "HOMER"
      
    # split.matrix.name <- unlist(strsplit(short.matrix.name, "_"))[1]
    # shorter.matrix.name <- split.matrix.name
    #if (grepl(split.matrix.name, "+")){
    #  shorter.matrix.name <- unlist(strsplit(split.matrix.name, "+"))[1]
    #}
    
    #uri <- paste0("www.uniprot.org/uniprot/?query=",idStr)
    if (nchar(id.pieces[1]) <=9){#!("+" %in% shorter.matrix.name)
      idStr <- paste0(id.pieces[1], "_", id.pieces[2])
      protIDURL <- geturlname(idStr) #gets the URL for the proteinID from the geneSymbol
      protID <- unlist(strsplit(protIDURL, "http://www.uniprot.org/uniprot/"))[-1]
      }else{
        protID <- rep(NA,1)
      }
      
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
                    pubmedID="26586801")
    printf("matrix.id: %s", matrix.id);
    tbl.md = rbind (tbl.md, data.frame (new.row, stringsAsFactors=FALSE))
    full.name = sprintf ('%s-%s-%s', organism, dataSource, matrix.id)
    rownames (tbl.md) [nrow (tbl.md)] = full.name
  } # for matrix.id
  browser()
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
  lines = strsplit (text, '\t')
  # Remove the arrow from sequence; save it and title separately
  consensus.sequence <- sub(">", "", lines [[1]][1])
  title <- lines[[1]][2]
  line.count = length(lines)
  
  # browser()
  result = matrix (nrow=line.count-1, ncol=4, dimnames=list(1:(line.count-1), c ('A','C','G','T')))
  row = 1
  for (line in lines [2:line.count]) {
    result [row,] = as.numeric (line)
    row = row + 1
  } # for line

  result <- t(result)
  
  return (list (title=title, consensus.sequence=consensus.sequence, matrix=result))
  
} # parsePwm
#----------------------------------------------------------------------------------------------------
