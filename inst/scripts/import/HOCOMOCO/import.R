# MotifDb/inst/scripes/import/HOCOMOCO/import.R
#------------------------------------------------------------------------------------------------------------------------
options (stringsAsFactors=FALSE)
printf <- function(...) print(noquote(sprintf(...)))
library(RCurl)
#------------------------------------------------------------------------------------------------------------------------
run = function (dataDir)
{
  dataDir <- file.path(dataDir)
  rawMatrixList <- readRawMatrices (dataDir)
  matrices <- extractMatrices (rawMatrixList)
  
  tbl.md <- createMetadataTable (dataDir, matrices,
                                 raw.metadata.filename="md-raw.tsv")
  matrices <- normalizeMatrices (matrices)
  matrices <- renameMatrices (matrices, tbl.md)
  
  serializedFile <- file.path(dataDir, "HOCOMOCO.RData")
  printf("writing %s to %s", "HOCOMOCO.RData", dataDir)
  
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
  
  
  filename <- file.path(dataDir, "HOCOMOCO", "HOCOMOCOv9_AD_PLAINTEXT_H_WPCM.txt") #old filename: "hoco.pcm"
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
createMetadataTable = function (dataDir, matrices, raw.metadata.filename)
{
  filename <- file.path(dataDir, "HOCOMOCO", "md-raw.tsv")
  printf("checking for readable metadata file:")
  printf("   %s", filename)
  stopifnot(file.exists(filename))
  
  tbl.raw <- read.table(filename, sep="\t", header=TRUE, as.is=TRUE)
  tbl.md = data.frame ()
  matrix.ids = names(matrices)
  geturlname <- function(name){
    h = getCurlHandle()
    z <- getURL(paste0("www.uniprot.org/uniprot/?query=",name),
                followlocation=TRUE, curl=h)
    getCurlInfo(h)$effective.url # catch the url redirect
  }
  for (matrix.id in matrix.ids) {
    matrix <- matrices[[matrix.id]]
    short.matrix.name <- sub("\\..*$", "", matrix.id)
    #stopifnot(length(grep(short.matrix.name, tbl.raw$symbol)) == 1)
    #md <- as.list(subset(tbl.raw, symbol==short.matrix.name))
    dataSource <- "HOCOMOCOv9_AD_PLAINTEXT_H_PWM_hg19"
    organism <- "Hsapiens"
    
    split.matrix.name <- unlist(strsplit(short.matrix.name, "_"))[1]
    shorter.matrix.name <- split.matrix.name
    #if (grepl(split.matrix.name, "+")){
    #  shorter.matrix.name <- unlist(strsplit(split.matrix.name, "+"))[1]
    #}
    
    #uri <- paste0("www.uniprot.org/uniprot/?query=",idStr)
    if (nchar(short.matrix.name) <=9){#!("+" %in% shorter.matrix.name)
      idStr <- paste0(shorter.matrix.name, "_HUMAN")
      protIDURL <- geturlname(idStr) #gets the URL for the proteinID from the geneSymbol
      protID <- unlist(strsplit(protIDURL, "http://www.uniprot.org/uniprot/"))[-1]
      }else{
        protID <- rep(NA,1)
      }
      
    new.row = list (providerName=matrix.id,
                    providerId=matrix.id, #"HOCOMOCO v8 and ChiPMunk 3.1"
                    dataSource="HOCOMOCOv9_AD_PLAINTEXT_H_PWM_hg19",
                    geneSymbol=shorter.matrix.name, #md$symbol
                    geneId="9606",
                    geneIdType="ENTREZ",
                    proteinId=protID,
                    proteinIdType="UNIPROT",
                    organism="Hsapiens",
                    sequenceCount=max(colSums(matrix)),
                    bindingSequence=NA_character_,
                    bindingDomain=NA,
                    tfFamily=NA, #family
                    experimentType="low- and high-throughput methods",
                    pubmedID="23175603")
    printf("matrix.id: %s", matrix.id);
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
  lines = strsplit (text, '\t')
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
