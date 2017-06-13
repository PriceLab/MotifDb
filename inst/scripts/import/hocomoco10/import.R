# MotifDb/inst/scripes/import/HOCOMOCO/import.R
#------------------------------------------------------------------------------------------------------------------------
options (stringsAsFactors=FALSE)
printf <- function(...) print(noquote(sprintf(...)))
library(RCurl)
#------------------------------------------------------------------------------------------------------------------------
run = function (dataDir)
{
  dataDir <- file.path(dataDir, "hocomoco10")
  rawMatrixList <- readRawMatrices (dataDir)
  matrices <- extractMatrices (rawMatrixList)
  
  tbl.md <- createMetadataTable(matrices) #(dataDir, matrices,
  #raw.metadata.filename="md-raw.tsv")
  matrices <- normalizeMatrices (matrices)
  matrices <- renameMatrices (matrices, tbl.md)
  
  serializedFile <- file.path(dataDir, "hocomoco10.RData")
  printf("writing %s to %s", "hocomoco10.RData", dataDir)
  
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
  
  filename <- file.path(dataDir, "HOCOMOCOv10_HUMAN_mono_jaspar_format.txt") #old filename: "hoco.pcm"
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

    start.line <- title.lines[title.line.count]
    end.line <- length(all.lines)
    new.pwm = parsePwm (all.lines [start.line:end.line])
    pwms = c (pwms, list (new.pwm))

    # Do the same for the mouse file
    filename <- file.path(dataDir, "HOCOMOCOv10_MOUSE_mono_jaspar_format.txt") #old filename: "hoco.pcm"
    printf("checking for readable mouse matrix file:")    
    printf("     %s", filename)    
    stopifnot(file.exists(filename))    
  
    all.lines = scan (filename, what=character(0), sep='\n', quiet=TRUE)    
    title.lines = grep ('^>', all.lines)    
    title.line.count <<- length (title.lines)   
    max = title.line.count - 1
  
  #loops through all motifs in the matrix file, one motif at a time
  for (i in 1:max) {
    start.line = title.lines [i]
    end.line = title.lines [i+1] - 1
    new.pwm = parsePwm (all.lines [start.line:end.line])
    pwms = c (pwms, list (new.pwm))
  } # for i

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
  matrix.names <- sub("^> ", "", matrix.names)
  names (matrices) <- matrix.names
  
  matrices
  
} # extractMatrices
#------------------------------------------------------------------------------------------------------------------------
createMetadataTable = function (matrices)
{
    # There is no file to import; create one     
    # filename <- file.path(dataDir, "hocomoco", "md-raw.tsv")
 # printf("checking for readable metadata file:")
 # printf("   %s", filename)
 # stopifnot(file.exists(filename)) 
#  tbl.raw <- read.table(filename, sep="\t", header=TRUE, as.is=TRUE)

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
      id.pieces <- unlist(strsplit(matrix.id, "_|\\."))
     
    # Check the organism and assign it
    if(grepl("HUMAN",matrix.id)){
        organism <- "Hsapiens"
    } else { organism <- "Mmusculus"
    }
    
      dataSource <- "HOCOMOCOv10"
      
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
                    providerId=matrix.id, #"HOCOMOCO v10"
                    dataSource= dataSource,
                    geneSymbol= id.pieces[1], #md$symbol
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

  # Remove the arrow from title
  title <- sub(">", "", lines [[1]][1])
  line.count = length(lines)
  
  # expect 4 rows, and a number of columns we can discern from
  # the incoming text.
  # browser()
  #secondLineParsed <- strsplit(lines[[2]], " ")[[1]]
  secondLineParsed <- lines[[2]]
  class(secondLineParsed) <- "numeric"
  cols <- length(secondLineParsed)
  result <- matrix (nrow=4, ncol=cols,
                    dimnames=list(c('A','C','G','T'),
                                  as.character(1:cols)))
  # loop over the four lines (for each base respectively)
  row = 1
  
  for(i in 2:line.count){
    #linesParsed <- strsplit(lines[[i]], " ")[[1]]

      linesParsed <- lines[[i]]
      class(linesParsed) <- "numeric"
    result [row,] = as.numeric (linesParsed)
    row = row + 1
  } # for i
  
  #return (list (title=title, consensus.sequence=consensus.sequence, matrix=result))
  return (list (title=title, matrix=result))
  
} # parsePwm
#----------------------------------------------------------------------------------------------------
