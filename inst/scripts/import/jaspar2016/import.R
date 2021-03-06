# MotifDb/inst/scripes/import/demo/import.R
#------------------------------------------------------------------------------------------------------------------------
options (stringsAsFactors=FALSE)
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
run = function (dataDir)
{
  dataDir <- file.path(dataDir, "jaspar2016")
  rawMatrixList <- readRawMatrices (dataDir)
  matrices <- extractMatrices (rawMatrixList)
  tbl.md <- createMetadataTable (dataDir, matrices,
                                 raw.metadata.filename="md-raw.tsv")
  matrices <- normalizeMatrices (matrices)
  matrices <- renameMatrices (matrices, tbl.md)
  
  serializedFile <- file.path(dataDir, "jaspar2016.RData")
  printf("writing %s to %s", "jaspar2016.RData", dataDir)

  save (matrices, tbl.md, file=serializedFile)
  printf("saved %d matrices to %s", length(matrices), serializedFile)
  printf("next step:  copy %s to <packageRoot>/MotifDb/inst/extdata, rebuild package",
         serializedFile)

} # run
#------------------------------------------------------------------------------------------------------------------------
readRawMatrices = function (dataDir)
{
    # our convention is that there is a shared "dataDir" visible to
    # the importer, and that within that directory there is one
    # subdirectory for each data source.
    # for this example importer, that directory will be <dataDir>/test
    # within which we will look for one small file "sample.pcm"
    

  filename <- file.path(dataDir, "sample.pcm")
  printf("checking for readable matrix file:")
  printf("     %s", filename)
  stopifnot(file.exists(filename))
  
  all.lines = scan (filename, what=character(0), sep='\n', quiet=TRUE)
  title.lines = grep ('^>', all.lines)
  title.line.count <<- length (title.lines)
    max = title.line.count -1

    # browser()
    # Matt's additions to fix brackets and delimiters
    all.lines <- gsub("[A,C,G,T]\\s+\\[\\s*","",all.lines)
    all.lines <- gsub("\\s*\\]","",all.lines)    
    all.lines <- gsub("\\s+","\t",all.lines)    

  pwms = list ()
  
    for (i in 1:max) {
       # print(i)
    start.line = title.lines [i]
    end.line = title.lines [i+1] - 1
    new.pwm = parsePwm (all.lines [start.line:end.line])
    pwms = c (pwms, list (new.pwm))
    } # for i

    # Add the last motif, which got missed
    start.line <- title.lines[title.line.count]
    end.line <- start.line + 4
    new.pwm = parsePwm (all.lines [start.line:end.line])
    pwms = c (pwms, list (new.pwm))
    
  
  invisible (pwms)

} # readRawMatrices
#------------------------------------------------------------------------------------------------------------------------
extractMatrices = function (pwm.list)
{
  matrices = sapply (pwm.list, function (element) element$matrix)
  matrix.names <- sapply (pwm.list, function (element) element$title)
  matrix.names <- sub("^>", "", matrix.names)
  names (matrices) <- matrix.names
  
  matrices

} # extractMatrices
#------------------------------------------------------------------------------------------------------------------------
createMetadataTable = function (dataDir, matrices, raw.metadata.filename)
{
  filename <- file.path(dataDir, "md-raw.tsv")
  printf("checking for readable metadata file:")
  printf("   %s", filename)
  stopifnot(file.exists(filename))

  # browser()
  tbl.raw <- read.table(filename, sep="\t", header=TRUE)
  tbl.md = data.frame ()
  matrix.ids = names(matrices)
  
  for (matrix.id in matrix.ids) {
    matrix <- matrices[[matrix.id]]
    #short.matrix.name <- sub("\\..*$", "", matrix.id)
    #stopifnot(length(grep(short.matrix.name, tbl.raw$ma.name)) == 1)
    # Instead of chopping off version, add it to ma.name and subset based on that
    id.pieces <- strsplit(matrix.id, split="\\.")
    md <- as.list(subset(tbl.raw,
                         ma.name==id.pieces[[1]][1] & unknown == id.pieces[[1]][2]))
    print(matrix.id)
    #browser()
    stopifnot(length(md$id) == 1)
    dataSource <- "jaspar2016"
    organism <- ncbiTaxonimicCodeToBiocLinnaean(md$ncbi.tax.code)
    new.row = list (providerName=matrix.id,
                    providerId=matrix.id,
                    dataSource=dataSource,
                    geneSymbol=md$gene.symbol,
                    geneId=NA,
                    geneIdType=NA,
                    proteinId=md$uniprot,
                    proteinIdType=guessProteinIdentifierType(md$uniprot),
                    organism=organism,
                    sequenceCount=max(colSums(matrix)),
                    bindingSequence=NA_character_,
                    bindingDomain=NA,
                    tfFamily=md$family,
                    experimentType=md$type,
                    pubmedID="24194598")
    tbl.md = rbind (tbl.md, data.frame (new.row, stringsAsFactors=FALSE))
    full.name = sprintf ('%s-%s-%s-%s',
                         organism,
                         dataSource,
                         md$gene.symbol,
                         matrix.id)
    rownames (tbl.md) [nrow (tbl.md)] = full.name
    } # for i

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
# an empirical and not altogether trustworthy solution to identifying identifier types.
guessProteinIdentifierType = function (moleculeName)
{
      if (is.na (moleculeName))
          return (NA_character_)      
  if (nchar (moleculeName) == 0)
    return (NA_character_)


  first.char = substr (moleculeName, 1, 1)

  if (first.char == 'Y')
    return ('SGD')

  if (first.char %in% c ('P', 'Q', 'O', 'A', 'C'))
    return ('UNIPROT')

  if (length (grep ('^NP_', moleculeName)) == 1)
    return ('REFSEQ')

   return (NA_character_)

} # guessProteinIdentifierType
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
  cols <- length(lines[[2]])
  result <- matrix (nrow=4, ncol=cols,
                   dimnames=list(c('A','C','G','T'),
                                 as.character(1:cols)))
  row = 1
  for(i in 2:line.count){
    result [row,] = as.numeric (lines[[i]])
    row = row + 1
    } # for i

  #result = t (result)
    
  #return (list (title=title, consensus.sequence=consensus.sequence, matrix=result))
  return (list (title=title, matrix=result))

} # parsePwm
#----------------------------------------------------------------------------------------------------
ncbiTaxonimicCodeToBiocLinnaean <- function(code)
{
  code <- as.character(code)
  
  lookup <- list("3702" = "Athaliana",
                 "3888" = "Psativum",
                 "4094" = "Nsp.",
                 "4102" = "Phybrida",
                 "4151" = "Amajus",
                 "4513" = "Hvulgare",
                 "4565" = "Taestivam",
                 "4577" = "Zmays",
                 "4932" = "Scerevisiae",
                 "6239" = "Celegans",
                 "7227" = "Dmelanogaster",
                 "7729" = "Hroretzi",
                 "7742" = "Vertebrata",
                 "8022" = "Omykiss",
                 "8355" = "Xlaevis",
                 "8364" = "Stropicalis",
                 "9031" = "Ggallus",
                 "9606" = "Hsapiens",
                 "9913" = "Btaurus",
                 "9986" = "Ocuniculus",
                 "10090" = "Mmusculus",
                 "10116" = "Rnorvegicus",
                 "10117" = "Rrattus")

  if (code %in% names(lookup))
      return(lookup[[code]])

  NA

} # ncbiTaxonimicCodeToLinnaean
#----------------------------------------------------------------------------------------------------
if(!interactive())
  run("..")
