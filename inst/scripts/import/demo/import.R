# MotifDb/inst/scripes/import/demo/import.R
#------------------------------------------------------------------------------------------------------------------------
options (stringsAsFactors=FALSE)
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
run = function (dataDir)
{
  dataDir <- file.path(dataDir, "demo")
  rawMatrixList <- readRawMatrices (dataDir)
  matrices <- extractMatrices (rawMatrixList)
  tbl.md <- createMetadataTable (dataDir, matrices,
                                 raw.metadata.filename="md-raw.tsv")
  matrices <- normalizeMatrices (matrices)
  matrices <- renameMatrices (matrices, tbl.md)

  serializedFile <- file.path(dataDir, "demo.RData")
  printf("writing %s to %s", "demo.RData", dataDir)

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
    

  filename <- file.path(dataDir, "jaspar2014", "sample.pcm")
  printf("checking for readable matrix file:")
  printf("     %s", filename)
  stopifnot(file.exists(filename))
  
  all.lines = scan (filename, what=character(0), sep='\n', quiet=TRUE)
  title.lines = grep ('^>', all.lines)
  title.line.count <<- length (title.lines)
  max = title.line.count - 1

  pwms = list ()
  
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
  matrix.names <- sub("^>", "", matrix.names)
  names (matrices) <- matrix.names
  
  matrices

} # extractMatrices
#------------------------------------------------------------------------------------------------------------------------
createMetadataTable = function (dataDir, matrices, raw.metadata.filename)
{
  filename <- file.path(dataDir, "jaspar2014", "md-raw.tsv")
  printf("checking for readable metadata file:")
  printf("   %s", filename)
  stopifnot(file.exists(filename))

  tbl.raw <- read.table(filename, sep="\t", header=TRUE, as.is=TRUE)
  tbl.md = data.frame ()
  matrix.ids = names (matrices)
  
  for (matrix.id in matrix.ids) {
    matrix <- matrices[[matrix.id]]
    short.matrix.name <- sub("\\..*$", "", matrix.id)
    stopifnot(length(grep(short.matrix.name, tbl.raw$ma.name)) == 1)
    md <- as.list(subset(tbl.raw, ma.name==short.matrix.name))
    dataSource <- "jaspar2014"
    #browser()
    organism <- ncbiTaxonimicCodeToLinnaean(md$ncbi.tax.code)
    new.row = list (providerName=matrix.id,
                    providerId=matrix.id,
                    dataSource=dataSource,
                    geneSymbol=md$gene.symbol,
                    geneId=NA,
                    geneIdType=NA,
                    proteinId=md$uniprot,
                    proteinIdType=guessProteinIdentifierType(md$uniprot),
                    organism=ncbiTaxonimicCodeToLinnaean(md$ncbi.tax.code),
                    sequenceCount=max(colSums(matrix)),
                    bindingSequence=NA_character_,
                    bindingDomain=NA,
                    tfFamily=md$family,
                    experimentType=md$type,
                    pubmedID="24194598")
    tbl.md = rbind (tbl.md, data.frame (new.row, stringsAsFactors=FALSE))
    full.name = sprintf ('%s-%s-%s', organism, dataSource, matrix.id)
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
  if (nchar (moleculeName) == 0)
    return (NA_character_)
  if (is.na (moleculeName))
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
   #browser()
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
ncbiTaxonimicCodeToLinnaean <- function(code)
{
  code <- as.character(code)
  
  lookup <- list("3702" = "Arabidopsis thaliana",
                 "3888" = "Pisum sativum",
                 "4094" = "Nicontia sp.",
                 "4102" = "Petunia hybrida",
                 "4151" = "Antirrhinum majus",
                 "4513" = "Hordeum vulgare",
                 "4565" = "Triticum aestivam",
                 "4577" = "Zea mays",
                 "4932" = "Saccaromyces cerevisiae",
                 "6239" = "Caenorhabditis elegans",
                 "7227" = "Drosophila melanogaster",
                 "7729" = "Halocynthia roretzi",
                 "7742" = "Vertebrata",
                 "8022" = "Onchorhynchus mykiss",
                 "8355" = "Xenopus laevis",
                 "8364" = "Silurana tropicalis",
                 "9031" = "Gallus gallus",
                 "9606" = "Homo sapiens",
                 "9913" = "Bos taurus",
                 "9986" = "Oryctolagus cuniculus",
                 "10090" = "Mus musculus",
                 "10116" = "Rattus norvegicus",
                 "10117" = "Rattus rattus")

  if (code %in% names(lookup))
      return(lookup[[code]])

  NA

} # ncbiTaxonimicCodeToLinnaean
#----------------------------------------------------------------------------------------------------
