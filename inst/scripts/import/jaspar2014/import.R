# jaspar2014/import.R
#-----------------------------------------------------------------------------------
library (org.Hs.eg.db)
library (org.Mm.eg.db)
#------------------------------------------------------------------------------------------------------------------------
options (stringsAsFactors=FALSE)
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
run = function (dataDir)
{
  dataDir <- file.path(dataDir,"jaspar2014")
  rawMatrixList <- readRawMatrices (dataDir)
  matrices <- extractMatrices (rawMatrixList)
  tbl.anno <- createAnnotationTable(dataDir)
  tbl.md <- createMetadataTable (tbl.anno, matrices)
  matrices <- normalizeMatrices (matrices)
  matrices <- renameMatrices (matrices, tbl.md)
  serializedFile <- file.path(dataDir, "jaspar2014.RData")
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
    

  filename <- file.path(dataDir, "matrix_data.txt")
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
  
  return (matrices)
} # extractMatrices
#------------------------------------------------------------------------------------------------------------------------
createAnnotationTable <- function(dataDir)
{
    file <- file.path(dataDir, "Matrix.txt")
    tbl.matrix <-  read.table(file, sep='\t', header=FALSE, as.is=TRUE)
    colnames(tbl.matrix) <- c('id', 'category', 'mID', 'version', 'binder')
      
    file <- file.path(dataDir, "MATRIX_PROTEIN.txt")
    stopifnot(file.exists(file))
    tbl.protein <- read.table(file, sep='\t', header=FALSE, as.is=TRUE)
    colnames(tbl.protein) <- c('id', 'proteinID')

    file <- file.path(dataDir, "MATRIX_SPECIES.txt")
    stopifnot(file.exists(file))
    tbl.species <- read.table(file, sep='\t', header=FALSE, as.is=TRUE)
    colnames(tbl.species) <- c('id', 'speciesID')
    
    file <- file.path(dataDir, "MATRIX_ANNOTATION.txt")
    stopifnot(file.exists(file))
    tbl.anno <- read.table(file, sep='\t', header=FALSE, as.is=TRUE, quote="")
    colnames(tbl.anno) <- c('id', 'attribute', 'value')

    file <- file.path(dataDir, "MATRIX_ANNOTATION.txt")
    tbl.family  <- subset(tbl.anno, attribute=='family') [, -2];   
    colnames(tbl.family) <- c('id', 'family')
       # create 5 2-column data.frames out of tbl.anno
       # which can all then be merged on the "id" column
    tbl.tax <- subset(tbl.anno, attribute=='tax_group') [,-2]; 
    colnames(tbl.tax) <- c('id', 'tax')

    tbl.class <- subset(tbl.anno, attribute=='class') [,-2];     
    colnames(tbl.class) <- c('id', 'class')

    tbl.comment <- subset(tbl.anno, attribute=='comment')[,-2];    
    colnames(tbl.comment) <- c('id', 'comment')

    tbl.pubmed  <- subset(tbl.anno, attribute=='medline')[,-2];    
    colnames(tbl.pubmed) <- c('id', 'pubmed')

    tbl.type    <- subset(tbl.anno, attribute=='type')[,-2];       
    colnames(tbl.type) <- c('id', 'type')

    tbl.md <- merge(tbl.matrix, tbl.species, all.x=TRUE)
    tbl.md <- merge(tbl.md, tbl.protein, all.x=TRUE)
    tbl.md <- merge(tbl.md, tbl.family, all.x=TRUE)
    tbl.md <- merge(tbl.md, tbl.tax, all.x=TRUE)
    tbl.md <- merge(tbl.md, tbl.class, all.x=TRUE)
    tbl.md <- merge(tbl.md, tbl.pubmed, all.x=TRUE)
    tbl.md <- merge(tbl.md, tbl.type, all.x=TRUE)

    fullID <- paste(tbl.md$mID, tbl.md$version, sep='.')
    tbl.md <- cbind(fullID, tbl.md, stringsAsFactors=FALSE)

    invisible(tbl.md)
    
} # createAnnotationTable
#-------------------------------------------------------------------------------
createMetadataTable = function (tbl.anno, matrices)
{
  options (stringsAsFactors=FALSE)
  tbl.md = data.frame ()
  matrix.ids = names (matrices) 
  dataSource = 'JASPAR_2014'
  
  for (matrix.id in matrix.ids) {
    matrix = matrices [[matrix.id]]
    tbl.sub = subset (tbl.anno, fullID==substr(matrix.id,0,8))
    if (nrow (tbl.sub) > 1) {  
        # the binder is a dimer, perhaps a homodimer, and two proteinIDs are provided. Arnt::Ahr
        # some others, a sampling:  P05412;P01100, P08047, P22814;Q24040, EAW80806;EAW53510
      dimer = paste (unique (tbl.sub$proteinID), collapse=';')
      tbl.sub [1, 'proteinID'] = dimer
      }
    anno = as.list (tbl.sub [1,])
    taxon.code = anno$speciesID
    geneId.info = assignGeneId (anno$proteinID)
    new.row = list (providerName=anno$binder,
                    providerId=anno$fullID,
                    dataSource=dataSource,
                    geneSymbol=anno$binder,
                    geneId=geneId.info$geneId,
                    geneIdType=geneId.info$type,
                    proteinId=anno$proteinID,
                    proteinIdType=guessProteinIdentifierType (anno$proteinID),
                    organism=convertTaxonCode(anno$speciesID),
                    sequenceCount=as.integer (mean (colSums (matrix))),
                    bindingSequence=NA_character_,
                    bindingDomain=anno$class,
                    tfFamily=anno$family,
                    experimentType=anno$type,
                    pubmedID=anno$pubmed)
    tbl.md = rbind (tbl.md, data.frame (new.row, stringsAsFactors=FALSE))
    full.name = sprintf ('%s-%s-%s-%s', convertTaxonCode(anno$speciesID), dataSource, anno$binder, anno$fullID)
    rownames (tbl.md) [nrow (tbl.md)] = full.name
    } # for i

      # Mmusculus-JASPAR_CORE-NF-kappaB-MA0061.1 has 'NA' for proteinID, not <NA>
    NA.string.indices = grep ('NA', tbl.md$proteinId)
    if (length (NA.string.indices) > 0)
      tbl.md$proteinId [NA.string.indices] = NA
   invisible (tbl.md)
}#createMetadataTable
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
#------------------------------------------------------------------------------------------------------------------------
# full names:   ('Mus musculus', 'Rattus norvegicus', 'Rattus rattus', 'Arabidopsis thaliana', 'Pisum sativum', 
#                'Nicotiana sylvestris', 'Petunia hybrida', 'Antirrhinum majus', 'Hordeum vulgare', 'Triticum aestivam',
#                'Zea mays', 'Saccharomyces cerevisiae', 'Caenorhabditis elegans', 'Drosophila melanogaster',
#                'Halocynthia roretzi', 'Vertebrata', 'Onchorhynchus mykiss', 'Xenopus laevis', 'Xenopus tropicalis', 
#                'Gallus gallus', 'Homo sapiens', 'Bos taurus', 'Oryctolagus cuniculus'),
convertTaxonCode = function (ncbi.code)
{
  #if (is.na (ncbi.code))  
  #  return (NA_character_)
  ncbi.code = as.character (ncbi.code)
  if (ncbi.code %in% c ('-', NA_character_, 'NA'))
    return ('Vertebrata')

  tbl = data.frame (code= c('10090', '10116', '10117', '3702', '3888', '4094', '4102',
                            '4151', '4513', '4565', '4577', '4932', '6239', '7227', '7729',
                            '7742', '8022', '8355', '8364', '9031', '9606', '9913', '9986'),
                    name=c ('Mmusculus', 'Rnorvegicus', 'Rrattus', 'Athaliana', 'Psativum', 
                            'Nsylvestris', 'Phybrida', 'Amajus', 'Hvulgare', 'Taestivam',
                            'Zmays', 'Scerevisiae', 'Celegans', 'Dmelanogaster',
                            'Hroretzi', 'Vertebrata', 'Omykiss', 'Xlaevis', 'Xtropicalis', 
                            'Gallus', 'Hsapiens', 'Btaurus', 'Ocuniculus'),
                    stringsAsFactors=FALSE)

  ncbi.code = as.character (ncbi.code)
  index = which (tbl$code == ncbi.code)
  if (length (index) == 1)
    return (tbl$name [index])
  if (nchar(ncbi.code)>6)
    return ('Vertebrata')
  else {
    browser()
    write (sprintf ("unable to map organism code |%s|", ncbi.code), stderr ())
    return (NA_character_)
    }

} # convertTaxonCode
#----------------------------------------------------------------------------------------------------
assignGeneId = function (proteinId)
{
  if (!exists ('id.uniprot.human')) {

    tbl = toTable (org.Hs.egUNIPROT)
    id.uniprot.human <<- as.character (tbl$gene_id)
    names (id.uniprot.human) <<- tbl$uniprot_id

    tbl = toTable (org.Hs.egREFSEQ)
    tbl = tbl [grep ('^NP_', tbl$accession),]
    id.refseq.human <<- as.character (tbl$gene_id)
    names (id.refseq.human) <<- tbl$accession

    tbl = toTable (org.Mm.egUNIPROT)
    id.uniprot.mouse <<- as.character (tbl$gene_id)
    names (id.uniprot.mouse) <<- tbl$uniprot_id

    tbl = toTable (org.Mm.egREFSEQ)
    tbl = tbl [grep ('^NP_', tbl$accession),]
    id.refseq.mouse <<- as.character (tbl$gene_id)
    names (id.refseq.mouse) <<- tbl$accession
    }

  proteinId = strsplit (proteinId, '\\.')[[1]][1]   # remove any trailing '.*'

  if (proteinId %in% names (id.uniprot.human))
    return (list (geneId=as.character (id.uniprot.human [proteinId]), type='ENTREZ'))

  if (proteinId %in% names (id.uniprot.mouse))
    return (list (geneId=as.character (id.uniprot.mouse [proteinId]), type='ENTREZ'))

  if (proteinId %in% names (id.refseq.human))
    return (list (geneId=as.character (id.refseq.human [proteinId]), type='ENTREZ'))

  if (proteinId %in% names (id.refseq.mouse))
    return (list (geneId=as.character (id.refseq.mouse [proteinId]), type='ENTREZ'))

  found.leading.Y = length (grep ("^Y", proteinId, perl=TRUE))

  if (found.leading.Y == 1)
    return (list (geneId=proteinId, type='SGD'))

  return (list (geneId=NA_character_, type=NA_character_))

} # assignGeneId
#-----------------------------------------------------------------------------
