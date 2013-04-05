# stamlab/import.R
#------------------------------------------------------------------------------------------------------------------------
library (org.Hs.eg.db)
library (org.Mm.eg.db)
#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
run = function (dataDir)
{
  dataDir <- file.path(dataDir, "stamlab")
  rawMatrixList <- readRawMatrices (dataDir)
  novels <- readNovelStatus (dataDir)
  matrices <- extractAndNormalizeMatrices (rawMatrixList)
  tbl.md <- createMetadataTable (matrices, novels)
  matrices <- renameMatrices (matrices, tbl.md)

  serializedFile <- "stamlab.RData"
  save (matrices, tbl.md, file=serializedFile)
  printf("saved %d matrices to %s", length(matrices), serializedFile)
  printf("next step:  copy %s to <packageRoot>/MotifDb/inst/extdata, rebuild package", serializedFile)

} # run
#------------------------------------------------------------------------------------------------------------------------
readRawMatrices = function (dataDir)
{
  filename <- file.path(dataDir, "de.novo.pwm")
  all.lines = scan (filename, what=character(0), sep='\n', quiet=TRUE)
  title.lines = grep ('UW.Motif', all.lines)
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
readNovelStatus = function (dataDir)
{
  filename <- file.path(dataDir, "novels.txt")
  novels = scan (filename, what=character(0), sep='\n', quiet=TRUE)
  all.names = sprintf ('UW.Motif.%04d', 1:683)
  status = rep (FALSE, length (all.names))
  true.novels = match (novels, all.names)
  status [true.novels] = TRUE
  names (status) = all.names
  return (status)
  
} # readNovelStatus
#------------------------------------------------------------------------------------------------------------------------
extractAndNormalizeMatrices = function (pwm.list)
{
  ms = sapply (pwm.list, function (element) element$matrix)
  nms = normalizeMatrices (ms)
  names (nms) = sapply (pwm.list, function (element) element$title)
  return (nms)

} # extractAndNormalizeMatrices
#------------------------------------------------------------------------------------------------------------------------
#  matrices = sapply (list.pwms, function (pwm) pwm$matrix)
#  matrix.names = sapply (list.pwms, function (pwm) pwm$title)
#  names (matrices) = matrix.names
convertRawMatricesToStandard = function (tbl.rmat)
{
  matrix.ids = unique (tbl.rmat$id)
  result =  vector ('list', length (matrix.ids))

  i = 1
  for (matrix.id in matrix.ids) {
    tbl.sub = subset (tbl.rmat, id == matrix.id)
      # sanity check this rather unusual representation of a position count matrix
    base.count = as.data.frame (table (tbl.sub$base))
    stopifnot (base.count$Var1 == c ('A', 'C', 'G', 'T'))
      # conservative length check.  actually expect sequence lengths of 6 - 20 bases
    if  (base.count$Freq [1] < 4 && base.count$Freq [1] > 30) {
      printf ('matrix.id %s has sequence of length %d', matrix.id, base.count$Freq [1])
      }
    stopifnot (all (base.count$Freq == base.count$Freq [1]))
    nucleotide.counts = tbl.sub$count
    row.names = c('A', 'C', 'G', 'T'); 
    col.names = 1:(nrow (tbl.sub) / 4)
    m = matrix (nucleotide.counts, byrow=TRUE, nrow=4, dimnames=list(row.names, col.names))
    result [[i]] = m
    i = i + 1
    } # for matrix.id

  names (result) = matrix.ids
  return (result)

} # convertRawMatricesToStandard 
#------------------------------------------------------------------------------------------------------------------------
createAnnotationTable = function ()
{
  tbl.matrix =  read.table ('MATRIX.txt', sep='\t', header=F, as.is=TRUE)
  colnames (tbl.matrix) = c ('id', 'category', 'mID', 'version', 'binder')

  tbl.protein = read.table ('MATRIX_PROTEIN.txt', sep='\t', header=F, as.is=TRUE)
  colnames (tbl.protein) =  c ('id', 'proteinID')

  tbl.species = read.table ('MATRIX_SPECIES.txt', sep='\t', header=F, as.is=TRUE)
  colnames (tbl.species) = c ('id', 'speciesID')

  tbl.anno = read.table ('MATRIX_ANNOTATION.txt', sep='\t', header=F, as.is=TRUE, quote="")
  colnames (tbl.anno) = c ('id', 'attribute', 'value')

  tbl.family  = subset (tbl.anno, attribute=='family') [, -2];   
  colnames (tbl.family) = c ('id', 'family')

  tbl.tax     = subset (tbl.anno, attribute=='tax_group') [,-2]; 
  colnames (tbl.tax) = c ('id', 'tax')

  tbl.class   = subset (tbl.anno, attribute=='class') [,-2];     
  colnames (tbl.class) = c ('id', 'class')

  tbl.comment = subset (tbl.anno, attribute=='comment')[,-2];    
  colnames (tbl.comment) = c ('id', 'comment')

  tbl.pubmed  = subset (tbl.anno, attribute=='medline')[,-2];    
  colnames (tbl.pubmed) = c ('id', 'pubmed')

  tbl.type    = subset (tbl.anno, attribute=='type')[,-2];       
  colnames (tbl.type) = c ('id', 'type')


  tbl.md = merge (tbl.matrix, tbl.species, all.x=TRUE)
  tbl.md = merge (tbl.md, tbl.protein, all.x=TRUE)
  tbl.md = merge (tbl.md, tbl.family, all.x=TRUE)
  tbl.md = merge (tbl.md, tbl.tax, all.x=TRUE)
  tbl.md = merge (tbl.md, tbl.class, all.x=TRUE)
  tbl.md = merge (tbl.md, tbl.pubmed, all.x=TRUE)
  tbl.md = merge (tbl.md, tbl.type, all.x=TRUE)

  fullID = paste (tbl.md$mID, tbl.md$version, sep='.')
  tbl.md = cbind (fullID, tbl.md, stringsAsFactors=FALSE)

  invisible (tbl.md)

} # createAnnotationTable
#------------------------------------------------------------------------------------------------------------------------
# assemble these columns:
#                      names=character(),                    # species-source-gene: stamlab-Hsapiens-UW.Motif.0001
#                      nativeNames=character(),              # UW.Motif.0001
#                      geneSymbols=character(),              # NA
#                      sequenceCounts=integer(),             # NA
#                      organisms=character(),                # Hsapiens
#                      bindingMolecules=character(),         # NA
#                      bindingMoleculeIdTypes=character(),   # NA
#                      bindingDomainTypes=character(),       # NA
#                      dataSources=character(),              # stamlab
#                      experimentTypes=character(),          # digital genomic footprinting
#                      pubmedIDs=character(),                # 22959076
#                      tfFamilies=character())               # NA
#
# from these
#
createMetadataTable = function (matrices, novels)
{
  options (stringsAsFactors=FALSE)
  tbl.md = data.frame ()
  matrix.ids = names (matrices) 
  
  for (matrix.id in matrix.ids) {
    matrix = matrices [[matrix.id]]
    taxon.code = 'Hsapiens'
    geneId.info = NA
    new.row = list (providerName=matrix.id,
                    providerId=matrix.id,
                    dataSource='stamlab',
                    geneSymbol=NA,
                    geneId=NA,
                    geneIdType=NA,
                    proteinId=NA,
                    proteinIdType=NA,
                    organism='Hsapiens',
                    sequenceCount=NA,
                    bindingSequence=NA_character_,
                    bindingDomain=NA,
                    tfFamily=NA,
                    experimentType='digital genomic footprinting',
                    pubmedID='22959076')
    tbl.md = rbind (tbl.md, data.frame (new.row, stringsAsFactors=FALSE))
    full.name = sprintf ('%s-%s-%s', 'Hsapiens', 'stamlab', matrix.id)
    rownames (tbl.md) [nrow (tbl.md)] = full.name
    } # for i

  novelPFM = rep ('knownMotif', nrow (tbl.md))
  novels.ordered = novels [tbl.md$providerName]  # make sure we follow the order in the tbl
  novelPFM [which (novels.ordered)] = 'novelMotif'
  tbl.md$geneId = novelPFM
  tbl.md$geneIdType = rep ('comment', nrow (tbl.md))

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
  else {
    write (sprintf (" unable to map organism code |%s|", ncbi.code), stderr ())
    return (NA_character_)
    }

} # convertTaxonCode
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

  if (length (grep ('^EAW', moleculeName)) == 1)
    return ('NCBI')

  if (length (grep ('^EAX', moleculeName)) == 1)
    return ('NCBI')

  if (length (grep ('^NP_', moleculeName)) == 1)
    return ('REFSEQ')

  if (length (grep ('^BAD', moleculeName)) == 1)
    return ('EMBL')

   return (NA_character_)

} # guessProteinIdentifierType
#------------------------------------------------------------------------------------------------------------------------
normalizeMatrices = function (matrices)
{
  mtx.normalized = sapply (matrices, function (mtx) apply (mtx, 2, function (colvector) colvector / sum (colvector)))
  invisible (mtx.normalized)

} # normalizeMatrices
#------------------------------------------------------------------------------------------------------------------------
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
#------------------------------------------------------------------------------------------------------------------------
parsePwm = function (text)
{
  #printf ('parsing pwm %s', text [1])
  lines = strsplit (text, '\t')
  title = lines [[1]][1]
  consensus.sequence = lines [[1]][2]
  line.count = length (lines)
  #printf ('%s: %s', title, consensus.sequence)

  result = matrix (nrow=line.count-1, ncol=4, dimnames=list(1:(line.count-1), c ('A','C','G','T')))  
  row = 1
  for (line in lines [2:line.count]) {
    result [row,] = as.numeric (line)
    row = row + 1
    } # for line

  result = t (result)
    
  return (list (title=title, consensus.sequence=consensus.sequence, matrix=result))

} # parsePwm
#------------------------------------------------------------------------------------------------------------------------
