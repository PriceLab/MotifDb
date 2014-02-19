# jaspar/import.R
# 
#------------------------------------------------------------------------------------------------------------------------
library (org.Hs.eg.db)
library (org.Mm.eg.db)
#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
run = function (dataDir)
{
  dataDir <- file.path(dataDir, "jaspar")
  tbl.rmat = readRawMatrices (dataDir)
  matrices = convertRawMatricesToStandard (tbl.rmat)
  tbl.anno = createAnnotationTable (dataDir)
  tbl.md = createMetadataTable (tbl.anno, matrices)
  matrices = renameMatrices (matrices, tbl.md)
  matrices = normalizeMatrices (matrices)
  serializedFile <- "jaspar.RData"
  save (matrices, tbl.md, file=serializedFile)
  printf("saved %d matrices to %s", length(matrices), serializedFile)
  printf("next step: copy %s to <packageRoot>/MotifDb/inst/extdata, rebuild package", serializedFile)

} # run
#------------------------------------------------------------------------------------------------------------------------

readRawMatrices = function (dataDir)
{
  file <- file.path(dataDir, 'MATRIX_DATA.txt')
  #tbl.matrices = read.table (file, sep='\t', header=FALSE, as.is=TRUE, colClasses=c ('character', 'character', 'numeric', 'numeric'))
  #colnames (tbl.matrices) = c ('id', 'base', 'pos', 'count')
  tbl.matrices = read.table (file, sep='\n', header=FALSE )
  
  invisible (tbl.matrices)

} # readRawMatrices
#------------------------------------------------------------------------------------------------------------------------
convertRawMatricesToStandard = function (tbl.rmat)
{
  #matrix.ids = unique (tbl.rmat$id)
  #result =  vector ('list', length (matrix.ids))

  #i = 1
  #for (matrix.id in matrix.ids) {
  #  tbl.sub = subset (tbl.rmat, id == matrix.id)
  #    # sanity check this rather unusual representation of a position count matrix
  #  base.count = as.data.frame (table (tbl.sub$base))
  #  stopifnot (base.count$Var1 == c ('A', 'C', 'G', 'T'))
      # conservative length check.  actually expect sequence lengths of 6 - 20 bases
  #  if  (base.count$Freq [1] < 4 && base.count$Freq [1] > 30) {
  #    printf ('matrix.id %s has sequence of length %d', matrix.id, base.count$Freq [1])
  #    }
  #  stopifnot (all (base.count$Freq == base.count$Freq [1]))
  # nucleotide.counts = tbl.sub$count
  #  row.names = c('A', 'C', 'G', 'T'); 
  #  col.names = 1:(nrow (tbl.sub) / 4)
  #  m = matrix (nucleotide.counts, byrow=TRUE, nrow=4, dimnames=list(row.names, col.names))
  #  result [[i]] = m
  #  i = i + 1
  #  } # for matrix.id

  #names (result) = matrix.ids

  result<- vector("list",length=593)
  matrices.id <-grep('>',tbl.matrices[[1]],value=TRUE)
  names (result) = matrices.id
  #Get a list of the ID's and then assign them as names 
  m.d <-grep('>',tbl.matrices[[1]],value=TRUE,invert=TRUE)
  data.list=list()
  x = 1
  y = 1
  while (x < length(m.d))
    {   
      matrices.position <- rbind(m.d[x],m.d[x+1],m.d[x+2],m.d[x+3])
      matrices.formatted <- gsub("\t"," ",matrices.position)
      #matrices.formatted2 <- as.numeric(strsplit(matrices.position, "\t"))
      row.names = c('A', 'C', 'G', 'T')
      col.names = 1:(ncol (matrices.formatted))
      #matrices.matrix <- matrix (matrices.formatted, byrow=TRUE, nrow=4,dimnames=list(row.names,))
      result[[y]] <- matrices.formatted
      y=y+1
      x=x+4
    }
  browser()
  return (result)

} # convertRawMatricesToStandard 
#------------------------------------------------------------------------------------------------------------------------
# read 'mysql' tables provide by jaspar: 
#          MATRIX.txt:  9229	CORE	MA0001	1	AGL3
#  MATRIX_PROTEIN.txt: 9229	P29383
#  MATRIX_SPECIES.txt: 9229	3702
#  MATRIX_ANNOTATION.txt: 
#     9229	class	Other Alpha-Helix
#     9229	comment	-
#     9229	family	MADS
#     9229	medline	7632923
#     9229	tax_group	plants
#     9229	type	SELEX
createAnnotationTable = function (dataDir)
{
  file <- file.path(dataDir, "MATRIX.txt")
  tbl.matrix =  read.table (file, sep='\t', header=F, as.is=TRUE)
  colnames (tbl.matrix) = c ('id', 'category', 'mID', 'version', 'binder')

  file <- file.path(dataDir, "MATRIX_PROTEIN.txt")
  tbl.protein = read.table (file, sep='\t', header=F, as.is=TRUE)
  colnames (tbl.protein) =  c ('id', 'proteinID')

  file <- file.path(dataDir, "MATRIX_SPECIES.txt")
  tbl.species = read.table (file, sep='\t', header=F, as.is=TRUE)
  colnames (tbl.species) = c ('id', 'speciesID')

  file <- file.path(dataDir, "MATRIX_ANNOTATION.txt")
  tbl.anno = read.table (file, sep='\t', header=F, as.is=TRUE, quote="")
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
#                      names=character(),                    # source-species-gene: UniPROBE-Mmusculus-Rhox11-306b; ScerTF-Scerevisiae-ABF2-e73a
#                      nativeNames=character(),              # badis.ABF2, Cell08/Rhox11_2205.1_pwm.txt
#                      geneSymbols=character(),              # ABF2, Rhox11
#                      sequenceCounts=integer(),             # often NA
#                      organisms=character(),                # Scerevisiae, Mmusculus
#                      bindingMolecules=character(),         # YMR072W, 194738
#                      bindingMoleculeIdTypes=character(),   # SGD, entrez gene
#                      bindingDomainTypes=character(),       # NA, Homeo
#                      dataSources=character(),              # ScerTF, UniPROBE
#                      experimentTypes=character(),          # NA, protein-binding microarray
#                      pubmedIDs=character(),                # 19111667, 1858359
#                      tfFamilies=character())               # NA, NA
#
# from these
#
createMetadataTable = function (tbl.anno, matrices)
{
  options (stringsAsFactors=FALSE)
  tbl.md = data.frame ()
  matrix.ids = names (matrices) 
  dataSource = 'JASPAR_CORE'
  
  for (matrix.id in matrix.ids) {
    matrix = matrices [[matrix.id]]
    stopifnot (length (intersect (matrix.id, tbl.anno$id)) == 1)
    tbl.sub = subset (tbl.anno, id==matrix.id)
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

} # createMetadataTable
#------------------------------------------------------------------------------------------------------------------------
renameMatrices = function (matrices, tbl.md)
{
  stopifnot (length (matrices) == nrow (tbl.md))
  names (matrices) = rownames (tbl.md)
  invisible (matrices)

} # renameMatrices
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
