# jolma2013/import.R
#------------------------------------------------------------------------------------------------------------------------
library (org.Hs.eg.db)
library (org.Mm.eg.db)
#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
#run = function (dataDir)
#{
#  dataDir <- file.path(dataDir, "jaspar")
#  tbl.rmat = readRawMatrices (dataDir)
#  matrices = convertRawMatricesToStandard (tbl.rmat)
#  tbl.anno = createAnnotationTable (dataDir)
#  tbl.md = createMetadataTable (tbl.anno, matrices)
#  matrices = renameMatrices (matrices, tbl.md)
#  matrices = normalizeMatrices (matrices)
#  serializedFile <- "jaspar.RData"
#  save (matrices, tbl.md, file=serializedFile)
#  printf("saved %d matrices to %s", length(matrices), serializedFile)
#  printf("next step: copy %s to <packageRoot>/MotifDb/inst/extdata, rebuild package", serializedFile)
#
#} # run
#-------------------------------------------------------------------------------
run = function (dataDir)
{
   dataDir <- file.path(dataDir, "jolma2013")
   x <- readRawMatrices(dataDir)
   matrices <- x$matrices
   tbl.mAnno <- x$anno
   tbl.mAnno$symbol[tbl.mAnno$symbol=="Tp53"] <- "Trp53"
   tbl.mAnno$symbol[tbl.mAnno$symbol=="Tp73"] <- "Trp73"

   tbl.anno <- readRawAnnotationTable(dataDir)
   mouse.syms <- subset(tbl.anno, organism=="mouse")$symbol
   organism <- rep("Hsapiens", nrow(tbl.mAnno))
   mouse.indices <- which(names(matrices) %in% mouse.syms)
   organism[mouse.indices] <- "Mmusculus"
   tbl.mAnno <- cbind(tbl.mAnno, organism, stringsAsFactors=FALSE)

   tbl.md <- createMetadataTable(tbl.mAnno, matrices)
   rownames(tbl.md) <- tbl.md$providerName

   matrices <- normalizeMatrices(matrices)
   names(matrices) <- tbl.md$providerName

   serializedFile <- "jolma2013.RData"
   save (matrices, tbl.md, file=serializedFile)
   printf("saved %d matrices to %s", length(matrices), serializedFile)
   printf("next step: copy %s to <packageRoot>/MotifDb/inst/extdata, and rebuild package", serializedFile)


   browser("rr")

} # run
#-------------------------------------------------------------------------------
normalizeMatrices = function (matrices)
{
  mtx.normalized = sapply (matrices, function (mtx) apply (mtx, 2, function (colvector) colvector / sum (colvector)))
  invisible (mtx.normalized)

} # normalizeMatrices
#------------------------------------------------------------------------------------------------------------------------
createMetadataTable = function (tbl.mAnno, matrices)
{
    options (stringsAsFactors=FALSE)
    tbl.md = data.frame ()
    dataSource = 'Jolma2013'

    provider.names <- createProviderNames(tbl.mAnno)

    stopifnot(length(which(duplicated(provider.names))) == 0)

    stopifnot(nrow(tbl.mAnno) == length(matrices))

    size <- length(provider.names)
    tbl <- data.frame(list(providerName=vector("character", size),
                           providerId=vector("character", size),
                           dataSource=vector("character",size),
                           geneSymbol=vector("character",size),
                           geneId=vector("character",size),
                           geneIdType=vector("character",size),
                           proteinId=vector("character",size),
                           proteinIdType=vector("character",size),
                           organism=vector("character",size),
                           sequenceCount=vector("integer", size),
                           bindingSequence=vector("character",size),
                           bindingDomain=vector("character",size),
                           tfFamily=vector("character",size),
                           experimentType=vector("character",size),
                           pubmedID=vector("character",size)))

    #rownames(tbl) <- provider.names
                           
    geneIDs <- determineGeneIDs(tbl.mAnno$symbol)
    counts = as.integer (sapply (matrices, function (mtx) as.integer (mean (round (colSums (mtx))))))

    for (r in 1:length(provider.names)){
        tbl[r, "providerName"] <- provider.names[r]
        tbl[r, "providerId"]   <- tbl.mAnno$symbol[r]
        tbl[r, "dataSource"]   <- "jolma2013"
        tbl[r, "geneSymbol"]   <- tbl.mAnno$symbol[r]
        tbl[r, "geneId"]   <- geneIDs[r]
        tbl[r, "geneIdType"]   <- "ENTREZ"
        tbl[r, "proteinId"]   <- NA_character_
        tbl[r, "proteinIdType"]   <- NA_character_
        tbl[r, "organism"]   <-  tbl.mAnno$organism[r]
        tbl[r, "sequenceCount"]   <- counts[r]
        tbl[r, "bindingSequence"]   <- tbl.mAnno$sequence[r]
        tbl[r, "bindingDomain"]   <- NA_character_
        tbl[r, "tfFamily"]   <-  tbl.mAnno$family[r]
        tbl[r, "experimentType"]   <- "SELEX"
        tbl[r, "pubmedID"]   <- "23332764"
        } # for r

   invisible (tbl)

} # createMetadataTable
##------------------------------------------------------------------------------------------------------------------------
readRawMatrices = function (dataDir)
{
    data.file <- file.path(dataDir, "TableS3-PWM-models.tsv")
    tbl.raw <- read.table(data.file, skip=15, header=FALSE, as.is=TRUE, fill=TRUE, sep="\t")

     # the tsv file, exported from the excel spreadsheet, \r converted to \n,
     # has an unconventional structure, with 5 rows per pwm, including 1 header line
     # followed by by 4 rows of base counts, in which the first column is A,C,T or G
     # each header line has 10 columns, apparently these

     # symbol family clone         type ligand          sequence batch Seed multinomial cycle
     #  BCL6B   C2H2   DBD  TGCGGG20NGA     AC TGCTTTCTAGGAATTMM     2    4   monomeric   yes
     #   CTCF   C2H2  full TAGCGA20NGCT     AJ NGCGCCMYCTAGYGGTN     2    4   monomeric   yes
     #   EGR1   C2H2   DBD  TCTCTT20NGA      Y    NMCGCCCMCGCANN     2    2   monomeric   yes
     #   EGR1   C2H2  full TACTAT20NATC     AA    NACGCCCACGCANN     2    4   monomeric    no
     #   EGR2   C2H2   DBD  TCGGCC20NGA      W       MCGCCCACGCA     2    3   monomeric    no

     # though there are some left-over column titles.  don't know quite what to make of them
     #         site:
     #         type:
     #  comment:
     #  Matrix is one of the representative PWMs:
  
    
    matrix.start.rows <- seq(3, nrow(tbl.raw), 5)
    tbl.anno <- tbl.raw[matrix.start.rows, 1:10]
    colnames(tbl.anno) <- c("symbol", "family", "clone", "type", "ligand",
                            "sequence", "batch", "Seed", "multinomial",
                            "cycle")
    tbl.indices <- data.frame(start = matrix.start.rows + 1, end = matrix.start.rows + 4)
    matrix.list = apply(tbl.indices, 1, function(x) tbl.raw[x[1]:x[2],])
    names(matrix.list) = tbl.anno$symbol

    matrix.list = lapply(matrix.list, function(mtx) {
       tmp = mtx[,-1]   # drop the ACTG column
          # identify empty or NA columns
       columns.to.discard = apply(tmp, 2, function(x) all(is.na(x) | all(x == "")))
       tmp = tmp[, ! columns.to.discard]
       tmp2 = as.numeric(unlist(tmp))
       tmp = matrix(tmp2, nrow = 4, ncol = length(tmp2)/4)
       rownames(tmp) = mtx[,1]
       colnames(tmp) = 1:ncol(tmp)
       tmp
       })

    unique.ids <- apply(tbl.anno, 1, function(s) paste(s, collapse="."))
    #names(matrix.list) <- unique.ids

       # insert some naming repairs discovered by diego diez
    
    tbl.anno$symbol[tbl.anno$symbol=="Tp53"] <- "Trp53"
    tbl.anno$symbol[tbl.anno$symbol=="Tp73"] <- "Trp73"

    invisible (list(matrices=matrix.list, anno=tbl.anno))

} # readRawMatrices
#------------------------------------------------------------------------------------------------------------------------
readRawAnnotationTable <- function(dataDir)
{
    data.file <- file.path(dataDir, "Cell_2013_Jolma_annotations.tsv")
    tbl.raw <- read.table(data.file, header=TRUE, as.is=TRUE, fill=TRUE, sep="\t")
        # no need for the protein sequence
    tbl.raw <- tbl.raw[,-3]

       # insert some naming repairs discovered by diego diez
    
    tbl.raw$symbol[tbl.raw$symbol=="Egr1_E410D_FARSDERtoFARSDDR"] <- "Egr1"

    invisible(tbl.raw)
    
} # readRawAnnotationTable
#------------------------------------------------------------------------------------------------------------------------
convertRawMatricesToStandard = function (matrices)
{
   matrix.ids = names(matrices)
   browser("cRMTS")

#  result =  vector ('list', length (matrix.ids))
#
#  i = 1
#  for (matrix.id in matrix.ids) {
#    tbl.sub = subset (tbl.rmat, id == matrix.id)
#      # sanity check this rather unusual representation of a position count matrix
#    base.count = as.data.frame (table (tbl.sub$base))
#    stopifnot (base.count$Var1 == c ('A', 'C', 'G', 'T'))
#      # conservative length check.  actually expect sequence lengths of 6 - 20 bases
#    if  (base.count$Freq [1] < 4 && base.count$Freq [1] > 30) {
#      printf ('matrix.id %s has sequence of length %d', matrix.id, base.count$Freq [1])
#      }
#    stopifnot (all (base.count$Freq == base.count$Freq [1]))
#    nucleotide.counts = tbl.sub$count
#    row.names = c('A', 'C', 'G', 'T'); 
#    col.names = 1:(nrow (tbl.sub) / 4)
#    m = matrix (nucleotide.counts, byrow=TRUE, nrow=4, dimnames=list(row.names, col.names))
#    result [[i]] = m
#    i = i + 1
#    } # for matrix.id
#
#  names (result) = matrix.ids
#  return (result)
#
} # convertRawMatricesToStandard 
##------------------------------------------------------------------------------------------------------------------------
#     anno.file <- file.path(dataDir, "Cell_2013_Jolma_annotations.tsv")
## read 'mysql' tables provide by jaspar: 
##          MATRIX.txt:  9229	CORE	MA0001	1	AGL3
##  MATRIX_PROTEIN.txt: 9229	P29383
##  MATRIX_SPECIES.txt: 9229	3702
##  MATRIX_ANNOTATION.txt: 
##     9229	class	Other Alpha-Helix
##     9229	comment	-
##     9229	family	MADS
##     9229	medline	7632923
##     9229	tax_group	plants
##     9229	type	SELEX
#createAnnotationTable = function (dataDir)
#{
#  file <- file.path(dataDir, "MATRIX.txt")
#  tbl.matrix =  read.table (file, sep='\t', header=F, as.is=TRUE)
#  colnames (tbl.matrix) = c ('id', 'category', 'mID', 'version', 'binder')
#
#  file <- file.path(dataDir, "MATRIX_PROTEIN.txt")
#  tbl.protein = read.table (file, sep='\t', header=F, as.is=TRUE)
#  colnames (tbl.protein) =  c ('id', 'proteinID')
#
#  file <- file.path(dataDir, "MATRIX_SPECIES.txt")
#  tbl.species = read.table (file, sep='\t', header=F, as.is=TRUE)
#  colnames (tbl.species) = c ('id', 'speciesID')
#
#  file <- file.path(dataDir, "MATRIX_ANNOTATION.txt")
#  tbl.anno = read.table (file, sep='\t', header=F, as.is=TRUE, quote="")
#  colnames (tbl.anno) = c ('id', 'attribute', 'value')
#
#  tbl.family  = subset (tbl.anno, attribute=='family') [, -2];   
#  colnames (tbl.family) = c ('id', 'family')
#
#  tbl.tax     = subset (tbl.anno, attribute=='tax_group') [,-2]; 
#  colnames (tbl.tax) = c ('id', 'tax')
#
#  tbl.class   = subset (tbl.anno, attribute=='class') [,-2];     
#  colnames (tbl.class) = c ('id', 'class')
#
#  tbl.comment = subset (tbl.anno, attribute=='comment')[,-2];    
#  colnames (tbl.comment) = c ('id', 'comment')
#
#  tbl.pubmed  = subset (tbl.anno, attribute=='medline')[,-2];    
#  colnames (tbl.pubmed) = c ('id', 'pubmed')
#
#  tbl.type    = subset (tbl.anno, attribute=='type')[,-2];       
#  colnames (tbl.type) = c ('id', 'type')
#
#
#  tbl.md = merge (tbl.matrix, tbl.species, all.x=TRUE)
#  tbl.md = merge (tbl.md, tbl.protein, all.x=TRUE)
#  tbl.md = merge (tbl.md, tbl.family, all.x=TRUE)
#  tbl.md = merge (tbl.md, tbl.tax, all.x=TRUE)
#  tbl.md = merge (tbl.md, tbl.class, all.x=TRUE)
#  tbl.md = merge (tbl.md, tbl.pubmed, all.x=TRUE)
#  tbl.md = merge (tbl.md, tbl.type, all.x=TRUE)
#
#  fullID = paste (tbl.md$mID, tbl.md$version, sep='.')
#  tbl.md = cbind (fullID, tbl.md, stringsAsFactors=FALSE)
#
#  invisible (tbl.md)
#
#} # createAnnotationTable
##------------------------------------------------------------------------------------------------------------------------
## assemble these columns:
##                      names=character(),                    # source-species-gene: UniPROBE-Mmusculus-Rhox11-306b; ScerTF-Scerevisiae-ABF2-e73a
##                      nativeNames=character(),              # badis.ABF2, Cell08/Rhox11_2205.1_pwm.txt
##                      geneSymbols=character(),              # ABF2, Rhox11
##                      sequenceCounts=integer(),             # often NA
##                      organisms=character(),                # Scerevisiae, Mmusculus
##                      bindingMolecules=character(),         # YMR072W, 194738
##                      bindingMoleculeIdTypes=character(),   # SGD, entrez gene
##                      bindingDomainTypes=character(),       # NA, Homeo
##                      dataSources=character(),              # ScerTF, UniPROBE
##                      experimentTypes=character(),          # NA, protein-binding microarray
##                      pubmedIDs=character(),                # 19111667, 1858359
##                      tfFamilies=character())               # NA, NA
##
## from these
##
#createMetadataTable = function (tbl.anno, matrices)
#{
#  options (stringsAsFactors=FALSE)
#  tbl.md = data.frame ()
#  matrix.ids = names (matrices) 
#  dataSource = 'JASPAR_CORE'
#  
#  for (matrix.id in matrix.ids) {
#    matrix = matrices [[matrix.id]]
#    stopifnot (length (intersect (matrix.id, tbl.anno$id)) == 1)
#    tbl.sub = subset (tbl.anno, id==matrix.id)
#    if (nrow (tbl.sub) > 1) {  
#        # the binder is a dimer, perhaps a homodimer, and two proteinIDs are provided. Arnt::Ahr
#        # some others, a sampling:  P05412;P01100, P08047, P22814;Q24040, EAW80806;EAW53510
#      dimer = paste (unique (tbl.sub$proteinID), collapse=';')
#      tbl.sub [1, 'proteinID'] = dimer
#      }
#    anno = as.list (tbl.sub [1,])
#    taxon.code = anno$speciesID
#    geneId.info = assignGeneId (anno$proteinID)
#    new.row = list (providerName=anno$binder,
#                    providerId=anno$fullID,
#                    dataSource=dataSource,
#                    geneSymbol=anno$binder,
#                    geneId=geneId.info$geneId,
#                    geneIdType=geneId.info$type,
#                    proteinId=anno$proteinID,
#                    proteinIdType=guessProteinIdentifierType (anno$proteinID),
#                    organism=convertTaxonCode(anno$speciesID),
#                    sequenceCount=as.integer (mean (colSums (matrix))),
#                    bindingSequence=NA_character_,
#                    bindingDomain=anno$class,
#                    tfFamily=anno$family,
#                    experimentType=anno$type,
#                    pubmedID=anno$pubmed)
#    tbl.md = rbind (tbl.md, data.frame (new.row, stringsAsFactors=FALSE))
#    full.name = sprintf ('%s-%s-%s-%s', convertTaxonCode(anno$speciesID), dataSource, anno$binder, anno$fullID)
#    rownames (tbl.md) [nrow (tbl.md)] = full.name
#    } # for i
#
#      # Mmusculus-JASPAR_CORE-NF-kappaB-MA0061.1 has 'NA' for proteinID, not <NA>
#    NA.string.indices = grep ('NA', tbl.md$proteinId)
#    if (length (NA.string.indices) > 0)
#      tbl.md$proteinId [NA.string.indices] = NA
#     
#   invisible (tbl.md)
#
#} # createMetadataTable
##------------------------------------------------------------------------------------------------------------------------
#renameMatrices = function (matrices, tbl.md)
#{
#  stopifnot (length (matrices) == nrow (tbl.md))
#  names (matrices) = rownames (tbl.md)
#  invisible (matrices)
#
#} # renameMatrices
##------------------------------------------------------------------------------------------------------------------------
## full names:   ('Mus musculus', 'Rattus norvegicus', 'Rattus rattus', 'Arabidopsis thaliana', 'Pisum sativum', 
##                'Nicotiana sylvestris', 'Petunia hybrida', 'Antirrhinum majus', 'Hordeum vulgare', 'Triticum aestivam',
##                'Zea mays', 'Saccharomyces cerevisiae', 'Caenorhabditis elegans', 'Drosophila melanogaster',
##                'Halocynthia roretzi', 'Vertebrata', 'Onchorhynchus mykiss', 'Xenopus laevis', 'Xenopus tropicalis', 
##                'Gallus gallus', 'Homo sapiens', 'Bos taurus', 'Oryctolagus cuniculus'),
#convertTaxonCode = function (ncbi.code)
#{
#  #if (is.na (ncbi.code))  
#  #  return (NA_character_)
#  ncbi.code = as.character (ncbi.code)
#  if (ncbi.code %in% c ('-', NA_character_, 'NA'))
#    return ('Vertebrata')
#
#  tbl = data.frame (code= c('10090', '10116', '10117', '3702', '3888', '4094', '4102',
#                            '4151', '4513', '4565', '4577', '4932', '6239', '7227', '7729',
#                            '7742', '8022', '8355', '8364', '9031', '9606', '9913', '9986'),
#                    name=c ('Mmusculus', 'Rnorvegicus', 'Rrattus', 'Athaliana', 'Psativum', 
#                            'Nsylvestris', 'Phybrida', 'Amajus', 'Hvulgare', 'Taestivam',
#                            'Zmays', 'Scerevisiae', 'Celegans', 'Dmelanogaster',
#                            'Hroretzi', 'Vertebrata', 'Omykiss', 'Xlaevis', 'Xtropicalis', 
#                            'Gallus', 'Hsapiens', 'Btaurus', 'Ocuniculus'),
#                    stringsAsFactors=FALSE)
#
#  ncbi.code = as.character (ncbi.code)
#  index = which (tbl$code == ncbi.code)
#  if (length (index) == 1)
#    return (tbl$name [index])
#  else {
#    write (sprintf (" unable to map organism code |%s|", ncbi.code), stderr ())
#    return (NA_character_)
#    }
#
#} # convertTaxonCode
##------------------------------------------------------------------------------------------------------------------------
## an empirical and not altogether trustworthy solution to identifying identifier types.
#guessProteinIdentifierType = function (moleculeName)
#{
#  if (nchar (moleculeName) == 0)
#    return (NA_character_)
#  if (is.na (moleculeName))
#    return (NA_character_) 
#
#  first.char = substr (moleculeName, 1, 1)
#
#  if (first.char == 'Y')
#    return ('SGD')
#
#  if (first.char %in% c ('P', 'Q', 'O', 'A', 'C'))
#    return ('UNIPROT')
#
#  if (length (grep ('^EAW', moleculeName)) == 1)
#    return ('NCBI')
#
#  if (length (grep ('^EAX', moleculeName)) == 1)
#    return ('NCBI')
#
#  if (length (grep ('^NP_', moleculeName)) == 1)
#    return ('REFSEQ')
#
#  if (length (grep ('^BAD', moleculeName)) == 1)
#    return ('EMBL')
#
#   return (NA_character_)
#
#} # guessProteinIdentifierType
##------------------------------------------------------------------------------------------------------------------------
#normalizeMatrices = function (matrices)
#{
#  mtx.normalized = sapply (matrices, function (mtx) apply (mtx, 2, function (colvector) colvector / sum (colvector)))
#  invisible (mtx.normalized)
#
#} # normalizeMatrices
##------------------------------------------------------------------------------------------------------------------------
#assignGeneId = function (proteinId)
#{
#  if (!exists ('id.uniprot.human')) {
#
#    tbl = toTable (org.Hs.egUNIPROT)
#    id.uniprot.human <<- as.character (tbl$gene_id)
#    names (id.uniprot.human) <<- tbl$uniprot_id
#
#    tbl = toTable (org.Hs.egREFSEQ)
#    tbl = tbl [grep ('^NP_', tbl$accession),]
#    id.refseq.human <<- as.character (tbl$gene_id)
#    names (id.refseq.human) <<- tbl$accession
#
#    tbl = toTable (org.Mm.egUNIPROT)
#    id.uniprot.mouse <<- as.character (tbl$gene_id)
#    names (id.uniprot.mouse) <<- tbl$uniprot_id
#
#    tbl = toTable (org.Mm.egREFSEQ)
#    tbl = tbl [grep ('^NP_', tbl$accession),]
#    id.refseq.mouse <<- as.character (tbl$gene_id)
#    names (id.refseq.mouse) <<- tbl$accession
#    }
#
#  proteinId = strsplit (proteinId, '\\.')[[1]][1]   # remove any trailing '.*'
#
#  if (proteinId %in% names (id.uniprot.human))
#    return (list (geneId=as.character (id.uniprot.human [proteinId]), type='ENTREZ'))
#
#  if (proteinId %in% names (id.uniprot.mouse))
#    return (list (geneId=as.character (id.uniprot.mouse [proteinId]), type='ENTREZ'))
#
#  if (proteinId %in% names (id.refseq.human))
#    return (list (geneId=as.character (id.refseq.human [proteinId]), type='ENTREZ'))
#
#  if (proteinId %in% names (id.refseq.mouse))
#    return (list (geneId=as.character (id.refseq.mouse [proteinId]), type='ENTREZ'))
#
#  found.leading.Y = length (grep ("^Y", proteinId, perl=TRUE))
#
#  if (found.leading.Y == 1)
#    return (list (geneId=proteinId, type='SGD'))
#
#  return (list (geneId=NA_character_, type=NA_character_))
#
#} # assignGeneId
#------------------------------------------------------------------------------- 
determineGeneIDs <- function(symbols)
{
    geneIDs <- mget(symbols, org.Hs.egSYMBOL2EG, ifnotfound=NA)
    failures <- unique(names(which(is.na(geneIDs))))
    geneIDs <- geneIDs[which(!is.na(geneIDs))]

    alias.failures <- c()
    
    if(length(failures) > 0){ #
       aliasGeneIDs <- mget(failures, org.Hs.egALIAS2EG, ifnotfound=NA)
       alias.failures <- names(which(is.na(aliasGeneIDs)))
       alias.success <- aliasGeneIDs[which(!is.na(aliasGeneIDs))]
       if(length(alias.success) > 0)
          geneIDs <- c(geneIDs, alias.success)
       if(length(alias.failures) > 0){
          mouseGeneIDs <- mget(alias.failures, org.Mm.egSYMBOL2EG, ifnotfound=NA)
          mouse.failures <- names(which(is.na(mouseGeneIDs)))
          mouse.success <- mouseGeneIDs[which(!is.na(aliasGeneIDs))]
          if(length(mouse.success) > 0)
             geneIDs <- c(geneIDs, mouse.success)
          if(length(mouse.failures) > 0){
             names(mouse.failures) <- mouse.failures
             mouse.failures [mouse.failures] <- NA
             geneIDs <- c(geneIDs, mouse.failures)
             } # if mouse.failures
          } # if alias.failures
       }# if failures
       
    result <- as.character(sapply (symbols,
                     function(sym) {target <- sprintf("^%s$", toupper(sym));
                                    as.character(geneIDs[grep(target, names(geneIDs))])[1]}))
    names(result) <- symbols
    

    result

} # determineGeneIDs
#-------------------------------------------------------------------------------
createProviderNames <- function(tbl.mAnno)
{
    x <- as.character(apply(tbl.mAnno, 1,
                      function(row) paste(c(row[11], "jolma2013", row[c(1)]), collapse="-")))

    dups <- which(duplicated(x))
    printf("incoming dups count: %d", length(dups))
    
    suffix <- 1
    
    while(length(dups) > 0) {
       suffix <- suffix + 1
       fixed.names <- as.character(sapply(x[dups],
                              function(dup) paste(c(dup, suffix), collapse="-")))
       x[dups] <- fixed.names
       dups <- which(duplicated(x))
       } # while dups

     x <- sub("2$", "2", x)
     x <- sub("2-3$", "3", x)
     x <- sub("2-3-4$", "4", x)
     x <- sub("2-3-4-5$", "5", x)
     x <- sub("2-3-4-5-6$", "6", x)
     x <- sub("2-3-4-5-6-7$", "7", x)
     x <- sub("2-3-4-5-6-7-8$", "8", x)
     x <- sub("2-3-4-5-6-7-8-9$", "9", x)

     x
    

} # createProviderNames
#-------------------------------------------------------------------------------
test.createProviderNames <- function()
{
   print("--- test.createProviderNames")
   dataDir <- file.path(dataDir, "jolma2013")
   x <- readRawMatrices(dataDir)
   tbl.mAnno <- x$anno
   provider.names <- createProviderNames(tbl.mAnno)
   checkEquals(length(duplicated(provider.names)), 0)
   
} # test.createProviderNames
#-------------------------------------------------------------------------------
