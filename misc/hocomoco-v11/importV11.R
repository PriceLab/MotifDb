# hocomoco.file <- "HOCOMOCOv11_core_HUMAN_mono_jaspar_format.txt"
# metadata.file <- "HOCOMOCOv11_core_annotation_HUMAN_mono.tsv"
# outputFile <- "hocomocov11-core.RData"

hocomoco.file <- "HOCOMOCOv11_full_HUMAN_mono_jaspar_format.txt"
metadata.file <- "HOCOMOCOv11_full_annotation_HUMAN_mono.tsv"
outputFile <- "hocomocov11-secondary.RData"
# MotifDb/inst/scripes/import/HOCOMOCO/import.R
#------------------------------------------------------------------------------------------------------------------------
library(RUnit)
library(RCurl)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
  test_readRawMatrices()
  test_extractMatrices()
  test_createMetaDataTable()
  test_renameMatrices()
  test_normalizeMatrices()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
readRawMatrices <- function (dataDir, count=-1)
{
  # our convention is that there is a shared "dataDir" visible to
  # the importer, and that within that directory there is one
  # subdirectory for each data source.
  # for this example importer, that directory will be <dataDir>/test
  # within which we will look for one small file "sample.pcm"

  filename <- file.path(dataDir, hocomoco.file)

  stopifnot(file.exists(filename))

  all.lines <- scan(filename, what=character(0), sep='\n', quiet=TRUE)
  title.lines <- grep('^>', all.lines)
  title.line.count <- length(title.lines)
  max <- title.line.count - 1
  if(count > 0)
     max <- count

  pwms <- list ()

  for (i in 1:max){ # loop through all motifs in the matrix file, one motif at a time
    start.line <- title.lines[i]
    end.line <- title.lines[i+1] - 1
    new.pwm <- parsePwm(all.lines [start.line:end.line])
    pwms <- c(pwms, list(new.pwm))
    } # for i

  invisible (pwms)

} # readRawMatrices
#------------------------------------------------------------------------------------------------------------------------
test_readRawMatrices <- function()
{
   printf("--- test_readRawMatrices")
   pwms <- readRawMatrices("./", count=3)
   checkEquals(length(pwms), 3)
   pwm.1 <- pwms[[1]]

   checkEquals(pwm.1$title, "AHR_HUMAN.H11MO.0.B")
   mtx <- pwm.1$matrix
   checkEquals(dim(mtx), c(4, 9))
   checkEquals(rownames(mtx), c("A", "C", "G", "T"))
   checkEquals(sum(mtx), 1386)

} # test_readRawMatrices
#------------------------------------------------------------------------------------------------------------------------
# simple operation: for each element in the pwm.list, create an intem in a new list where
# rather than a sublist of title and matrix, each element is simply a matrix, with the element named
extractMatrices <- function(pwm.list)
{
  matrices <- sapply(pwm.list, function (element) element$matrix)
  matrix.names <- sapply(pwm.list, function (element) element$title)
  matrix.names <- sub("^> ", "", matrix.names)
  names(matrices) <- matrix.names
  matrices

} # extractMatrices
#------------------------------------------------------------------------------------------------------------------------
test_extractMatrices <- function()
{
   printf("--- test_extractMatrices")
   pwms <- readRawMatrices("./", count=3)
   matrices <- extractMatrices(pwms)
   checkEquals(length(matrices), 3)
   checkEquals(names(matrices), c("AHR_HUMAN.H11MO.0.B", "AIRE_HUMAN.H11MO.0.C", "ALX1_HUMAN.H11MO.0.B"))

} # test_extractMatrices
#------------------------------------------------------------------------------------------------------------------------
createMetaDataTable <- function(dataDir, matrices, raw.metadata.filename)
{
  filename <- file.path(dataDir, raw.metadata.filename)
    #printf("checking for readable metadata file:")
    #printf("   %s", filename)
  stopifnot(file.exists(filename))

  tbl.raw <- read.table(filename, sep="\t", header=TRUE, as.is=TRUE)
  tbl.md <- data.frame ()
  matrix.ids <- names(matrices)

  for (matrix.id in matrix.ids) {
    tokens <- strsplit(matrix.id, ".", fixed=TRUE)[[1]]
    reliability.score <- tokens[length(tokens)]
    matrix <- matrices[[matrix.id]]
    geneSymbol <- sub("_HUMAN.*$", "", matrix.id)
    tbl.sub <- subset(tbl.raw, Model==matrix.id)
    entrez.gene <- tbl.sub$EntrezGene[1]
    uniprot.id <- tbl.sub$UniProt.ID[1]
    tf.family <- paste(tbl.sub$TF.family[1], tbl.sub$TF.subfamily[1], sep=": ")
    experiment.type <- tbl.sub$Data.source
    dataSource <- sprintf("HOCOMOCOv11%s", reliability.score)
    organism <- "Hsapiens"

    new.row <- list (providerName=matrix.id,
                     providerId=matrix.id, #"HOCOMOCO v8 and ChiPMunk 3.1"
                     dataSource=dataSource,
                     geneSymbol=geneSymbol,
                     geneId=entrez.gene,
                     geneIdType="ENTREZ",
                     proteinId=uniprot.id,
                     proteinIdType="UNIPROT",
                     organism="Hsapiens",
                     sequenceCount=max(colSums(matrix)),
                     bindingSequence=NA_character_,
                     bindingDomain=NA,
                     tfFamily=tf.family,
                     experimentType=experiment.type,
                     pubmedID="23175603")
    #printf("matrix.id: %s", matrix.id);
       # inefficient but acceptable:
    tbl.md <- rbind(tbl.md, data.frame (new.row, stringsAsFactors=FALSE))
    full.name <- sprintf ('%s-%s-%s', organism, dataSource, matrix.id)
    rownames (tbl.md) [nrow (tbl.md)] <- full.name
    } # for matrix.id

  invisible (tbl.md)

} # createMetaDataTable
#------------------------------------------------------------------------------------------------------------------------
test_createMetaDataTable <- function()
{
   printf("--- test_createMetaDataTable")

   pwms <- readRawMatrices("./", count=3)
   matrices <- extractMatrices(pwms)

   tbl.md <- createMetaDataTable("./", matrices, raw.metadata.filename=metadata.file)
   checkEquals(dim(tbl.md), c(3, 15))
     # make sure the reliability score (A-D) is searchably found in the matrix rowname
   checkEquals(length(grep("HOCOMOCOv11[ABCD]", rownames(tbl.md))), 3)
   checkEquals(length(grep("HOCOMOCOv11[ABCD]", tbl.md$dataSource)), 3)

   pwms <- readRawMatrices("./", count=-1)
   matrices <- extractMatrices(pwms)
   tbl.md <- createMetaDataTable("./", matrices, raw.metadata.filename=metadata.file)
   checkEquals(length(grep("HOCOMOCOv11[ABCD]", rownames(tbl.md))), nrow(tbl.md))

} # test_createMetaDataTable
#------------------------------------------------------------------------------------------------------------------------
# we now have metadata rownames, one per matrix, each of which has all (or at least at lot) of information
# used in querying.  For instance, the hocomoco give name for the first matrix is
#   "AHR_HUMAN.H11MO.0.B"
# we transform that here to
#   "Hsapiens-HOCOMOCOv11B-AHR_HUMAN.H11MO.0.B"
renameMatrices <- function (matrices, tbl.md)
{
  stopifnot(length(matrices) == nrow(tbl.md))
  names(matrices) <- rownames(tbl.md)
  invisible(matrices)

} # renameMatrices
#------------------------------------------------------------------------------------------------------------------------
test_renameMatrices <- function()
{
   printf("--- test_renameMatrices")
   pwms <- readRawMatrices("./", count=3)
   matrices <- extractMatrices(pwms)
   checkEquals(names(matrices), c("AHR_HUMAN.H11MO.0.B",
                                  "AIRE_HUMAN.H11MO.0.C",
                                  "ALX1_HUMAN.H11MO.0.B"))
   tbl.md <- createMetaDataTable("./", matrices, raw.metadata.filename=metadata.file)
   matrices.renamed <- renameMatrices(matrices, tbl.md)

   checkEquals(names(matrices.renamed), c("Hsapiens-HOCOMOCOv11B-AHR_HUMAN.H11MO.0.B",
                                          "Hsapiens-HOCOMOCOv11C-AIRE_HUMAN.H11MO.0.C",
                                          "Hsapiens-HOCOMOCOv11B-ALX1_HUMAN.H11MO.0.B"))

} # test_renameMatrices
#------------------------------------------------------------------------------------------------------------------------
normalizeMatrices <- function(matrices)
{
  mtx.normalized <- sapply (matrices,
                           function (mtx) apply (mtx, 2, function (colvector) colvector / sum (colvector)))

  invisible (mtx.normalized)

} # normalizeMatrices
#------------------------------------------------------------------------------------------------------------------------
test_normalizeMatrices <- function()
{
   printf("--- test_normalizeMatrices")

   pwms <- readRawMatrices("./", count=3)
   matrices <- extractMatrices(pwms)
   matrices.norm <- normalizeMatrices(matrices)

     # after normalization, each matrix column should total to 1.0
     # so sum(matrix) should be equals to the number of columns
     # after normalization, the summed row counts should be in the same order,
     # indicating that the structure has been preserved

   for(i in 1:3){
     checkEquals(sum(matrices.norm[[i]]), ncol(matrices.norm[[i]]))
     checkEquals(order(rowSums(matrices[[i]])), order(rowSums(matrices.norm[[i]])))
     } # for i

} # test_normalizeMatrices
#------------------------------------------------------------------------------------------------------------------------
parsePwm <- function (text)
{
  lines <- strsplit (text, '\t')
  stopifnot(length(lines)==5) # title line, one line for each base
  title <- sub(">", "", lines [[1]][1], fixed=TRUE)
  line.count <- length(lines)

  # expect 4 rows, and a number of columns we can discern from
  # the incoming text.
  secondLineParsed <- unlist(strsplit(lines[[2]], " "))
  cols <- length(secondLineParsed)
  result <- matrix (nrow=4, ncol=cols,
                    dimnames=list(c('A','C','G','T'),
                                  as.character(1:cols)))
  # loop over the four lines (for each base respectively)
  row <- 1

  for(i in 2:line.count){
    lineParsed <- unlist(strsplit(lines[[i]], "\t"))
    #print(lineParsed)
    result [row,] <- as.numeric(lineParsed)
    row <- row + 1
    } # for i

  #browser()

  #return (list (title=title, consensus.sequence=consensus.sequence, matrix=result))
  return (list (title=title, matrix=result))

} # parsePwm
#----------------------------------------------------------------------------------------------------
readMetaDataTable <- function()
{
   tbl <- read.table(metadata.file, sep="\t", header=TRUE, as.is=TRUE)
   dim(tbl)

    # Model                                AHR_HUMAN.H11MO.0.B
    # Transcription.factor                                 AHR
    # Model.length                                           9
    # Quality                                                B
    # Model.rank                                             0
    # Consensus                                      dKhGCGTGh
    # Model.release                                 HOCOMOCOv9
    # Data.source                                  Integrative
    # Best.auROC..human.                                  <NA>
    # Best.auROC..mouse.                                  <NA>
    # Peak.sets.in.benchmark..human.                      <NA>
    # Peak.sets.in.benchmark..mouse.                      <NA>
    # Aligned.words                                        157
    # TF.family                      PAS domain factors{1.2.5}
    # TF.subfamily                   Ahr-like factors{1.2.5.1}
    # HGNC                                                 348
    # EntrezGene                                           196
    # UniProt.ID                                     AHR_HUMAN
    # UniProt.AC                                        P35869


    #                     Mmusculus-jaspar2014-Arnt-MA0004
    #     providerName    "MA0004.1 Arnt"
    #     providerId      "MA0004.1 Arnt"
    #     dataSource      "jaspar2014"
    #     geneSymbol      "Arnt"
    #     geneId          NA
    #     geneIdType      NA
    #     proteinId       "P53762"
    #     proteinIdType   "UNIPROT"
    #     organism        "Mmusculus"
    #     sequenceCount   "20"
    #     bindingSequence NA
    #     bindingDomain   NA
    #     tfFamily        "Helix-Loop-Helix"
    #     experimentType  "SELEX"
    #     pubmedID        "24194598"




} # readMetaDataTable
#----------------------------------------------------------------------------------------------------
run <- function (dataDir=".")
{
  rawMatrixList <- readRawMatrices("./", dataDir)
  length(rawMatrixList)
  matrices <- extractMatrices (rawMatrixList)
  length(matrices)
  matrices[1]
  tbl.md <- createMetaDataTable(dataDir, matrices, raw.metadata.filename=metadata.file)
  matrices <- normalizeMatrices(matrices)
  matrices <- renameMatrices(matrices, tbl.md)

  serializedFile <- file.path(dataDir, outputFile)
  printf("writing %s to %s", outputFile, dataDir)

  save(matrices, tbl.md, file=serializedFile)
  printf("saved %d matrices to %s", length(matrices), serializedFile)
  printf("next step:  copy %s to <packageRoot>/MotifDb/inst/extdata, rebuild package", serializedFile)

} # run
#------------------------------------------------------------------------------------------------------------------------
blendCoreAndSecondaryDataSets <- function()
{
   load("hocomocov11-full.Rdata")
   matrices.secondary <- matrices
   tbl.secondary <- tbl.md
   length(matrices.secondary);   # 768
   dim(tbl.secondary)            # 768 15

   checkTrue(all(names(matrices.secondary) %in% rownames(tbl.secondary)))

   load("hocomocov11-core.RData")
   matrices.core <- matrices
   tbl.core <- tbl.md
   checkTrue(all(names(matrices.core) %in% rownames(tbl.core)))
   length(matrices.core);  # 400
   dim(tbl.core)           # 400 15

   checkTrue(all(names(matrices.core) %in% rownames(tbl.secondary)))
   checkTrue(all(names(matrices.core) %in% names(matrices.secondary)))

      # above tests establish that tbl.secondary and matrices.secondary are a proper superset of the core set
      # all we need is to now mark the tbl.secondary$dataSource to distinguish core & secondary

   mapping <- match(rownames(tbl.core), rownames(tbl.secondary))
   checkEquals(length(mapping), 400)
   unmapped.secondary.only <- setdiff(seq_len(nrow(tbl.secondary)), mapping)

   checkEquals(length(unmapped.secondary.only), 368)

     # are all core matrices in secondary, with the same names
   checkTrue(all(names(matrices.core) == names(matrices.secondary)[mapping]))
   checkTrue(all(rownames(tbl.core) == rownames(tbl.secondary)[mapping]))

   nrow(subset(tbl.secondary, geneId==""))  # 3
   rownames(subset(tbl.secondary, geneId==""))
      # [1] "Hsapiens-HOCOMOCOv11D-FOXO6_HUMAN.H11MO.0.D"
      # [2] "Hsapiens-HOCOMOCOv11D-KLF14_HUMAN.H11MO.0.D"
      # [3] "Hsapiens-HOCOMOCOv11D-ZN713_HUMAN.H11MO.0.D"

      # FOXO6 forkhead box O6 [Homo sapiens (human)]
      #   Gene ID: 100132074, updated on 13-Mar-2020
      # KLF14 Kruppel like factor 14 [Homo sapiens (human)]
      #   Gene ID: 136259, updated on 13-Mar-2020
      # ZNF713 zinc finger protein 713 [Homo sapiens (human)],
      #   Gene ID: 349075, updated on 13-Mar-2020

    tbl.secondary["Hsapiens-HOCOMOCOv11D-FOXO6_HUMAN.H11MO.0.D", "geneId"] <- "100132074"
    tbl.secondary["Hsapiens-HOCOMOCOv11D-KLF14_HUMAN.H11MO.0.D", "geneId"] <- "136259"
    tbl.secondary["Hsapiens-HOCOMOCOv11D-ZN713_HUMAN.H11MO.0.D", "geneId"] <- "349075"

    nrow(subset(tbl.secondary, geneId==""))  # was 3, should now be zero

    length(mapping)
    length(unmapped.secondary.only)

     # revise all entries in the dataSource column
     # the dataSource starts off: HOCOMOCOv11[ABCD]
     # we want, in order to support MotifDb::query
     #  HOCOMOCOv11-[core|secondary]-[ABCD]

   x <- tbl.secondary$dataSource
     # pick off just the category: A, B, C, D
   class <- substr(x, 12, 12)

   tbl.secondary$dataSource[unmapped.secondary.only] <- paste0("HOCOMOCOv11-secondary-", class[unmapped.secondary.only])
   tbl.secondary$dataSource[mapping] <-  paste0("HOCOMOCOv11-core-", class[mapping])

     # now revise all of the tbl.md rownames: they need "core" and "secondary"
     # inserted also

   x <- rownames(tbl.md)


   matrices <- matrices.secondary
   tbl.md <- tbl.secondary
   save(matrices, tbl.md, file="hocomoco11.RData")


} # blendCoreAndSecondaryDataSets
#------------------------------------------------------------------------------------------------------------------------
