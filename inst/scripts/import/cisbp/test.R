# MotifDb/inst/scripes/import/cispb/test.R
#------------------------------------------------------------------------------------------------------------------------
# Determination and Inference of Eukaryotic Transcription Factor Sequence Specificity
# Weirauch et al, Cell sep 2014
# pmid 25215497
#
# CIS-BP Database: Catalog of Inferred Sequence Binding Preferences
# http://cisbp.ccbr.utoronto.ca/
# discovered that 1460 matrices were empty except for the colnames, leaving 5099 good matrix files
# only 
# 
#------------------------------------------------------------------------------------------------------------------------
library (RUnit)
#------------------------------------------------------------------------------------------------------------------------
source("import.R")
#------------------------------------------------------------------------------------------------------------------------
run.tests = function ()
{
  dataDir <- "/Users/pshannon/s/data/public/TFBS";
    # chosen at random from the 1143 non-duplicate motifs among the 6559 pwm files provided by cisbp
  test.parsePwm(dataDir, filename="M0316_1.02.txt")
  md.sample <- test.getMetadata(dataDir, motifName="M0316_1.02")
  test.translateMetadataToMotifDbStandardForm(md.sample)
  test.createMotifDbArchiveFile(dataDir)

} # run.tests
#------------------------------------------------------------------------------------------------------------------------
# create representative text here, a list of 5 lines, and make sure parsePWM can translate it into
# a one-element list, a 4 x 6 matrix (4 nucleotides by 6 positions) appropriately named.
test.parsePwm <- function(dataDir, filename)
{
   print("--- test.parsePwm")
   title <- sub(".txt", "", filename, fixe=TRUE)
   full.path <- file.path(dataDir, "cisbp", "pwms", filename)
   checkTrue(file.exists(full.path))

   text <- scan(full.path, sep="\n", quiet=TRUE, what=character(0)) 
   pwm <- parsePwm(title, text)
   checkEquals(names(pwm), c("title", "matrix"))
   checkEquals(pwm$title, title)
   m <- pwm$matrix
   checkEquals(dim(m),c(4, 10))
   checkEqualsNumeric(colSums(m), rep(1,10))

} # test.parsePwm
#------------------------------------------------------------------------------------------------------------------------
test.getMetadata <- function(dataDir, motifName)
{
   metadata.file <- file.path(dataDir, "cisbp", "cisbp-metadata-6559-matrices.Rdata")
   checkTrue(file.exists(metadata.file))
   load(metadata.file)
   checkTrue(motifName %in% tbl$motif_id)
   row <- grep(motifName, tbl$motif_id)
   md <- as.list(tbl[row,])
   checkEquals(md, list(ma_id="MA0116870_1.02",
                        tf_id="T025973_1.02",
                        motif_id="M0316_1.02",
                        species="Mus_musculus",
                        TF_Name="Nfil3",
                        PMID="25215497",
                        Family_Name="bZIP",
                        DBID="ENSMUSP00000065363"))
   printf("%s", TRUE)
   invisible(md)

} # test.getMetadata
#------------------------------------------------------------------------------------------------------------------------
test.translateMetadataToMotifDbStandardForm <- function(metadata.element)
{
   checkTrue(is.list(metadata.element))
   x <- translateMetadataToMotifDbStandardForm(metadata.element)
   checkTrue(is.list(x))
   checkTrue(all(c("providerName", "providerId", "dataSource", "geneSymbol", "geneId", "geneIdType",
                   "proteinId", "proteinIdType", "organism", "sequenceCount", "bindingSequence",
                   "bindingDomain", "tfFamily", "experimentType", "pubmedID") %in% names(x)))

   checkEquals(x$providerName,  "M0316_1.02")
   checkEquals(x$providerId,    "MA0116870_1.02")
   checkEquals(x$dataSource,    "cispb_1.02")
   checkEquals(x$geneSymbol,    "Nfil3")
   checkEquals(x$proteinId,     "ENSMUSP00000065363")
   checkEquals(x$proteinIdType, "ensembl")
   checkEquals(x$organism,      "Mmusculus")
   checkEquals(x$tfFamily,      "bZIP")
   checkEquals(x$pubmedID,      "25215497")

   checkTrue(is.na(x$geneId))
   checkTrue(is.na(x$geneIdType))
   checkTrue(is.na(x$sequenceCount))
   checkTrue(is.na(x$bindingSequence))
   checkTrue(is.na(x$bindingDomain))
   checkTrue(is.na(x$experimentType))

} # test.translateMetadataToMotifDbStandardForm
#------------------------------------------------------------------------------------------------------------------------
test.createMotifDbArchiveFile <- function(dataDir)
{
   print("--- test.createMotifDbArchiveFile")
   
     # we can test with a few filenames from the the specified <dataDir>, or
     # (if <count> is not given, use all of them.
     # be forewarned: of the 5099 non-empty matrix files offered by cisbp, 

   filename <- "tmp10.RData"
   createMotifDbArchiveFile(dataDir, filename, count=10)
   load(filename)
   checkEquals(dim(tbl.md), c(10,15))
   checkEquals(length(matrices), 10)
   checkTrue(all(rownames(tbl.md) == names(matrices)))
   checkTrue(all(as.numeric(colSums(matrices[[1]])) - 1.0 < 1e-8))

        # spot check two  matrices.  as pwms, all columns should sum to very close to 1.0
   checkTrue(all(as.numeric(colSums(matrices[[1]])) - 1.0 < 1e-8))
   checkTrue(all(as.numeric(colSums(matrices[[7]])) - 1.0 < 1e-8))

} # test.createMotifDbArchiveFile
#------------------------------------------------------------------------------------------------------------------------
## test.readRawMatrices = function (dataDir)
## {
##   print ('--- test.readRawMatrices')
##   list.pwms = readRawMatrices (dataDir)
##   checkEquals (length (list.pwms), 4)
##   checkEquals (names (list.pwms [[1]]), c ("title", "matrix"))
##   checkEquals (rownames (list.pwms[[1]]$matrix),  c ("A", "C", "G", "T"))
##   invisible (list.pwms)
## 
## } # test.readRawMatrices
## #------------------------------------------------------------------------------------------------------------------------
## test.extractMatrices = function (matrices.raw)
## {
##   print ('--- test.extractMatrices')
## 
##   matrices.fixed <- extractMatrices(matrices.raw)
##   checkEquals(length(matrices.fixed), 4)
## 
##   checkEquals(names(matrices.fixed),
##             c("MA0004.1 Arnt", "MA0006.1 Arnt::Ahr",  "MA0008.1 HAT5",
##               "MA0009.1 T"))
## 
##     # make sure all columns in all matrices sum to 1.0
##   #checkTrue (all(sapply(matrices.fixed,
##   #                       function(m)all(abs(colSums(m)-1.0) < 1e-10))))
## 
##   invisible (matrices.fixed)
## 
## } # test.extractMatrices
## #------------------------------------------------------------------------------------------------------------------------
## test.normalizeMatrices = function (matrices)
## {
##   print ('--- test.normalizeMatrices')
## 
##   matrices.fixed <- normalizeMatrices(matrices)
##   checkEquals(length(matrices.fixed), 4)
## 
##   checkEquals(names(matrices.fixed),
##             c("MA0004.1 Arnt", "MA0006.1 Arnt::Ahr",  "MA0008.1 HAT5",
##               "MA0009.1 T"))
## 
##     # make sure all columns in all matrices sum to 1.0
##   checkTrue (all(sapply(matrices.fixed,
##                          function(m)all(abs(colSums(m)-1.0) < 1e-10))))
## 
##   invisible (matrices.fixed)
## 
## } # test.extracNormalizeMatrices
## #------------------------------------------------------------------------------------------------------------------------
## test.createMetadataTable = function (dataDir, matrices, raw.metadata.filename)
## {
##     print ('--- test.createMetadataTable')
## 
##        # try it with just the first matrix
##     tbl.md = createMetadataTable (dataDir, matrices, raw.metadata.filename)
## 
##     checkEquals (dim (tbl.md), c (length(matrices), 15))
##     checkEquals (colnames (tbl.md), c ("providerName", "providerId", "dataSource", "geneSymbol", "geneId", "geneIdType", 
##                                        "proteinId", "proteinIdType", "organism", "sequenceCount", "bindingSequence",
##                                        "bindingDomain", "tfFamily", "experimentType", "pubmedID"))
## 
##     with(tbl.md[1,], checkEquals(providerName,   "MA0004.1 Arnt"),
##                      checkEquals(providerId,     "MA0004.1 Arnt"),
##                      checkEquals(dataSource,     "jaspar2014"),
##                      checkEquals(geneSymbol,     "Arnt"),
##                      checkEquals(proteinId,      "P53762"),
##                      checkEquals(proteinIdType,  "uniprot"),
##                      checkEquals(organism,       "Mus musculus"),
##                      checkEquals(sequenceCount,  20),
##                      checkEquals(tfFamily,       "Helix-Loop-Helix"),
##                      checkEquals(experimentType, "SELEX"),
##                      checkEquals(pubmedID,       "24194598"),
##                      checkTrue(is.na(geneId)),
##                      checkTrue(is.na(geneIdType)),
##                      checkTrue(is.na(bindingSequence)),
##                      checkTrue(is.na(bindingDomain)))
##              
##     invisible (tbl.md)
## 
## } # test.createMetadataTable
## #------------------------------------------------------------------------------------------------------------------------
## test.renameMatrices = function (matrices, tbl.md)
## {
##   print("--- test.renameMatrices")
##   
##     # try it with just the first two matrices
##   matrix.pair <- matrices[1:2]
##   tbl.pair <- tbl.md[1:2,]
##   matrix.pair.renamed <- renameMatrices (matrix.pair, tbl.pair)
##   checkEquals (names (matrix.pair.renamed), 
##                c("Mus musculus-jaspar2014-MA0004.1 Arnt",
##                  "Mus musculus-jaspar2014-MA0006.1 Arnt::Ahr"))
## 
## } # test.renameMatrices
## #------------------------------------------------------------------------------------------------------------------------
## test.guessProteinIdentifierType <- function()
## {
##   print ('--- test.guessProteinIdentifierType')
## 
##   checkEquals (guessProteinIdentifierType ('P29383'), 'UNIPROT')
##   checkEquals (guessProteinIdentifierType ('NP_234234'), 'REFSEQ')
## 
## } # test.guessProteinIdentifierType
## #------------------------------------------------------------------------------------------------------------------------
## test.normalizeMatrices = function (matrices)
## {
##   print ('--- test.normalizeMatrices')
## 
##   colsums = as.integer (sapply (matrices, function (mtx) as.integer (mean (round (colSums (mtx))))))
##   #checkTrue (all (colsums > 1))
## 
##   matrices.norm = normalizeMatrices (matrices)
## 
##   colsums = as.integer (sapply (matrices.norm, function (mtx) as.integer (mean (round (colSums (mtx))))))
##   checkTrue (all (colsums == 1))
## 
##   invisible (matrices.norm)
## 
## } # test.normalizeMatrices
## #------------------------------------------------------------------------------------------------------------------------
## test.ncbiTaxonimicCodeToLinnaean <- function()
## {
##     print("--- test.ncbiTaxonimicCodeToLinnaean")
## 
##     checkTrue(is.na(ncbiTaxonimicCodeToLinnaean("purplePeopleEater32")))
## 
##     checkEquals(ncbiTaxonimicCodeToLinnaean("9606"), "Homo sapiens")
##     checkEquals(ncbiTaxonimicCodeToLinnaean(9606), "Homo sapiens")
##     checkEquals(ncbiTaxonimicCodeToLinnaean("10090"), "Mus musculus")
## 
## } # test.ncbiTaxonimicCodeToLinnaean
## #-------------------------------------------------------------------------------
## 
