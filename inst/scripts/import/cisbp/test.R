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
run.tests = function (dataDir=NA)
{
  
  printf("testing cisbp matrix and metadata import");
  if(is.na(dataDir))
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

   filename <- tempfile()
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
