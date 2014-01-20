# MotifDb/inst/scripes/import/demo/test.R
#------------------------------------------------------------------------------------------------------------------------
library (RUnit)
#------------------------------------------------------------------------------------------------------------------------
source("import.R")
#------------------------------------------------------------------------------------------------------------------------
run.tests = function ()
{
  test.parsePwm()
  matrices.raw <- test.readRawMatrices (dataDir="./")
  matrices <- test.extractMatrices(matrices.raw)

  test.ncbiTaxonimicCodeToLinnaean()
  test.guessProteinIdentifierType()

  tbl.md <- test.createMetadataTable(dataDir="./", matrices, "md-raw.tsv")
  matrices <- test.normalizeMatrices(matrices)
  matrices.renamed <- test.renameMatrices(matrices, tbl.md)

} # run.tests
#------------------------------------------------------------------------------------------------------------------------
# create representative text here, a list of 5 lines, and make sure parsePWM can translate it into
# a one-element list, a 4 x 6 matrix (4 nucleotides by 6 positions) appropriately named.
test.parsePwm <- function()
{
    print("--- test.parsePwm")
      # a sample matrix, obtained as if from
      # all.lines <- scan (filename, what=character(0), sep='\n', quiet=TRUE)
      # text <- all.lines[1:5]
    text <-  c(">MA0006.1 Arnt::Ahr",
               "3\t0\t0\t0\t0\t0",
               "8\t0\t23\t0\t0\t0",
               "2\t23\t0\t23\t0\t24",
               "11\t1\t1\t1\t24\t0")
    pwm <- parsePwm(text)
    checkEquals(names(pwm), c("title", "matrix"))
    checkEquals(pwm$title, ">MA0006.1 Arnt::Ahr")
    m <- pwm$matrix
    checkEquals(dim(m),c(4, 6))
    checkTrue(all(as.numeric(colSums(m))==24))

} # test.parsePwm
#------------------------------------------------------------------------------------------------------------------------
test.readRawMatrices = function (dataDir)
{
  print ('--- test.readRawMatrices')
  list.pwms = readRawMatrices (dataDir)
  checkEquals (length (list.pwms), 4)
  checkEquals (names (list.pwms [[1]]), c ("title", "matrix"))
  checkEquals (rownames (list.pwms[[1]]$matrix),  c ("A", "C", "G", "T"))
  invisible (list.pwms)

} # test.readRawMatrices
#------------------------------------------------------------------------------------------------------------------------
test.extractMatrices = function (matrices.raw)
{
  print ('--- test.extractMatrices')

  matrices.fixed <- extractMatrices(matrices.raw)
  checkEquals(length(matrices.fixed), 4)

  checkEquals(names(matrices.fixed),
            c("MA0004.1 Arnt", "MA0006.1 Arnt::Ahr",  "MA0008.1 HAT5",
              "MA0009.1 T"))

    # make sure all columns in all matrices sum to 1.0
  #checkTrue (all(sapply(matrices.fixed,
  #                       function(m)all(abs(colSums(m)-1.0) < 1e-10))))

  invisible (matrices.fixed)

} # test.extractMatrices
#------------------------------------------------------------------------------------------------------------------------
test.normalizeMatrices = function (matrices)
{
  print ('--- test.normalizeMatrices')

  matrices.fixed <- normalizeMatrices(matrices)
  checkEquals(length(matrices.fixed), 4)

  checkEquals(names(matrices.fixed),
            c("MA0004.1 Arnt", "MA0006.1 Arnt::Ahr",  "MA0008.1 HAT5",
              "MA0009.1 T"))

    # make sure all columns in all matrices sum to 1.0
  checkTrue (all(sapply(matrices.fixed,
                         function(m)all(abs(colSums(m)-1.0) < 1e-10))))

  invisible (matrices.fixed)

} # test.extracNormalizeMatrices
#------------------------------------------------------------------------------------------------------------------------
test.createMetadataTable = function (dataDir, matrices, raw.metadata.filename)
{
    print ('--- test.createMetadataTable')

       # try it with just the first matrix
    tbl.md = createMetadataTable (dataDir, matrices, raw.metadata.filename)

    checkEquals (dim (tbl.md), c (length(matrices), 15))
    checkEquals (colnames (tbl.md), c ("providerName", "providerId", "dataSource", "geneSymbol", "geneId", "geneIdType", 
                                       "proteinId", "proteinIdType", "organism", "sequenceCount", "bindingSequence",
                                       "bindingDomain", "tfFamily", "experimentType", "pubmedID"))

    with(tbl.md[1,], checkEquals(providerName,   "MA0004.1 Arnt"),
                     checkEquals(providerId,     "MA0004.1 Arnt"),
                     checkEquals(dataSource,     "jaspar2014"),
                     checkEquals(geneSymbol,     "Arnt"),
                     checkEquals(proteinId,      "P53762"),
                     checkEquals(proteinIdType,  "uniprot"),
                     checkEquals(organism,       "Mus musculus"),
                     checkEquals(sequenceCount,  20),
                     checkEquals(tfFamily,       "Helix-Loop-Helix"),
                     checkEquals(experimentType, "SELEX"),
                     checkEquals(pubmedID,       "24194598"),
                     checkTrue(is.na(geneId)),
                     checkTrue(is.na(geneIdType)),
                     checkTrue(is.na(bindingSequence)),
                     checkTrue(is.na(bindingDomain)))
             
    invisible (tbl.md)

} # test.createMetadataTable
#------------------------------------------------------------------------------------------------------------------------
test.renameMatrices = function (matrices, tbl.md)
{
  print("--- test.renameMatrices")
  
    # try it with just the first two matrices
  matrix.pair <- matrices[1:2]
  tbl.pair <- tbl.md[1:2,]
  matrix.pair.renamed <- renameMatrices (matrix.pair, tbl.pair)
  checkEquals (names (matrix.pair.renamed), 
               c("Mus musculus-jaspar2014-MA0004.1 Arnt",
                 "Mus musculus-jaspar2014-MA0006.1 Arnt::Ahr"))

} # test.renameMatrices
#------------------------------------------------------------------------------------------------------------------------
test.guessProteinIdentifierType <- function()
{
  print ('--- test.guessProteinIdentifierType')

  checkEquals (guessProteinIdentifierType ('P29383'), 'UNIPROT')
  checkEquals (guessProteinIdentifierType ('NP_234234'), 'REFSEQ')

} # test.guessProteinIdentifierType
#------------------------------------------------------------------------------------------------------------------------
test.normalizeMatrices = function (matrices)
{
  print ('--- test.normalizeMatrices')

  colsums = as.integer (sapply (matrices, function (mtx) as.integer (mean (round (colSums (mtx))))))
  #checkTrue (all (colsums > 1))

  matrices.norm = normalizeMatrices (matrices)

  colsums = as.integer (sapply (matrices.norm, function (mtx) as.integer (mean (round (colSums (mtx))))))
  checkTrue (all (colsums == 1))

  invisible (matrices.norm)

} # test.normalizeMatrices
#------------------------------------------------------------------------------------------------------------------------
test.ncbiTaxonimicCodeToLinnaean <- function()
{
    print("--- test.ncbiTaxonimicCodeToLinnaean")

    checkTrue(is.na(ncbiTaxonimicCodeToLinnaean("purplePeopleEater32")))

    checkEquals(ncbiTaxonimicCodeToLinnaean("9606"), "Homo sapiens")
    checkEquals(ncbiTaxonimicCodeToLinnaean(9606), "Homo sapiens")
    checkEquals(ncbiTaxonimicCodeToLinnaean("10090"), "Mus musculus")

} # test.ncbiTaxonimicCodeToLinnaean
#-------------------------------------------------------------------------------
