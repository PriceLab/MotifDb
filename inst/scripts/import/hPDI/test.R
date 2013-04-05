# hPDI/test.R
#------------------------------------------------------------------------------------------------------------------------
library (RUnit)
#------------------------------------------------------------------------------------------------------------------------
source("import.R")
#------------------------------------------------------------------------------------------------------------------------
run.tests = function (dataDir)
{
  dataDir <- file.path(dataDir, "hPDI")
  #freshStart ()

  x.filenames <- test.getMatrixFilenames (dataDir)
  x.matrices <- test.readMatrices (x.filenames)
  x.tbl.anno.raw <- test.readRawAnnotation (dataDir)
  x.refseqs  <- test.uniprotToRefSeq (x.tbl.anno.raw)
  test.extractProteinIDs (dataDir)
  x.tbl.annop.raw <- test.addProteinColumn (x.tbl.anno.raw)
  x.tbl.annof.raw <- test.trimAnnoTable (x.tbl.annop.raw)

    # now, one summary test for the above tests combined into one method
  x.tbl.anno <- test.readAnnotation (dataDir)

  x.tbl.md <- test.createMetadata (dataDir)
  x.mat5f <- test.normalizeMatrices (x.matrices [1:5])
  x.mat5fr <- test.renameMatrices (dataDir, x.filenames) 

} # run.tests
#------------------------------------------------------------------------------------------------------------------------
test.getMatrixFilenames = function  (dataDir)
{
  print("--- test.getMatrixFilenames")
  filenames = getMatrixFilenames (dataDir)
  checkTrue (length (filenames) == 437)
  checkTrue(file.exists(filenames[1]))
  checkTrue(file.exists(filenames[437]))

  invisible (filenames)

} # test.getMatrixFilenames
#------------------------------------------------------------------------------------------------------------------------
test.readMatrices = function (filenames)
{
  print ('--- test.readMatrices')
  matrices = readMatrices (filenames[1:10])
  checkEquals (length (matrices), 10)
  checkEquals(names(matrices)[1], "ABCF2")
  checkEquals(dim(matrices[["ABCF2"]]), c(4,6))
  checkEquals(rownames(matrices[["ABCF2"]]), c("A", "C", "G", "T"))
  invisible (matrices)

} # test.readMatrices
#------------------------------------------------------------------------------------------------------------------------
test.extractProteinIDs = function (dataDir)
{
  print ('--- test.extractProteinIDs')
  
  filename <- file.path(dataDir, "protein_annotation.txt")
  raw.list = tail (read.table (filename, sep='\t', header=T, fill=T, stringsAsFactors=FALSE, quote="")$ensembl_description, n=16)
  protein.ids <<- extractProteinIDs (raw.list)
  checkEquals (length (protein.ids), length (raw.list))
  checkEquals (protein.ids [2], 'NP_001012677')
  checkEquals (protein.ids [5], 'NP_001074306')
  checkTrue (is.na (protein.ids [10]))   # Q9GND9 -> NA
  checkEquals (protein.ids [11], 'NP_001139192') # Q8NAP8 -> ZBT8B ->  NP_001139192.1

} # test.extractProteinIDs
#------------------------------------------------------------------------------------------------------------------------
test.readRawAnnotation = function (dataDir)
{
  print ('--- test.readRawAnnotation')
  tbl.anno = readRawAnnotation (dataDir)
  checkEquals (dim (tbl.anno), c (4191, 16))
    # are all possible proteins mapped to refseq?
  
  invisible (tbl.anno)

} # test.readRawAnnotation
#------------------------------------------------------------------------------------------------------------------------
test.addProteinColumn = function (tbl.anno)
{
  print ('--- test.addProteinColumn')
  tbl.small = tail (tbl.anno, n=16)
  tbl.sp = addProteinColumn (tbl.small)
  checkEquals (ncol (tbl.small) + 1, ncol (tbl.sp))
  checkEquals (nrow (tbl.small), nrow (tbl.sp))
  checkEquals (colnames (tbl.sp) [length (colnames (tbl.sp))], 'protein')
  proteins = tbl.sp$protein
  checkEquals (length (grep ('NP_', proteins)), 3)
  checkEquals (length (which (is.na (proteins))), 13)

  invisible (tbl.sp)

} # test.addProteinColumn
#------------------------------------------------------------------------------------------------------------------------
test.trimAnnoTable = function (tbl.anno)
{
  print ('--- test.trimAnnoTable')
  checkEquals (ncol (tbl.anno), 17)
  tbl.trimmed = trimAnnoTable (tbl.anno)
  checkEquals (ncol (tbl.trimmed), 4)
  checkEquals (colnames (tbl.trimmed), c ('geneSymbol', 'geneID', 'pfamDomain', 'protein'))
  invisible (tbl.trimmed)

} # test.trimAnnoTable
#------------------------------------------------------------------------------------------------------------------------
test.readAnnotation = function (dataDir)
{
  print("--- test.readAnnotation")
  tbl.anno = readAnnotation (dataDir)
  checkEquals (dim (tbl.anno), c (4191, 4))
  checkEquals (colnames (tbl.anno), c ("geneSymbol", "geneID", "pfamDomain", "protein"))

  invisible (tbl.anno)

} # test.readAnnotation
#------------------------------------------------------------------------------------------------------------------------
test.createMetadata = function (dataDir)
{
  print ('--- test.createMetadata')
  tbl.anno = readAnnotation (dataDir)
  filenames = getMatrixFilenames (dataDir)
  matrices.all = readMatrices (filenames)
    # take the top 3, and to be sure of a good test, add two with no protein annotation provided

  names.of.na.protein.matrices = head (subset (tbl.anno, geneSymbol %in% names (matrices.all) & is.na (protein)), n=2)$geneSymbol
  matrices = c (matrices.all [1:3], matrices.all [names.of.na.protein.matrices])
  tbl.md = createMetadata (matrices, tbl.anno)

  checkEquals (nrow (tbl.md), 5) 
  checkEquals (colnames (tbl.md), c ("providerName", "providerId", "dataSource", "geneSymbol", "geneId", "geneIdType", 
                                     "proteinId", "proteinIdType", "organism", "sequenceCount", "bindingSequence",
                                     "bindingDomain", "tfFamily", "experimentType", "pubmedID"))

  checkTrue (all (!is.na (tbl.md$providerName)))
  checkTrue (all (!is.na (tbl.md$providerId)))
  checkTrue (all (!is.na (tbl.md$dataSource)))
  checkTrue (all (!is.na (tbl.md$organism)))
  checkTrue (all (!is.na (tbl.md$sequenceCount)))

  checkTrue (all (is.na (tbl.md$bindingSequence)))

  checkTrue (all (tbl.md$dataSource=='hPDI'))
  checkTrue (is.integer (tbl.md$sequenceCount))
  checkTrue (is.character (tbl.md$geneId))

  syms = tbl.md$geneSymbol
  ids = tbl.md$geneId
  predicted.ids = as.character (mget (syms, org.Hs.egSYMBOL2EG))
  checkEquals (ids, predicted.ids)

  invisible (tbl.md)

} # test.createMetadata
#------------------------------------------------------------------------------------------------------------------------
test.uniprotToRefSeq = function (tbl.anno)
{
  print ('--- test.uniprotToRefSeq')
  samples = tbl.anno$protein [sample (1:nrow (tbl.anno), size=10, replace=FALSE)]
  samples = c ("Q13523", "Q9UER7", "P52738", "Q96Q11", "Q8NCF5", "Q15784", "Q01664", "O14980", "P0C1Z6")
  result = uniprotToRefSeq (samples)
  checkEquals (names (result), samples)
  refseqs = as.character (result)
  checkEquals (length (grep ('^NP_', as.character  (refseqs))), length (samples))

  samples.2 = c ('foo', "Q9UER7")
  result.2 = uniprotToRefSeq (samples.2)
  checkEquals (names (result.2), samples.2)
  checkEquals (as.vector (result.2), c (NA_character_, "NP_001135441"))

  invisible (result)

} # test.uniprotToRefSeq
#------------------------------------------------------------------------------------------------------------------------
test.normalizeMatrices = function (matrices)
{
  print ('--- test.normalizeMatrices')

  colsums = as.integer (sapply (matrices, function (mtx) as.integer (mean (round (colSums (mtx))))))
  checkTrue (all (colsums > 1))

  matrices.norm = normalizeMatrices (matrices)
  checkEquals (names (matrices.norm), names (matrices))

  colsums = as.integer (sapply (matrices.norm, function (mtx) as.integer (mean (round (colSums (mtx))))))
  checkTrue (all (colsums == 1))
  invisible (matrices.norm)

} # test.normalizeMatrices
#------------------------------------------------------------------------------------------------------------------------
test.renameMatrices = function (dataDir, filenames)
{
  print ('--- test.renameMatrices')

  matrices = readMatrices (filenames)
  tbl.anno = readAnnotation (dataDir)
  tbl.md = createMetadata (matrices, tbl.anno)
   
  old.names = names (matrices)
  mtxr = renameMatrices (matrices, tbl.md)
  new.names = names (mtxr)

  checkEquals (new.names, paste ('Hsapiens-hPDI-', old.names, sep=''))
  invisible (mtxr)

} # test.renameMatrices
#------------------------------------------------------------------------------------------------------------------------
