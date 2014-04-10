# jaspar2014/test.R
#------------------------------------------------------------------------------------------------------------------------
library (RUnit)
source("import2.R")
#------------------------------------------------------------------------------------------------------------------------
run.tests = function (dataDir)
{
  dataDir <- file.path(dataDir, "jaspar2014")
  x.tbl.rmat <<- test.readRawMatrices (dataDir)
  x.matrices <<- test.extractMatrices (x.tbl.rmat)
  x.tbl.anno <<- test.createAnnotationTable (dataDir)
  test.assignGeneId (dataDir)
  x.tbl.md <<- test.createMetadataTable (x.tbl.anno, x.matrices)
  x.matrices.renamed <<- test.renameMatrices (x.matrices, x.tbl.md, x.tbl.anno)
  x.matrices.normalized <<- test.normalizeMatrices (x.matrices.renamed)


} # run.tests
#------------------------------------------------------------------------------------------------------------------------
test.readRawMatrices = function (dataDir)
{
  print ('--- test.readMatrices')
  tbl.rmat = readRawMatrices (dataDir)
  browser()
  checkTrue (length (tbl.rmat) == 592)
  testMatrice <-tbl.rmat[[2]]
  checkTrue(testMatrice[1]==">MA0006.1 Arnt::Ahr")
    
  invisible (tbl.rmat)

} # test.readRawMatrices
#------------------------------------------------------------------------------------------------------------------------
test.extractMatrices = function (tbl.rmat)
{
  print ('--- test.convertRawMatricesToStandard') 
  matrices = extractMatrices (tbl.rmat)

  checkTrue(length(matrices) == length (tbl.rmat))
  checkTrue(typeof(names(matrices)) == "character")
  checkTrue(sum(colSums(matrices[[1]])) ==120)
  checkTrue(sum(colSums(matrices[[2]])) ==144)
    invisible (matrices)
  
} # test.convertRawMatricesToStandard 
#------------------------------------------------------------------------------------------------------------------------
test.createAnnotationTable = function (dataDir)
{
  print ('--- test.createAnnotationTable')
  tbl.anno = createAnnotationTable (dataDir)
  expected = c ("fullID", "id", "category", "mID", "version", "binder", "speciesID", "proteinID", "family", "tax", "class", "pubmed", "type")
  checkEquals (colnames (tbl.anno), expected)

  checkEquals (head (tbl.anno$fullID),  c ("MA0001.1", "MA0002.1", "MA0003.1", "MA0004.1", "MA0005.1", "MA0006.1"))
  invisible (tbl.anno)

} # test.createAnnotationTable
#------------------------------------------------------------------------------------------------------------------------
test.createMetadataTable = function (tbl.anno, matrices)
{
  print ('--- test.createMetadataTable')
   # try it first with just two matrices
  tbl.md = createMetadataTable (tbl.anno, matrices [1:2])
  checkEquals (dim (tbl.md), c (2, 15))
  checkEquals (colnames (tbl.md), c ("providerName", "providerId", "dataSource", "geneSymbol", "geneId", "geneIdType", 
                                     "proteinId", "proteinIdType", "organism", "sequenceCount", "bindingSequence",
                                     "bindingDomain", "tfFamily", "experimentType", "pubmedID"))
  checkEquals (tbl.md$proteinId, c ('P53762', 'P30561, P53762'))
  checkEquals (tbl.md$proteinIdType, c ('UNIPROT', 'UNIPROT'))

    # now use the whole table
  tbl.md = createMetadataTable (tbl.anno, matrices)
  checkEquals (dim (tbl.md), c (length (matrices), 15))
    # test for proper conversion of speciesID = NA or '-' to Vertebrata

  checkEquals (which (is.na (tbl.md$organism)), integer (0))
  checkEquals (grep ('-', tbl.md$organism), integer (0))
    # Mmusculus-JASPAR_CORE-NF-kappaB-MA0061.1 had 'NA' for proteinID, not <NA>. fixed?
  checkEquals (grep ('NA', tbl.md$proteinId), integer (0))
  invisible (tbl.md)

} # test.createMetadataTable
#------------------------------------------------------------------------------------------------------------------------
test.renameMatrices = function (matrices, tbl.md, tbl.anno)
{
    # try it with just the first two matrices
  matrix.pair = matrices [1:2]
  tbl.md = createMetadataTable (tbl.anno, matrix.pair)
  checkEquals (dim (tbl.md), c (2, 15))
  old.matrix.names = names (matrix.pair)
  matrices.renamed = renameMatrices (matrix.pair, tbl.md)
    # test:  the old name is an id, '9229'.  find, in tbl.anno, the fullID, 'MA0001.1'.  then make sure 'MA000.1' is
    # in the new name of that same matrix

  for (i in 1:length (matrix.pair)) {
    fullID = subset(x.tbl.anno, fullID== substring(old.matrix.names[i],0,8))$fullID
    checkTrue (length (grep (fullID, names (matrices.renamed) [i])) == 1)
    } # for i

    # now try it for the whole set, with selective focused tests
 
  tbl.md = createMetadataTable (tbl.anno, matrices)
  checkEquals (nrow (tbl.md), length (matrices))
  old.matrix.names = names (matrices)
  matrices.renamed = renameMatrices (matrices, tbl.md)
 
  checkEquals (nrow (tbl.md), length (matrices.renamed))
  checkEquals (length (grep ('-MA0', names (matrices.renamed))), length (matrices.renamed))

  invisible (matrices.renamed)

} # test.renameMatrices
#------------------------------------------------------------------------------------------------------------------------
test.convertTaxonCode = function ()
{
  print ('--- test.convertTaxonCode')

  checkEquals (convertTaxonCode ('9606'), 'Hsapiens')
  checkEquals (convertTaxonCode (9606), 'Hsapiens')
     # anomalous codes, which an examination of the jaspar website reveals as 'vertebrates'
  checkEquals (convertTaxonCode (NA), 'Vertebrata')
  checkEquals (convertTaxonCode ('NA'), 'Vertebrata')
  checkEquals (convertTaxonCode (NA_character_), 'Vertebrata')
  checkEquals (convertTaxonCode ('-'), 'Vertebrata')

} # test.convertTaxonCode
#------------------------------------------------------------------------------------------------------------------------
test.guessProteinIdentifierType = function (moleculeName)
{
  print ('--- test.guessProteinIdentifierType')
  checkEquals (guessProteinIdentifierType ('P29383'), 'UNIPROT')

  all.types = sapply (x.tbl.anno$proteinID, guessProteinIdentifierType)
  checkTrue (length (which (is.na (all.types))) < 12)   # got most of them.

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
test.assignGeneId = function (dataDir, proteinId)
{
  print ('--- test.assignGeneId')
  uniprot.ids = c ('Q9GRA5', 'P31314', 'AAC18941', 'O49397')
  refseq.ids  = c ('NP_995315.1', 'NP_032840', 'NP_599022')
  yeast.ids   = c ('YKL112W', 'YMR072W', 'YLR131C')

  checkEquals (assignGeneId ('NP_995315.1'), list (geneId='4782', type='ENTREZ'))
  checkEquals (assignGeneId ('NP_599022'),   list (geneId='6095', type='ENTREZ'))

  checkEquals (assignGeneId ('P31314'),      list (geneId='3195', type='ENTREZ'))

  checkEquals (assignGeneId ('YKL112W'),     list (geneId='YKL112W', type='SGD'))

 

} # test.assignGeneId
#------------------------------------------------------------------------------------------------------------------------
