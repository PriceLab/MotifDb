# jaspar/test.R
# 
#------------------------------------------------------------------------------------------------------------------------
library (RUnit)
source("import.R")
#------------------------------------------------------------------------------------------------------------------------
run.tests = function (dataDir)
{
  dataDir <- file.path(dataDir, "jaspar")
  x.tbl.rmat <<- test.readRawMatrices (dataDir)
  x.matrices <<- test.convertRawMatricesToStandard (x.tbl.rmat)
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
  checkEquals (ncol (tbl.rmat), 4)
  checkEquals (colnames (tbl.rmat), c ('id', 'base', 'pos', 'count'))
  checkEquals (class (tbl.rmat$id),    'character')
  checkEquals (class (tbl.rmat$base),  'character')
  checkEquals (class (tbl.rmat$pos),   'numeric')
  checkEquals (class (tbl.rmat$count), 'numeric')

  checkTrue (nrow (tbl.rmat) > 18000)   # about 450 motifs, each represeted by 4 rows (ACGT) and about 10 positions
  checkTrue (length (unique (tbl.rmat$id)) > 450)
  invisible (tbl.rmat)

} # test.readRawMatrices
#------------------------------------------------------------------------------------------------------------------------
test.convertRawMatricesToStandard = function (tbl.rmat)
{
  print ('--- test.convertRawMatricesToStandard')
   # get just the first two raw matrices

  first.two.ids = head (unique (tbl.rmat$id), n=2)
  rows = nrow (subset (tbl.rmat, id %in% first.two.ids))
  matrices = convertRawMatricesToStandard (tbl.rmat [1:rows,])
  checkEquals (length (matrices), 2)
  checkEquals (names (matrices), first.two.ids)

    # it will not always be true, but IS true for the first two matrices, currently "9229" and  "9231", that there
    # are an equal number of nucleotides at each position.
  checkTrue (all (colSums (matrices [[1]]) == 97))
  checkTrue (all (colSums (matrices [[2]]) == 185))

    # now run all the matrices through
  matrices = convertRawMatricesToStandard (tbl.rmat)
  checkEquals (length (matrices), 459)
  checkEquals (names (matrices)[1:2], first.two.ids)
  
  invisible (matrices)
  
} # test.convertRawMatricesToStandard 
#------------------------------------------------------------------------------------------------------------------------
test.createAnnotationTable = function (dataDir)
{
  print ('--- test.createAnnotationTable')
  tbl.anno = createAnnotationTable (dataDir)
  checkEquals (dim (tbl.anno), c (513, 13))
  expected = c ("fullID", "id", "category", "mID", "version", "binder", "speciesID", "proteinID", "family", "tax", "class", "pubmed", "type")
  checkEquals (colnames (tbl.anno), expected)

  checkEquals (head (tbl.anno$fullID),  c ("MA0001.1", "MA0003.1", "MA0004.1", "MA0005.1", "MA0006.1", "MA0006.1"))
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
  checkEquals (tbl.md$proteinId, c ('P29383', 'P05549'))
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
    fullID = subset (x.tbl.anno, id==old.matrix.names [i])$fullID
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

    # see how successful this is over all 513 proteinIds

  tbl.anno = createAnnotationTable (dataDir)
  mtx.geneId = as.data.frame (t (sapply (tbl.anno$proteinID, assignGeneId)))
  tbl.types = as.data.frame (table (as.character (mtx.geneId$type), useNA='always'), stringsAsFactors=FALSE)
  checkEquals (tbl.types$Var1, c ("ENTREZ", "SGD", NA))
  checkEquals (tbl.types$Freq, c (142, 177, 194))

} # test.assignGeneId
#------------------------------------------------------------------------------------------------------------------------
