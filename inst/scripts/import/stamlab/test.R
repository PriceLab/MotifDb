# stamlab/test.R
# notes:  3 matrices come w/o speciesID, tax = 'vertebrates'.  not our problem to fix, at least not yet.
#         TBP, HNF4A and CEBPA (MA0108.2, MA0114.1, MA0102.2)
#------------------------------------------------------------------------------------------------------------------------
library (RUnit)
#------------------------------------------------------------------------------------------------------------------------
source("import.R")
#------------------------------------------------------------------------------------------------------------------------
run.tests = function (dataDir)
{
  dataDir <- file.path(dataDir, "stamlab")
  x.rawMatrixList <<- test.readRawMatrices (dataDir)
  x.novels <<- test.readNovelStatus (dataDir)
  x.matrices <<- test.extractAndNormalizeMatrices (x.rawMatrixList)
  x.tbl.md <<- test.createMetadataTable (x.matrices, x.novels)
  x.matrices.renamed <<- test.renameMatrices (x.matrices, x.tbl.md)

} # run.tests
#------------------------------------------------------------------------------------------------------------------------
test.readRawMatrices = function (dataDir)
{
  print ('--- test.readMatrices')
  list.pwms = readRawMatrices (dataDir)
  checkEquals (length (list.pwms), 683)
  checkEquals (names (list.pwms [[1]]), c ("title", "consensus.sequence", "matrix"))
  checkEquals (rownames (list.pwms[[1]]$matrix),  c ("A", "C", "G", "T"))
  invisible (list.pwms)

} # test.readRawMatrices
#------------------------------------------------------------------------------------------------------------------------
test.readNovelStatus = function (dataDir)
{
  print ('--- test.readNovelStatus')
  novel.status = readNovelStatus (dataDir)
  checkEquals (length (novel.status), 683)
  checkEquals (length (which (novel.status == TRUE)), 289)
    # do a spot check around first novel in novels.txt
  checkEquals (as.logical (novel.status [c ('UW.Motif.0010', 'UW.Motif.0011', 'UW.Motif.0012')]),
                             c (FALSE, FALSE, TRUE))
  invisible (novel.status)

} # readNovelStatus
#------------------------------------------------------------------------------------------------------------------------
test.extractAndNormalizeMatrices = function (x.rawMatrixList)
{
  print ('--- test.extractAndNormalizeMatrices')
  matrices.fixed <<- extractAndNormalizeMatrices (x.rawMatrixList)
    # make sure a UW.Motif.0nnn name accompanies each matrix
  checkEquals (length (grep ('UW.Motif.0', names (matrices.fixed))), length (matrices.fixed))
    # make sure all columns in all matrices sum to 1.0
  checkTrue (all (sapply (matrices.fixed, function (m) all (abs (colSums (m) - 1.0) < 1e-10))))
  invisible (matrices.fixed)

} # test.extractAndNormalizeMatrices
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
test.createAnnotationTable = function ()
{
  print ('--- test.createAnnotationTable')
  tbl.anno = createAnnotationTable ()
  checkEquals (dim (tbl.anno), c (513, 13))
  expected = c ("fullID", "id", "category", "mID", "version", "binder", "speciesID", "proteinID", "family", "tax", "class", "pubmed", "type")
  checkEquals (colnames (tbl.anno), expected)

  checkEquals (head (tbl.anno$fullID),  c ("MA0001.1", "MA0003.1", "MA0004.1", "MA0005.1", "MA0006.1", "MA0006.1"))
  invisible (tbl.anno)

} # test.createAnnotationTable
#------------------------------------------------------------------------------------------------------------------------
test.createMetadataTable = function (x.matrices, x.novels)
{
   print ('--- test.createMetadataTable')
    # try it first with just two matrices
   tbl.md = createMetadataTable (x.matrices [1:12], x.novels [1:12])
   checkEquals (dim (tbl.md), c (12, 15))
   checkEquals (colnames (tbl.md), c ("providerName", "providerId", "dataSource", "geneSymbol", "geneId", "geneIdType", 
                                      "proteinId", "proteinIdType", "organism", "sequenceCount", "bindingSequence",
                                      "bindingDomain", "tfFamily", "experimentType", "pubmedID"))
   checkEquals (tbl.md$providerName [1:2], c ('UW.Motif.0001', 'UW.Motif.0002'))
   checkEquals (tbl.md$providerId [1:2], c ('UW.Motif.0001', 'UW.Motif.0002'))
   checkEquals (tbl.md$pubmedID [1:2], c ('22959076', '22959076'))
   checkEquals (tbl.md$dataSource [1:2], c ('stamlab', 'stamlab'))
   checkEquals (tbl.md$organism [1:2], c ('Hsapiens', 'Hsapiens'))
   checkEquals (tbl.md$experimentType [1:2], c ('digital genomic footprinting', 'digital genomic footprinting'))
   checkEquals (tbl.md$geneId, c (rep ('knownMotif', 11), 'novelMotif'))
   checkTrue(all(is.na(tbl.md$geneIdType)))

   invisible (tbl.md)

} # test.createMetadataTable
#------------------------------------------------------------------------------------------------------------------------
test.renameMatrices = function (matrices, tbl.md)
{
  print("--- test.renameMatrices")
  
    # try it with just the first two matrices
  matrix.pair = matrices [1:2]
  tbl.pair = tbl.md [1:2,]
  matrix.pair.renamed = renameMatrices (matrix.pair, tbl.pair)
  checkEquals (names (matrix.pair.renamed), c ("Hsapiens-stamlab-UW.Motif.0001", "Hsapiens-stamlab-UW.Motif.0002"))

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
test.assignGeneId = function (proteinId)
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

  tbl.anno = createAnnotationTable ()
  mtx.geneId = as.data.frame (t (sapply (tbl.anno$proteinID, assignGeneId)))
  tbl.types = as.data.frame (table (as.character (mtx.geneId$type), useNA='always'), stringsAsFactors=FALSE)
  checkEquals (tbl.types$Var1, c ("ENTREZ", "SGD", NA))
  checkEquals (tbl.types$Freq, c (141, 177, 195))

} # test.assignGeneId
#------------------------------------------------------------------------------------------------------------------------
test.parsePwm = function ()
{
  print ('--- test.parsePwm')
  lines = c ('UW.Motif.0006	aggaaatg',
             '0.890585	0.007855	0.051323	0.050237',
             '0.060732	0.004506	0.894170	0.040593',
             '0.072765	0.037935	0.860704	0.028596',
             '0.929585	0.024037	0.034914	0.011464',
             '0.931220	0.023231	0.029078	0.016471',
             '0.857934	0.044211	0.072594	0.025261',
             '0.065840	0.013777	0.058013	0.862370',
             '0.049937	0.036238	0.861871	0.051953')
  m6 = parsePwm (lines)
  checkEquals (names (m6), c ("title", "consensus.sequence", "matrix"))
  pwm = m6$matrix
  checkEquals (dim (pwm), c (4, 8))
  checkEquals (rownames (pwm), c ('A', 'C', 'G', 'T'))
  invisible (m6)

} # test.parsePwm
#------------------------------------------------------------------------------------------------------------------------
