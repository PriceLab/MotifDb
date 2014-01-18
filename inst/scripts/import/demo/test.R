# example/test.R
#------------------------------------------------------------------------------------------------------------------------
library (RUnit)
#------------------------------------------------------------------------------------------------------------------------
source("import.R")
#------------------------------------------------------------------------------------------------------------------------
run.tests = function ()
{
  matrices.raw <- test.readRawMatrices (dataDir="./")
  matrices <- test.extractAndNormalizeMatrices (matrices.raw)
  #tbl.md <- test.createMetadataTable (matrices)
  #matrices.renamed <- test.renameMatrices (matrices, tbl.md)

} # run.tests
#------------------------------------------------------------------------------------------------------------------------
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
  print ('--- test.readMatrices')
  list.pwms = readRawMatrices (dataDir)
  checkEquals (length (list.pwms), 4)
  checkEquals (names (list.pwms [[1]]), c ("title", "matrix"))
  checkEquals (rownames (list.pwms[[1]]$matrix),  c ("A", "C", "G", "T"))
  invisible (list.pwms)

} # test.readRawMatrices
#------------------------------------------------------------------------------------------------------------------------
test.extractAndNormalizeMatrices = function (matrices.raw)
{
  print ('--- test.extractAndNormalizeMatrices')

  matrices.fixed <- extractAndNormalizeMatrices(matrices.raw)
  checkEquals(length(matrices.fixed), 4)

  checkEquals(names(matrices.fixed),
            c("MA0004.1 Arnt", "MA0006.1 Arnt::Ahr",  "MA0008.1 HAT5",
              "MA0009.1 T"))

    # make sure all columns in all matrices sum to 1.0
  checkTrue (all(sapply(matrices.fixed,
                         function(m)all(abs(colSums(m)-1.0) < 1e-10))))

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
