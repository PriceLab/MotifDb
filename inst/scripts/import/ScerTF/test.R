# ScerTF/test.R
#------------------------------------------------------------------------------------------------------------------------
library (RUnit)
library (org.Sc.sgd.db)
#------------------------------------------------------------------------------------------------------------------------
source("import.R")
#------------------------------------------------------------------------------------------------------------------------
run.tests = function (dataDir)
{
  dataDir <- file.path(dataDir, "ScerTF")
  freshStart ()
  x.filenames <<- test.getMatrixFilenames (dataDir)
  txxa <<- test.createMatrixNameUniqifier ()
  txx0 <<- test.toOrf ()
  txx1 <<- test.toUniprot ()
  txx2 <<- test.createExperimentRefTable ()
  txx3 <<- test.parsePWMfromText (dataDir)
  x.matrices <<- test.readAndParse (dataDir)
  x.tbl.md  <<- test.createMetadata (dataDir)
  x.matrices.renamed <<- test.renameMatrices (dataDir, x.matrices, x.tbl.md)

} # run.tests
#------------------------------------------------------------------------------------------------------------------------
test.createExperimentRefTable = function ()
{
  print ('--- test.createExperimentRefTable')
  tbl.ref = createExperimentRefTable ()
  checkEquals (dim (tbl.ref), c (10, 5))
  checkEquals (colnames (tbl.ref), c ('author', 'year', 'pmid', 'organism', 'titles'))

     # make sure neither author nor pmid repeat
  checkEquals (length (unique (tbl.ref$author)), nrow (tbl.ref))
  checkEquals (length (unique (tbl.ref$pmid)), nrow (tbl.ref))

  invisible (tbl.ref)

} # test.createExperiementsTable
#------------------------------------------------------------------------------------------------------------------------
test.parsePWMfromText = function (dataDir)
{
  print ('--- test.parsePWMfromText')

  file <- file.path(dataDir, 'macisaac.ABF1')
  checkTrue(file.exists(file))
  
  lines.of.text = scan (file, what=character(0), sep='\n', quiet=TRUE)
  pwm.abf1 = parsePWMfromText (lines.of.text)
  checkEquals (dim (pwm.abf1), c (4, 15))
  checkEquals (colnames (pwm.abf1), as.character (1:15))
  checkEquals (rownames (pwm.abf1), c ('A', 'C', 'G', 'T'))
  
} # test.parsePWMfromText
#-----------------------------------------------------------------------------------------------------------------------
test.readAndParse = function (dataDir)
{
  print ('--- test.readAndParse')

  all.files = getMatrixFilenames (dataDir)
  sample.file.1 = grep ('badis.ABF2', all.files)
  sample.file.2 = grep ('badis.CAT8', all.files)
  checkEquals (length (sample.file.1), 1)
  checkEquals (length (sample.file.2), 1)
  
  mtx.test = readAndParse (file.path(dataDir, all.files [c (sample.file.1, sample.file.2)]))

  checkEquals (length (mtx.test), 2)
  checkEquals (names (mtx.test), c ("badis.ABF2", "badis.CAT8"))

  checkTrue (all (colSums (mtx.test [[1]]) == 1))
  checkTrue (all (colSums (mtx.test [[2]]) == 1))

  checkEquals (dim (mtx.test [[1]]), c (4,6))
  checkEquals (dim (mtx.test [[2]]), c (4,6))

  invisible (mtx.test)

} # test.readAndParse
#-----------------------------------------------------------------------------------------------------------------------
test.toOrf = function ()
{
  print ('--- test.toOrf')
     # SUT1 is the proper gene symbol for YGL162W, and an alias for YMR125W.  mget returns both.  we want just the first
  checkEquals (toOrf ('SUT1'), 'YGL162W') 
  checkEquals (toOrf (c ( "STB5", "SUT1", "THO2")), c ('YHR178W', 'YGL162W', 'YNL139C'))
  
  checkEquals (toOrf ('bogus', quiet=TRUE), 'bogus')
  checkEquals (toOrf (c ( "STB5", "SUT1", "bogus", "THO2"), quiet=TRUE), c ('YHR178W', 'YGL162W', 'bogus', 'YNL139C'))

} # test.toOrf 
#------------------------------------------------------------------------------------------------------------------------
test.toUniprot = function ()
{
  print ('--- test.toUniprot')
    # demonstrate the problem
  checkEquals (unlist (unname (mget ('YCR039C', org.Sc.sgdUNIPROT))), c ("P0CY08", "P0CY09"))

    # but we want only one
  checkEquals (toUniprot (c ('YCR039C')), "P0CY08")

    # make sure it works embedded in a list of orfs
  checkEquals (toUniprot (c ('YOL108C', 'YCR039C', 'YML065W')), c ("P13902", "P0CY08", "P54784"))

    # an unrecognized orf should return NA
  checkTrue (is.na (toUniprot ('bogus')))

} # test.toUniprot
#------------------------------------------------------------------------------------------------------------------------
test.createMetadata = function (dataDir)
{
  print ('--- test.createMetadata')

  matrices = readAndParse (file.path(dataDir, getMatrixFilenames (dataDir)))

  tbl.md = createMetadata (matrices, createExperimentRefTable ())
  checkEquals (dim (tbl.md), c (length (matrices), 15))
  expected.columns = c ("providerName", "providerId", "dataSource", "geneSymbol", "geneId", "geneIdType", "proteinId", 
                        "proteinIdType", "organism", "sequenceCount", "bindingSequence", "bindingDomain", "tfFamily",
                        "experimentType", "pubmedID")

  checkEquals (colnames (tbl.md), expected.columns)
  checkEquals (unique (tbl.md$organism), 'Scerevisiae')
  ecm = 'Scerevisiae-ScerTF-ECM23-badis'
  checkTrue (ecm %in% rownames (tbl.md))
  x = tbl.md [ecm,]
  checkEquals (x$providerName, "badis.ECM23")
  checkEquals (x$geneSymbol, "ECM23")
  checkEquals (x$geneId,  "YPL021W")
  checkEquals (x$geneIdType, "SGD")
  checkEquals (x$proteinId, "Q02710")
  checkEquals (x$proteinIdType, "UNIPROT")
  checkEquals (x$pubmedID, "19111667")

  checkEquals (nrow (subset (tbl.md, is.na (proteinId) & !is.na (proteinIdType))), 0)

  invisible (tbl.md)

} # test.createMetadata
#-----------------------------------------------------------------------------------------------------------------------
test.getMatrixFilenames = function (dataDir)
{
  print ('--- test.getMatrixFilenames')

  checkEquals (length (getMatrixFilenames (dataDir)), 196)

} # test.getMatrixFilenames
#-----------------------------------------------------------------------------------------------------------------------
test.renameMatrices = function (dataDir, matrices, tbl.md, tbl.anno)
{
  print ('--- test.renameMatrices')

    # try it with just the first two matrices
  checkEquals (dim (tbl.md), c (196, 15))
  old.matrix.names = names (matrices)
  matrices.renamed = renameMatrices (matrices, tbl.md [1:2,])
  new.matrix.names = names (matrices.renamed)

  #print (old.matrix.names)
  #print (names (matrices.renamed))

  gene.names = sapply (strsplit (old.matrix.names, '\\.'), function (tokens) return (tokens [2]))
  author.names = sapply (strsplit (old.matrix.names, '\\.'), function (tokens) return (tokens [1]))

     # though order of gene and author is reversed, both should be found, in boh old and new matrix names
  for (m in 1:length (matrices)) {
    checkTrue (length (grep (gene.names [m], old.matrix.names [m])) == 1)
    checkTrue (length (grep (author.names [m], old.matrix.names [m])) == 1)
    checkTrue (length (grep (gene.names [m], new.matrix.names [m])) == 1)
    checkTrue (length (grep (author.names [m], new.matrix.names [m])) == 1)
    }


    # now check them all
  
  tbl.ref = createExperimentRefTable ()
  all.files = getMatrixFilenames (dataDir)
  matrices = readAndParse (file.path(dataDir, all.files))
  tbl.md = createMetadata (matrices, tbl.ref)

  old.matrix.names = names (matrices)
  matrices.renamed = renameMatrices (matrices, tbl.md)
  new.matrix.names = names (matrices.renamed)

  for (m in 1:length (matrices)) {
    checkTrue (length (grep (gene.names [m], old.matrix.names [m])) == 1)
    checkTrue (length (grep (author.names [m], old.matrix.names [m])) == 1)
    checkTrue (length (grep (gene.names [m], new.matrix.names [m])) == 1)
    checkTrue (length (grep (author.names [m], new.matrix.names [m])) == 1)
    }

  #printf ('validated %d new matrix names', length (matrices))

  invisible (matrices.renamed)

} # test.renameMatrices
#------------------------------------------------------------------------------------------------------------------------
test.createMatrixNameUniqifier = function ()
{
  print ('--- test.createMatrixNameUniqifier')
  
  data = c (8,8,8,7,2,2,5,6,5,3,10,5,8,6,9,5,8,8,4,5,8,6,9,2,1,0,5,7,7,2,4,4,3,7,7,9,9,6,1,3)
  test.matrix = matrix (data=data, nrow=4, ncol=10)
  uniqifier = createMatrixNameUniqifier (test.matrix)
  xxx <<- uniqifier
  checkEquals (uniqifier, "b42f")

} # test.createMatrixNameUniqifier
#------------------------------------------------------------------------------------------------------------------------
