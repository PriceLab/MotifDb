# flyFactorSurvey/test.R
#-------------------------------------------------------------------------------
library (RUnit)
#-------------------------------------------------------------------------------
source("import.R")
#-------------------------------------------------------------------------------
run.tests = function (dataDir)
{
    dataDir <- file.path(dataDir, "flyFactorSurvey")
    freshStart ()
    test.fbgnToIDs()
    x.list.BD <- test.createXrefBindingDomain (dataDir)
    x.xref <- test.createXref ()
    x.filenames <- test.getMatrixFilenames (dataDir)
    x.tbl.ref <- test.createExperimentRefTable ()
    x.pwm <- test.parsePWMfromText (dataDir)
    x.mat3 <- test.readAndParse (dataDir)
    x.tbl.md <- test.createMetadata (x.mat3, x.tbl.ref, x.xref, x.list.BD)
    x.mat3f <- test.normalizeMatrices (x.mat3)
    x.mat3fr <- test.renameMatrices (x.mat3f, x.tbl.md [1:3,])

} # run.tests
#-------------------------------------------------------------------------------
freshStart = function ()
{
    output.files.easy.to.regenerate = grep ('RData$', dir (), value=TRUE)
    if (length (output.files.easy.to.regenerate) > 0)
        unlink (output.files.easy.to.regenerate)

} # freshStart
#-------------------------------------------------------------------------------
test.createXrefBindingDomain = function (dataDir)
{
    print ('--- test.createXrefBindingDomain')
    list.BD = createXrefBindingDomain (dataDir)
    checkEquals (length (list.BD), 777)
    checkTrue (all (grepl ('FBgn', names (list.BD))))
        # most abundant domain is zf-C2H2
    checkEquals(names (sort (table (as.character (list.BD)), decreasing=TRUE)[1]),
                "zf-C2H2")
    invisible(list.BD)

} # test.createXrefBindingDomain
#-------------------------------------------------------------------------------
test.createXref = function ()
{
    print ('--- test.createXref')
    xref = createXref ()
    checkEquals (xref [['FBgn0261819']], 'F0J881')
    #checkEquals (xref [['FBgn0259750']], 'E1JHF4') # gone with org.Dm.eg.db_2.9.0
    checkEquals (xref [['FBgn0000014']], 'E1JIM9')
    invisible (xref)

} # test.createXref
#-------------------------------------------------------------------------------
test.getMatrixFilenames = function (dataDir)
{
    print ('--- test.getMatrixFilenames')

    filenames = getMatrixFilenames (dataDir)
       # as of (03 apr 2013): 614 matrix files, one binding domains extra
       # file, TFfile2b.tsv
    
    checkTrue(length (filenames) > 600)
    sample.files <- c("Abd-A_FlyReg_FBgn0000014.pfm", "Abd-B_FlyReg_FBgn0000015.pfm",
                      "AbdA_Cell_FBgn0000014.pfm")
    checkTrue(all(sample.files %in% filenames))

    invisible (filenames)

} # test.getMatrixFilenames
#-------------------------------------------------------------------------------
test.createExperimentRefTable = function ()
{
    print ('--- test.createExperimentRefTable')
    tbl.ref = createExperimentRefTable ()
    checkEquals (nrow (tbl.ref), 6)

    invisible (tbl.ref)

} # test.createExperimentRefTable
#-------------------------------------------------------------------------------
test.parsePWMfromText = function (dataDir)
{
    print ('--- test.parsePWMfromText')
    path <- file.path(dataDir, 'tup_SOLEXA_10_FBgn0003896.pfm')
    checkTrue(file.exists(path))
    
    lines.of.text =  scan (path, sep='\n',
                           what=character(0), comment='#', quiet=TRUE)
    pwm.tup =  parsePWMfromText (lines.of.text)
    checkEquals (dim (pwm.tup), c (4, 9))
    checkEquals (colnames (pwm.tup), as.character (1:9))
    checkEquals (rownames (pwm.tup), c ('A', 'C', 'G', 'T'))

  invisible (pwm.tup)
  
} # test.parsePWMfromText
#-------------------------------------------------------------------------------
test.readAndParse = function (dataDir)
{
    print ('--- test.readAndParse')

    all.files = file.path(dataDir,
                          getMatrixFilenames (dataDir))
    
    sample.files <- c("Abd-A_FlyReg_FBgn0000014.pfm",
                      "Abd-B_FlyReg_FBgn0000015.pfm",
                      "AbdA_Cell_FBgn0000014.pfm")

    checkTrue(all(sample.files %in% list.files(dataDir)))
    sample.files <- file.path(dataDir, sample.files)
    checkTrue(all(file.exists(sample.files)))

      # by sorting these filenames, we know the order of the three returned
      # matrices, and results can easily be checked
    
    mtx.test = readAndParse (sample.files)

    checkEquals (length (mtx.test), 3)
        # names are alphabetized, and sort differently in linux and macos
        # work around that by 
    expected.names = c("Abd-A_FlyReg_FBgn0000014.pfm",
                       "Abd-B_FlyReg_FBgn0000015.pfm",
                       "AbdA_Cell_FBgn0000014.pfm")
    actual.names = sort(names(mtx.test))
    checkEquals(expected.names, actual.names)

    checkEquals (dim (mtx.test [[1]]), c (4, 8))
    checkEquals (dim (mtx.test [[2]]), c (4, 8))
    checkEquals (dim (mtx.test [[3]]), c (4, 7))

    invisible (mtx.test)

} # test.readAndParse
#-------------------------------------------------------------------------------
test.createMetadata = function (matrices, tbl.ref, xref, list.BD)
{
    print ('--- test.createMetadata')

    tbl.md = createMetadata (matrices, tbl.ref, xref, list.BD)
    checkEquals (nrow (tbl.md), length (matrices))
    checkEquals (ncol (tbl.md), 15)
    expected <- c("providerName", "providerId", "dataSource", "geneSymbol",
                  "geneId", "geneIdType", "proteinId", "proteinIdType",
                  "organism", "sequenceCount", "bindingSequence",
                  "bindingDomain", "tfFamily", "experimentType", "pubmedID")

    checkEquals (colnames (tbl.md), expected)
    tester =  "Dmelanogaster-FlyFactorSurvey-Abd.A_FlyReg_FBgn0000014"
    checkTrue (tester %in% rownames (tbl.md))
    x = as.list (tbl.md [tester,])
    checkEquals (x$providerName, "Abd-A_FlyReg_FBgn0000014")
    checkEquals (x$providerId, "FBgn0000014")
    checkEquals (x$geneSymbol, "abd-A")
    checkEquals (x$geneId, "42037")
    checkEquals (x$geneIdType, "ENTREZ")
    checkEquals (x$proteinId, "E1JIM9")
    checkEquals (x$proteinIdType, "UNIPROT")
    checkEquals (x$organism, "Dmelanogaster")
    checkEquals (x$sequenceCount, 37)
    checkTrue   (is.na (x$bindingSequence))
    checkEquals (x$bindingDomain, 'Homeobox')
    checkTrue   (is.na (x$tfFamily))
    checkEquals (x$experimentType, "DNASE I footprinting")
    checkTrue   (is.na (x$pubmedID))
  
    invisible (tbl.md)

} # test.createMetadata
#-------------------------------------------------------------------------------
test.normalizeMatrices = function (matrices)
{
  print ('--- test.normalizeMatrices')

  colsums = as.integer (sapply (matrices, function (mtx) as.integer (mean (round (colSums (mtx))))))
  checkTrue (all (colsums > 1))

  matrices.norm = normalizeMatrices (matrices)

  colsums = as.integer (sapply (matrices.norm, function (mtx) as.integer (mean (round (colSums (mtx))))))
  checkTrue (all (colsums == 1))
  invisible (matrices.norm)

} # test.normalizeMatrices
#------------------------------------------------------------------------------------------------------------------------
test.renameMatrices = function (matrices, tbl.md)
{
    print ('--- test.renameMatrices')
    stopifnot (length (matrices) == 3)
    stopifnot (nrow (tbl.md) == 3)

    old.names = names (matrices)
    mtxr = renameMatrices (matrices, tbl.md)

      # make sure the dash in the incoming, eg,
      #     "Abd-A_FlyReg_FBgn0000014.pfm"
      # is converted to a '.', reserving dashes for 
      # the separation of dataSource-species-nativeGeneName


    checkEquals (names (mtxr) [1],
                 "Dmelanogaster-FlyFactorSurvey-Abd.A_FlyReg_FBgn0000014")

    checkEquals (names (mtxr) [2],
                 "Dmelanogaster-FlyFactorSurvey-Abd.B_FlyReg_FBgn0000015")

    checkEquals (names (mtxr) [3],
                 "Dmelanogaster-FlyFactorSurvey-AbdA_Cell_FBgn0000014")

    invisible (mtxr)

} # renameMatrices
#-------------------------------------------------------------------------------
test.fbgnToIDs <- function()
{
    print("--- test.fbgnToIDs")
       # first, make sure that unmappable FBgn ids are handled properly
    bogus <- fbgnToIDs("bogus")
    checkEquals(dim(bogus), c(1,4))
    checkEquals(rownames(bogus), "bogus")
    checkEquals(as.character(as.list(bogus[1,],)), c(rep("bogus", 3), "FLYBASE"))

       # disable the translation of NA to rowname
    bogus <- fbgnToIDs("bogus", useInputForMissingValues=FALSE)
    checkEquals(dim(bogus), c(1,4))
    checkEquals(rownames(bogus), "bogus")
    checkTrue(all(is.na(bogus[1,1:3])))

    fbgns <- c ("FBgn0003267", "FBgn0000611", "FBgn0004394", "FBgn0027364",
                "FBgn0000413")
    tbl.ids <- fbgnToIDs(fbgns)
    checkEquals(tbl.ids$flybase_id, fbgns)
    checkEquals(tbl.ids$symbol, c("ro", "exd", "pdm2", "Six4", "da"))
    checkEquals(tbl.ids$geneIdType, rep("ENTREZ", 5))

         # 0259750 is obsolete, replace by 0264442
    tbl.ids <- fbgnToIDs("FBgn0259750")
    checkEquals(as.list(tbl.ids),
                list(flybase_id="FBgn0264442", symbol="ab", gene_id="34560",
                     geneIdType="ENTREZ"))

       # http://flybase.org/static_pages/downloads/IDConv.html
       # FBgn0003986 	 FBgn0261930 	FBgn0261930	 vnd 
    tbl.ids <- fbgnToIDs("FBgn0003986")
    checkEquals(as.character(as.list(tbl.ids[1,],)),
                c(rep("FBgn0003986", 3), "FLYBASE"))
    
    serializedFbgnsForTesting <- system.file(package="MotifDb", "scripts",
                                             "import", "fbgns.RData")
    checkTrue(file.exists(serializedFbgnsForTesting))
    load(serializedFbgnsForTesting)
    checkEquals(length(fbgns), 326)
    tbl.ids <- fbgnToIDs(fbgns)
    checkEquals(dim(tbl.ids), c(length(fbgns), 4))

       # no NA's should survive
    checkTrue(!any(is.na(tbl.ids)))

       # what percentage of the 326 ids fail to map?  do NOT convert NAs
    tbl.ids <- fbgnToIDs(fbgns, useInputForMissingValues=FALSE)
    failure.count <- length(which(is.na(tbl.ids$flybase_id)))
    checkEquals(failure.count, 18)

} # test.fbgnToIDsp
#-------------------------------------------------------------------------------
# robert stojnic reported that MotifDb 1.1.6 (and before) had some incorrect
# gene symbols, as seen here:
#   library(MotifDb)
#   mdb <- MotifDb
#   x <- query(mdb, "ab_SANGER_10_FBgn0259750")
#   t(as.matrix(mcols(x)))
#                 [,1]                                   
#   providerName    "ab_SANGER_10_FBgn0259750"             
#   providerId      "FBgn0259750"
#   dataSource      "FlyFactorSurvey"                      
#   geneSymbol      "Ab"                                   
#   geneId          "FBgn0259750"                          
#   geneIdType      "FLYBASE"                              
#   proteinId       "E1JHF4"                               
#   proteinIdType   "UNIPROT"                              
#   organism        "Dmelanogaster"                        
#   sequenceCount   "20"                                   
#   bindingSequence NA                                     
#   bindingDomain   "zf-C2H2"                              
#   tfFamily        NA                                     
#   experimentType  "bacterial 1-hybrid, SANGER sequencing"
#   pubmedID        NA
#
# FBgn0259750 is obsolete.  flybase provides a converter:
#    http://flybase.org/static_pages/downloads/IDConv.html shows this, and provide the current FBgn:
#
#           Submitted ID          Current ID    Converted ID    Related record 
#            FBgn0259750          FBgn0264442    FBgn0264442                ab 
#
# make sure we get this right
# the matrix files are:
#
#     dataDir/ab_SANGER_10_FBgn0259750.pfm
#     dataDir/ab_SOLEXA_5_FBgn0259750.pfm
#
test.geneNaming <- function(dataDir)
{

   dataDir <- "/Users/pshannon/s/data/public/TFBS/flyFactorSurvey/unpacked"
   filenames = getMatrixFilenames (dataDir)

   removers <- grep (bindingDomainXrefSourceFile(), filenames)
   if(length(removers) > 0)
       filenames <- filenames[-removers]
  
   keepers <- grep("FBgn0259750", filenames)
   checkTrue(length(keepers) == 2)
   filenames <- filenames[keepers]
   full.filenames = file.path(dataDir, filenames)
   list.BD = createXrefBindingDomain (dataDir)
   xref = createXref ()  # maps flybase ids to uniprot
   tbl.ref = createExperimentRefTable ()
   matrices = readAndParse (full.filenames)
   tbl.md = createMetadata (matrices, tbl.ref, xref, list.BD)
     # now normalize, not in readAndParse, because we need the counts to go into tbl.md$sequenceCount
   checkEquals(tbl.md$geneSymbol, c("ab", "ab"))
   matrices = normalizeMatrices (matrices)
   matrices = renameMatrices (matrices, tbl.md)

} # test.geneNaming
#-------------------------------------------------------------------------------
