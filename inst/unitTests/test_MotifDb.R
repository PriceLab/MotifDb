library (MotifDb)
library (RUnit)
library (MotIV)
library (seqLogo)
#------------------------------------------------------------------------------------------------------------------------
run.tests = function ()
{
  test.emptyCtor ()
  test.nonEmptyCtor ()
  test.MotifDb.normalMode ()
  test.MotifDb.emptyMode ()
  test.allMatricesAreNormalized ()
  test.providerNames ()
  test.geneSymbols ()
  test.geneIdsAndTypes ()
  test.proteinIds ()
  test.sequenceCount ()
  test.longNames ()
  test.organisms ()
  test.bindingDomains ()
  test.experimentTypes ()
  test.tfFamilies ()
  test.bindingSequences ()
  test.flyBindingDomains ()
  test.pubmedIDs ()
  test.allFullNames ()
  test.subset ()
  test.subsetWithVariables ()
  test.query ()
  test.transformMatrixToMemeRepresentation ()
  test.matrixToMemeText ()
  test.export_memeFormatStdOut ()
  test.export_memeFormatToFile ()
  test.export_memeFormatToFileDuplication ()
  test.export_memeFormatToFile_run_tomtom ()
  test.run_MotIV ()
  test.MotIV.toTable ()
  test.flyFactorGeneSymbols()

} # run.tests
#------------------------------------------------------------------------------------------------------------------------
test.emptyCtor = function ()
{
  print ('--- test.emptyCtor')
  motif.list = MotifDb:::MotifList ()
  checkEquals (length (motif.list), 0)  

} # test.emptyCtor
#------------------------------------------------------------------------------------------------------------------------
test.nonEmptyCtor = function ()
{
  print ('--- test.nonEmptyCtor')
  mtx = matrix (runif (20), nrow=4, ncol=5, byrow=T, dimnames=list(c ('A', 'C', 'G', 'T'), as.character (1:5)))
  mtx.normalized = apply (mtx, 2, function (colvector) colvector / sum (colvector))   
  matrixList = list (mtx.normalized)

  tbl.md = data.frame (providerName='',
                       providerId='',
                       dataSource='',
                       geneSymbol='',
                       geneId='',
                       geneIdType='',
                       proteinId='',
                       proteinIdType='',
                       organism='',
                       sequenceCount='',
                       bindingSequence='',
                       bindingDomain='',
                       tfFamily='',
                       experimentType='',
                       pubmedID='', 
                       stringsAsFactors=FALSE)
  names (matrixList) = 'test'
  rownames (tbl.md) = 'test'
  motif.list = MotifDb:::MotifList (matrixList, tbl.md)
  checkEquals (length (motif.list), 1)  

} # test.nonEmptyCtor
#------------------------------------------------------------------------------------------------------------------------
# 'normal' in that all included data sources are already loaded
test.MotifDb.normalMode = function ()
{
  print ('--- test.MotifDb.normalMode')

  mdb = MotifDb #  (quiet=TRUE)
    # (5 jun 2012)
    # JASPAR_CORE: 459
    # ScerTF: 196
    # UniPROBE: 380
    # FlyFactorSurvey: 614
    # hPDI: 437
  checkTrue (length (mdb) > 2080)

} # test.MotifDb.normalMode
#------------------------------------------------------------------------------------------------------------------------
# this mode is not intended for users, but may see use in the future.  
test.MotifDb.emptyMode = function ()
{
  print ('--- test.MotifDb.emptyMode')

  mdb = MotifDb:::.MotifDb (loadAllSources=FALSE, quiet=TRUE)
  checkTrue (length (mdb) == 0) 

} # test.MotifDb.emptyMode
#------------------------------------------------------------------------------------------------------------------------
# 3 jaspar matrices arrive with organism = NA.  find and fix.  make sure here.
# NA-JASPAR_CORE-TBP-MA0108.2: JASPAR gives <NA> for speciesID
# NA-JASPAR_CORE-HNF4A-MA0114.1: JASPAR gives <NA> for speciesID
# NA-JASPAR_CORE-CEBPA-MA0102.2: JASPAR gives '-' for speciesID, website says 'vertebrates'

test.noNAorganisms = function ()

{
  print ('--- test.noNAorganisms')
  checkEquals (which (is.na (mcols(MotifDb)$organism)), integer (0))

} # test.noNAorganisms
#------------------------------------------------------------------------------------------------------------------------
test.allMatricesAreNormalized = function ()
{
  print ('--- test.allMatricesAreNormalized')
  mdb = MotifDb# (quiet=TRUE)
  matrices = mdb@listData
    # a lenient test required by "Cparvum-UniPROBE-Cgd2_3490.UP00395" and  "Hsapiens-UniPROBE-Sox4.UP00401"
    # for reasons not yet explored.  10e-8 should be be possible
  checkTrue (all (sapply (matrices, function (m) all (abs (colSums (m) - 1.0) < 0.02))))
             
} # test.allMatricesAreNormalized
#------------------------------------------------------------------------------------------------------------------------
test.providerNames = function ()
{
  print ('--- test.getProviderNames')
  mdb = MotifDb # ()
  pn = mcols(mdb)$providerName
  checkEquals (length (which (is.na (pn))), 0)
  checkEquals (length (which (pn == '')), 0)

} # test.providerNames
#------------------------------------------------------------------------------------------------------------------------
test.geneSymbols = function ()
{
  print ('--- test.getGeneSymbols')
  mdb = MotifDb # ()
  syms = mcols(mdb)$geneSymbol
  checkEquals (length (which (is.na (syms))), 683)  # no symols yet for the dgf stamlab motifs
  checkEquals (length (which (syms == '')), 0)

} # test.geneSymbols
#------------------------------------------------------------------------------------------------------------------------
test.geneIdsAndTypes = function ()
{
  print ('--- test.getGeneIdsAndTypes')
  mdb = MotifDb
  geneIds = mcols(mdb)$geneId
  geneIdTypes = mcols(mdb)$geneIdType
  typeCounts = as.list (table (geneIdTypes))
  checkEquals(typeCounts, list(ENTREZ=2197, FLYBASE=30, SGD=453, comment=683))

  na.count = length (which (is.na (geneIds)))
  checkEquals (na.count, 304) 
  empty.count = length (which (geneIds == ''))
  checkEquals (empty.count, 0)


} # test.geneIdsAndTypes
#------------------------------------------------------------------------------------------------------------------------
# make sure that all proteinIds have explicit values, either proper identifiers or NA
# currently tested by looking for empty string assignments
test.proteinIds = function ()
{
  print ('--- test.proteinIds')
  mdb = MotifDb # (quiet=TRUE)
  NA.string.count = length (grep ('NA', mcols(mdb)$proteinId))
  checkEquals (NA.string.count, 0)
  
  empty.count = length (which (mcols(mdb)$proteinId==""))
  if (empty.count > 0)
    browser ('test.proteinIds')     

  checkEquals (empty.count, 0)

     # FlyFactorSurvey, as digested by me, had a blanket assigment of UNIPROT to all proteinIds
     # Herve' pointed out that this applied also to entries with no proteinId.
     # make sure this is fixed

  x = mcols(mdb)
  checkEquals (nrow (subset (x, !is.na (proteinIdType) & is.na (proteinId))), 0)

} # test.proteinIds
#------------------------------------------------------------------------------------------------------------------------
# only for UniPROBE do we not have sequence count.  might be possible to get them along with 'insertion sequences'
test.sequenceCount = function ()
{
  print ('--- test.sequenceCount')
  mdb = MotifDb # ()
  x = mcols(mdb)
  if (interactive ()) {
    x.up = subset (x, dataSource == 'UniPROBE')
    checkTrue (all (is.na (x.up$sequenceCount)))
    }
  else {
    uniprobe.indices = which (x$dataSource == 'UniPROBE')
    checkTrue (all (is.na (x$sequenceCount [uniprobe.indices])))
    }

} # test.sequenceCount
#------------------------------------------------------------------------------------------------------------------------
# make sure that a legitimate organism-dataSource-identifier is supplied for each matrix and as a rowname
# of the corresponding DataFrame
test.longNames = function ()
{
  print ('--- test.longNames')
  mdb = MotifDb 
  longNames = strsplit (names (mdb), '-')
  organisms = unique (sapply (longNames, '[', 1))
  
  dataSources = unique (sapply (longNames, '[', 2))

  recognized.dataSources = unique (mcols(mdb)$dataSource)
  recognized.organisms = unique (mcols(mdb)$organism)
    # a few (3) matrices from JASPAR core have NA organism.  make this into a character
    # so that it can be matched up against the 'NA' extracted from longNames just above
  na.indices = which (is.na (recognized.organisms))
  if (length (na.indices) > 0)
  recognized.organisms [na.indices] = 'NA'

  checkTrue (all (organisms %in% recognized.organisms))
  checkTrue (all (dataSources %in% recognized.dataSources))

} # test.longNames
#------------------------------------------------------------------------------------------------------------------------
# make sure that a legitimate organism is specified for each matrix
test.organisms = function ()
{
  print ('--- test.organisms')
  mdb = MotifDb # (quiet=TRUE)
  organisms = mcols(mdb)$organism

     # jaspar_core has 3 NA speciesId: TBP, HNF4A and CEBPA (MA0108.2, MA0114.1, MA0102.2)
     # their website shows these as vertebrates, which I map to 'Vertebrata'.  An organismID of '-'
     # gets the same treatment, matching website also.
  checkEquals (which (is.na (mcols(MotifDb)$organism)), integer (0))

  empty.count = length (which (mcols(mdb)$organism==""))
  checkEquals (empty.count, 0)

} # test.organisms
#------------------------------------------------------------------------------------------------------------------------
test.bindingDomains = function ()
{
  print ('--- test.bindingDomains')
  mdb = MotifDb # (quiet=TRUE)
  checkTrue (length (unique (mcols(mdb)$bindingDomain)) > 1)

} # test.bindingDomains
#------------------------------------------------------------------------------------------------------------------------
test.flyBindingDomains = function ()
{
  print ('--- test.flyBindingDomains')

  x = mcols(MotifDb)
  tmp = as.list (head (sort (table (subset (x, organism=='Dmelanogaster')$bindingDomain), decreasing=TRUE), n=3))

    # these counts will likely change with a fresh load of data from FlyFactorSurvey.

  checkEquals (tmp$Homeobox, 212)
  checkEquals (tmp[['zf-C2H2']], 160)
  checkEquals (tmp$HLH, 109)
  checkEquals (length (which (is.na (subset (x, organism=='Dmelanogaster')$bindingDomain))), 24)

} # test.flyBindingDomains
#------------------------------------------------------------------------------------------------------------------------
test.experimentTypes = function ()
{
  print ('--- test.experimentTypes')
  mdb = MotifDb # (quiet=TRUE)
  x = mcols(mdb)
  checkTrue (length (unique (x$experimentType)) >= 18)
  checkEquals (length (which (x$experimentType=='')), 0)

} # test.experimentTypes
#------------------------------------------------------------------------------------------------------------------------
test.tfFamilies = function ()
{
  print ('--- test.tfFamilies')
  mdb = MotifDb # (quiet=TRUE)
  checkTrue (length (unique (mcols(mdb)$tfFamily)) > 1)

} # test.tfFamilies
#------------------------------------------------------------------------------------------------------------------------
test.bindingSequences = function ()
{
  print ('--- test.bindingSequences')
  mdb = MotifDb # (quiet=TRUE)
  checkTrue (length (unique (mcols(mdb)$bindingSequence)) > 1)

} # test.tfFamilies
#------------------------------------------------------------------------------------------------------------------------
test.pubmedIDs = function ()
{
  print ('--- test.pubmedIDs')
  x = mcols(MotifDb) # (quiet=TRUE))
  checkTrue (length (unique (x$pubmedID)) >= 139)
  checkEquals (length (which (x$pubmedID == '')), 0)

} # test.pubmedIDs
#------------------------------------------------------------------------------------------------------------------------
# every matrix, and every corresponding metadata table row, has a name formed thus:
#
#    dataSource-organism-nativeGeneName
#
#  where nativeGeneName is how the matrix is named in the dataSource from which it comes.
#  thus:
#          FlyFactorSurvey-Dmelanogaster-Abd.A_FlyReg_FBgn0000014
#          UniPROBE-Scerevisiae-Asg1-UP00350
#          ScerTF-Scerevisiae-ABF2-badis
#          JASPAR_CORE-Rrattus-Ar-MA0007.1
# 
test.allFullNames = function ()
{
  print ('--- test.allFullNames')
  mdb = MotifDb # (quiet=TRUE)
  matrices = mdb@listData
  fullNames = names (matrices)

  all.dataSources = unique (mcols(mdb)$dataSource)
  checkTrue (length (all.dataSources) >= 4)
  
  for (source in all.dataSources) {
     this.dataSource <<- source
     matrices.by.source = subset (mdb, dataSource==this.dataSource)
     matrix.name = names (matrices.by.source)[1]
        #  FlyFactorSurvey: Dmelanogaster-FlyFactorSurvey-ab_SANGER_10_FBgn0259750
        #             hPDI: Hsapiens-hPDI-ABCF2
        #      JASPAR_CORE: JASPAR_CORE-Athaliana-AGL3-MA0001.1
        #         UniPROBE: Hsapiens-UniPROBE-Sox4.UP00401
     checkTrue (grep (this.dataSource, matrix.name) == 1)
     #printf ('%20s: %s', source, matrix.name)
     }

  return (TRUE)

} # test.allFullNames
#------------------------------------------------------------------------------------------------------------------------
test.subset = function ()
{
  if (interactive ()) {
    print ('--- test.subset')
    mdb = MotifDb
    checkTrue ('geneSymbol' %in% colnames (elementMetadata (mdb)))
    mdb.sub = subset (mdb, geneSymbol=='ABCF2')
    checkEquals (length (mdb.sub), 1)
    } # if interactive  
  
} # test.subset
#------------------------------------------------------------------------------------------------------------------------
test.subsetWithVariables = function ()
{
  if (interactive ()) {
    print ('--- test.subsetWithVariables')
  
    mdb = MotifDb # ()
    checkTrue ('geneSymbol' %in% colnames (elementMetadata (mdb)))
    target.gene <<- 'ABCF2'
    mdb.sub = subset (mdb, geneSymbol==target.gene)
    checkEquals (length (mdb.sub), 1)
    } # if interactive  
  
} # test.subsetWithVariables
#------------------------------------------------------------------------------------------------------------------------
test.query = function ()
{
  print ('--- test.query')
  mdb = MotifDb

    # do queries on dataSource counts match those from a contingency table?
  sources.list = as.list (table (mcols(mdb)$dataSource))
  checkEquals (length (query (mdb, 'flyfactorsurvey')), sources.list$FlyFactorSurvey)
  checkEquals (length (query (mdb, 'uniprobe')), sources.list$UniPROBE)
  checkEquals (length (query (mdb, 'UniPROBE')), sources.list$UniPROBE)

    # gene symbols which begin with 'sox' are quite common.  can we them?
    # there are currently (19 jul 2012) 18, but since this may change, our test is approximate

  sox.entries = query (mdb, '^sox')
  checkTrue (length (sox.entries) > 10)
  checkTrue (length (sox.entries) < 100)

    # manual inspection reveals that some of these genes have names which are all capitalized.  test that.
  checkTrue (length (query (mdb, '^sox', ignore.case=TRUE)) > length (query (mdb, '^SOX', ignore.case=FALSE)))

    # make sure that queries can be stacked up, and that order of call does not affect the outcome
   uniprobe.sox.matrices = query (query (mdb, 'uniprobe'), '^sox')
   sox.uniprobe.matrices = query (query (mdb, '^sox'), 'uniprobe')

   checkEquals (length (uniprobe.sox.matrices), length (sox.uniprobe.matrices))

     # an approximate count check
  checkTrue (length (uniprobe.sox.matrices) > 10)
  checkTrue (length (uniprobe.sox.matrices) < 30)

   checkEquals (unique (mcols(uniprobe.sox.matrices)$dataSource), 'UniPROBE')
   checkEquals (unique (mcols(sox.uniprobe.matrices)$dataSource), 'UniPROBE')
   gene.symbols = sort (unique (mcols(uniprobe.sox.matrices)$geneSymbol))
  
} # test.query
#------------------------------------------------------------------------------------------------------------------------
test.transformMatrixToMemeRepresentation = function ()
{
  print ('--- test.transformMatrixToMemeRepresentation')
  mdb = MotifDb # ()
  matrix.agl3 = subset (mdb, dataSource=='JASPAR_CORE' & organism=='Athaliana' & geneSymbol=='AGL3')[[1]]
  checkEquals (dim (matrix.agl3), c (4, 10))

    # a transposed frequency is needed for meme.  we store normalized frequency matrices, to just transposition is needed
  tfm = MotifDb:::transformMatrixToMemeRepresentation (matrix.agl3)
  checkEquals (dim (tfm), c (10, 4))
  checkTrue (all (round (rowSums (tfm)) == 1.0))

} # test.transformMatrixToMemeRepresentation
#------------------------------------------------------------------------------------------------------------------------
# translate a motif matrix from MotifList into text suitable for meme and tomtom input.
# choose the top two from UniPROBE.  they reliably have high counts, and easily distinguished normalized frequencies
test.matrixToMemeText = function ()
{
  print ('--- test.matrixToMemeText')
  sox4 = subset (MotifDb, dataSource=='UniPROBE' & organism=='Hsapiens' & geneSymbol=='Sox4')
  rtg3 = subset (MotifDb, dataSource=='UniPROBE' & organism=='Scerevisiae' & geneSymbol=='Rtg3')

  text.sox4 = MotifDb:::matrixToMemeText (sox4)
     # check these with t (sox4 [[1]])
  line1.sox4  = " 0.2457457130  0.1950426500  0.2287887620  0.3304228750"
  line14.sox4 = " 0.2821643030  0.2286132160  0.1585395830  0.3306828990"
  
  checkEquals (length (text.sox4), 29)
  checkEquals (text.sox4 [1], "MEME version 4")
  checkEquals (text.sox4 [10], "MOTIF Hsapiens-UniPROBE-Sox4.UP00401")
  checkEquals (grep (line1.sox4, text.sox4), 12)
  checkEquals (grep (line14.sox4, text.sox4), 25)

  text.rtg3 = MotifDb:::matrixToMemeText (rtg3)
  line1.rtg3  = " 0.3935122858  0.1453016447  0.3308830322  0.1303030373"
  line20.rtg3 = " 0.2490417648  0.3966478493  0.1083586569  0.2459517291"
  checkEquals (length (text.rtg3), 35)         # 4 trailing empty lines
  checkEquals (text.rtg3 [1], "MEME version 4")
  checkEquals (text.rtg3 [10], "MOTIF Scerevisiae-UniPROBE-Rtg3.UP00356")
  checkEquals (grep (line1.rtg3, text.rtg3), 12)
  checkEquals (grep (line20.rtg3, text.rtg3), 31)

    # now call with both matrices, and see if the right results are returned
  text.both = MotifDb:::matrixToMemeText (c (sox4, rtg3))
  checkEquals (text.both [1], "MEME version 4")
  checkEquals (grep (line1.sox4, text.both), 12)
  checkEquals (grep (line14.sox4, text.both), 25)
  checkEquals (grep (line1.rtg3, text.both),  29)
  checkEquals (grep (line20.rtg3, text.both), 48)

} # test.matrixToMemeText
#------------------------------------------------------------------------------------------------------------------------
# translate a motif matrix from MotifList into text suitable for meme and tomtom input.
# choose the top two from UniPROBE.  they reliably have high counts, and easily distinguished normalized frequencies
#test.matrixToMemeText_mapVersion = function ()
#{
#  print ('--- test.matrixToMemeText_mapVersion')
#  sox4 = subset (MotifDb, dataSource=='UniPROBE' & organism=='Hsapiens' & geneSymbol=='Sox4')
#  rtg3 = subset (MotifDb, dataSource=='UniPROBE' & organism=='Scerevisiae' & geneSymbol=='Rtg3')
#
#  text.sox4 = MotifDb:::mapped.broken.matrixToMemeText (sox4)
#  text.rtg3 = MotifDb:::mapped.broken.matrixToMemeText (rtg3)
#  text.both = MotifDb:::mapped.broken.matrixToMemeText (c (sox4, rtg3))
#
#     # check these with t (sox4 [[1]])
#  line1.sox4  = " 0.2457457130  0.1950426500  0.2287887620  0.3304228750"
#  line14.sox4 = " 0.2821643030  0.2286132160  0.1585395830  0.3306828990"
#
#  checkEquals (length (text.sox4), 29)
#  checkEquals (text.sox4 [1], "MEME version 4")
#  checkEquals (text.sox4 [10], "MOTIF Hsapiens-UniPROBE-Sox4.UP00401")
#  checkEquals (grep (line1.sox4, text.sox4), 12)
#  checkEquals (grep (line14.sox4, text.sox4), 25)
#
#  line1.rtg3  = " 0.3935122858  0.1453016447  0.3308830322  0.1303030373"
#  line20.rtg3 = " 0.2490417648  0.3966478493  0.1083586569  0.2459517291"
#  checkEquals (length (text.rtg3), 35)         # 4 trailing empty lines
#  checkEquals (text.rtg3 [1], "MEME version 4")
#  checkEquals (text.rtg3 [10], "MOTIF Scerevisiae-UniPROBE-Rtg3.UP00356")
#  checkEquals (grep (line1.rtg3, text.rtg3), 12)
#  checkEquals (grep (line20.rtg3, text.rtg3), 31)
#
#    # now call with both matrices, and see if the right results are returned
#  checkEquals (text.both [1], "MEME version 4")
#  checkEquals (grep (line1.sox4, text.both), 12)
#  checkEquals (grep (line14.sox4, text.both), 25)
#  checkEquals (grep (line1.rtg3, text.both),  29)
#  checkEquals (grep (line20.rtg3, text.both), 48)
#
#} # test.matrixToMemeText_mapVersion
#------------------------------------------------------------------------------------------------------------------------
test.export_memeFormatStdOut = function ()
{
  print ('--- test.export_memeFormatStdOut')
  mdb = MotifDb # ()
  mdb.chicken = subset (mdb, organism=='Gallus')
  checkEquals (length (mdb.chicken), 2)
    # text is cat-ed to stdout, so not avaialable here to check.
    # but just like print, export also returns the text invisibly.
    # so that CAN be checked.
  
  meme.text = export (mdb.chicken, format='meme')
  checkEquals (length (meme.text), 1)   # just one long string
  checkTrue (is.character (meme.text))
  checkTrue (nchar (meme.text) > 800)   # 1002 as of (10 aug 2012)
  return (TRUE)

} # test.exportMemeFormatToStdOut
#------------------------------------------------------------------------------------------------------------------------
test.export_memeFormatToFile = function ()
{
  print ('--- test.export_memeFormatToFile')
  mdb = MotifDb # ()
  mdb.chicken = subset (mdb, organism=='Gallus')
  checkEquals (length (mdb.chicken), 2)
  output.file = tempfile ()
  meme.text = export (mdb.chicken, output.file, 'meme')
  retrieved = scan (output.file, what=character (0), sep='\n', quiet=TRUE)
  invisible (retrieved)

} # test.exportMemeFormatToFile
#------------------------------------------------------------------------------------------------------------------------
test.export_memeFormatToFileDuplication = function ()
{
  print ('--- test.export_memeFormatToFileDuplication')
  mdb = MotifDb # ()
  mdb.mouse = subset (mdb, organism=='Mmusculus')
  checkEquals (length (mdb.mouse), 462)
  output.file = 'mouse.txt' # tempfile ()
  max = 3
  meme.text = export (mdb.mouse [1:max], output.file, 'meme')
  retrieved = scan (output.file, what=character (0), sep='\n', quiet=TRUE)
  invisible (retrieved)

} # test.exportMemeFormatToFileDuplication
#------------------------------------------------------------------------------------------------------------------------
test.export_memeFormatToFile_run_tomtom = function (max=50)
{
  if (interactive ()) {
    print ('--- test.export_memeFormatToFile_run_tomtom')
    mdb = MotifDb
    sox4.mouse = subset (mdb, organism=='Mmusculus' & geneSymbol=='Sox4')
    all.human = subset (mdb, organism=='Hsapiens')

    tomtom.tmp.dir = tempfile ()
    print (tomtom.tmp.dir)
    stopifnot (dir.create (tomtom.tmp.dir))
    sox4.file.path = file.path (tomtom.tmp.dir, 'sox4.mouse.text')
    all.human.file.path = file.path (tomtom.tmp.dir,'all.human.text')
    export (sox4.mouse, sox4.file.path, 'meme')
    export (all.human, all.human.file.path, 'meme')

       # find similarity of motif #1 to all the motifs in mdbMany

    cmd = sprintf ('tomtom -no-ssc -oc %s -verbosity 3 -min-overlap 5 -mi 1 -dist pearson -evalue -thresh 10 %s %s',
                    tomtom.tmp.dir, sox4.file.path, all.human.file.path)
    system (cmd)
    cmd = sprintf ('open %s/tomtom.html', tomtom.tmp.dir)
    system (cmd)
    } # if interactive

} # test.export_memeFormatToFile_run_tomtom
#------------------------------------------------------------------------------------------------------------------------
test.run_MotIV = function ()
{
  library (MotIV)
  print ('--- test.run_MotIV')
  mdb = MotifDb # ()

  db.tmp = mdb@listData

     # match motif 1 against the entire MotifDb  collection
  motif.hits = motifMatch (db.tmp [1], database=db.tmp)
     # the long way to extract the matrix name.  see MotIV.toTable below for more convenient way
  checkEquals (motif.hits@bestMatch[[1]]@aligns[[1]]@TF@name, names (db.tmp)[1])

     # jaspar, uniprobe and ScerTF contribute a total of 1035 matrices
  checkTrue (length (mdb) >= 1035)

  motif.hits =  motifMatch (db.tmp [1035], database=db.tmp)
  tbl.hits = MotIV.toTable (motif.hits)
  checkEquals (names (db.tmp)[1035], tbl.hits [1, 'name'])   # first hit should be the target itself
  invisible (tbl.hits)
  
} # test.run_MotIV
#------------------------------------------------------------------------------------------------------------------------
MotIV.toTable = function (match)
{
  stopifnot (length (match@bestMatch) >= 1)
  alignments = match@bestMatch[[1]]@aligns

  df = data.frame (stringsAsFactors=FALSE)
  for (alignment in alignments) {
    x = alignment
    name = x@TF@name
    eVal = x@evalue
    sequence = x@sequence
    match = x@match
    strand = x@strand
    df = rbind (df, data.frame (name=name, eVal=eVal, sequence=sequence, match=match, strand=strand, stringsAsFactors=FALSE))
    } # for alignment

  return (df)

} # MotIV.toTable 
#------------------------------------------------------------------------------------------------------------------------
test.MotIV.toTable = function ()
{
  print ('--- test.MotIVtoTable')
  mdb = MotifDb # ()
  test.hits = motifMatch (mdb[1]@listData, database=jaspar)
  tbl.hits =  MotIV.toTable (test.hits)
  checkEquals (dim (tbl.hits), c (5, 5))
  checkEquals (colnames (tbl.hits), c ("name", "eVal", "sequence", "match", "strand"))

} # test.MotIV.toTable 
#------------------------------------------------------------------------------------------------------------------------
pwmMatch.toTable = function (motifMatch) {
   if (length (motifMatch@bestMatch) == 0)
   return (NA)

   df.list = vector("list", length(motifMatch@bestMatch))
   for (k in seq(length(motifMatch@bestMatch)))
   {
       alignments = motifMatch@bestMatch[[k]]@aligns
       df = data.frame (stringsAsFactors=FALSE)
       for (alignment in alignments) {
           x = alignment
           name = x@TF@name
           eVal = x@evalue
           sequence = x@sequence
           match = x@match
           strand = x@strand
           df = rbind (df, data.frame (name=name, eVal=eVal, sequence=sequence, match=match, strand=strand, stringsAsFactors=FALSE))
       } # for alignment
       df.list[[k]]=df
   }
   names(df.list) <- names(motifMatch)
   return (df.list)

} # pwmMatch.toTable
#------------------------------------------------------------------------------
# Robert Stojnic reports incorrect gene symbols for matrices obtained from
# flyFactorSurvey.
# the solution was to abandon the original strategy of extracting the
# symbol from the matrix (and file) name.
# now, the flybase importer ("inst/scripts/import/flyFactorSurvey/import.R")
# uses FBgn id (which can be reliably extracted) and uses indpendent
# data sources to learn the gene symbol.
#
# robert's email:
#  I'm working on using MotifDb motifs in my PWMEnrich package and I
#  have noticed that there is a slight problem with gene symbols for
#  Drosophila. In particular, the gene symbols do not always correspond
#  to the gene ID and are frequently mis-capitalized. In Drosophila z
#  and Z are two different genes and capitalization does matter if
#  someone is to use the gene symbols. Also, in some cases the symbols
#  are missing hyphens or parenthesis. I have used the gene IDs and the
#  Flybase annotation database to set the correct gene symbols for
#  Drosophila, please find attached the result of my re-annotation.
#
#  looking at his correctedMotifDbDmel.csv 
#
#    head(read.table("correctedMotifDbDmel.csv", sep=",", header=TRUE, stringsAsFactors=FALSE))
#                  providerName oldGeneSymbol newGeneSymbol
#    1 ab_SANGER_10_FBgn0259750            Ab            ab
#    2  ab_SOLEXA_5_FBgn0259750            Ab            ab
#    3 Abd-A_FlyReg_FBgn0000014         Abd-a         abd-A
#    4 Abd-B_FlyReg_FBgn0000015         Abd-b         Abd-B
#    5    AbdA_Cell_FBgn0000014          Abda         abd-A
#    6  AbdA_SOLEXA_FBgn0000014          Abda         abd-A
#
test.flyFactorGeneSymbols <- function()
{
    print ("--- test.flyFactorGeneSymbols")
    mdb = MotifDb
    checkEquals(mcols(query(mdb, "FBgn0259750"))$geneSymbol, c("ab", "ab"))
    checkEquals(mcols(query(mdb, "FBgn0000014"))$geneSymbol, rep("abd-A", 3))
    checkEquals(mcols(query(mdb, "FBgn0000015"))$geneSymbol, rep("Abd-B", 3))

} # test.flyFactorGeneSymbols
#-------------------------------------------------------------------------------

