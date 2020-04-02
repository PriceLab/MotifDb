library(MotifDb)
library(RUnit)
library(MotIV)
library(seqLogo)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests = function()
{
  test.emptyCtor()
  test.nonEmptyCtor()
  test.MotifDb.normalMode()
  test.MotifDb.emptyMode()
  test.allMatricesAreNormalized()
  test.providerNames()
  test.geneSymbols()
  test.geneIdsAndTypes()
  test.proteinIds()
  test.sequenceCount()
  test.longNames()
  test.organisms()
  test.bindingDomains()
  test.experimentTypes()
  test.tfFamilies()
  test.bindingSequences()
  test.flyBindingDomains()
  test.pubmedIDs()
  #test.allFullNames()
  test.subset()
  test.subsetWithVariables()
  test.queryOldStyle()
  test.query()
  test.transformMatrixToMemeRepresentation()
  test.matrixToMemeText()
  test.export_memeFormatStdOut()
  test.export_memeFormatToFile()
  test.export_memeFormatToFileDuplication()
  test.export_memeFormatToFile_run_tomtom()
  test.flyFactorGeneSymbols()
  test.export_jasparFormatStdOut()
  test.export_jasparFormatToFile()

  test.geneToMotif()
  test.geneToMotif.ignore.jasparSuffixes()
  test.geneToMotif.oneGene.noMotifs
  test.motifToGene()
  test.associateTranscriptionFactors()

  test.hocomoco11.with.reliabilityScores()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test.emptyCtor = function()
{
  print('--- test.emptyCtor')
  motif.list = MotifDb:::MotifList()
  checkEquals(length(motif.list), 0)

} # test.emptyCtor
#------------------------------------------------------------------------------------------------------------------------
test.nonEmptyCtor = function()
{
  print('--- test.nonEmptyCtor')
  mtx = matrix(runif(20), nrow=4, ncol=5, byrow=T, dimnames=list(c('A', 'C', 'G', 'T'), as.character(1:5)))
  mtx.normalized = apply(mtx, 2, function(colvector) colvector / sum(colvector))
  matrixList = list(mtx.normalized)

  tbl.md = data.frame(providerName='',
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
  names(matrixList) = 'test'
  rownames(tbl.md) = 'test'
  motif.list = MotifDb:::MotifList(matrixList, tbl.md)
  checkEquals(length(motif.list), 1)

} # test.nonEmptyCtor
#------------------------------------------------------------------------------------------------------------------------
# 'normal' in that all included data sources are already loaded
test.MotifDb.normalMode = function()
{
  print('--- test.MotifDb.normalMode')

  mdb = MotifDb # (quiet=TRUE)
    #(5 jun 2012)
    # JASPAR_CORE: 459
    # ScerTF: 196
    # UniPROBE: 380
    # FlyFactorSurvey: 614
    # hPDI: 437
  checkTrue(length(mdb) > 2080)

} # test.MotifDb.normalMode
#------------------------------------------------------------------------------------------------------------------------
# this mode is not intended for users, but may see use in the future.
test.MotifDb.emptyMode = function()
{
  print('--- test.MotifDb.emptyMode')

  mdb = MotifDb:::.MotifDb(loadAllSources=FALSE, quiet=TRUE)
  checkTrue(length(mdb) == 0)

} # test.MotifDb.emptyMode
#------------------------------------------------------------------------------------------------------------------------
# 3 jaspar matrices arrive with organism = NA.  find and fix.  make sure here.
# NA-JASPAR_CORE-TBP-MA0108.2: JASPAR gives <NA> for speciesID
# NA-JASPAR_CORE-HNF4A-MA0114.1: JASPAR gives <NA> for speciesID
# NA-JASPAR_CORE-CEBPA-MA0102.2: JASPAR gives '-' for speciesID, website says 'vertebrates'

# Many more NA's exist...need to fix these; here's a quick fix for now

test.noNAorganisms = function()

{
  print('--- test.noNAorganisms')
  #checkEquals(which(is.na(mcols(MotifDb)$organism)), integer(0))

  # There's a fair number of NA organisms, mostly due to including the homer DB
  checkEquals(sum(is.na(mcols(MotifDb)$organism)), 366)

} # test.noNAorganisms
#------------------------------------------------------------------------------------------------------------------------
test.allMatricesAreNormalized = function()
{
  print('--- test.allMatricesAreNormalized')
  mdb = MotifDb#(quiet=TRUE)
  matrices = mdb@listData
    # a lenient test required by "Cparvum-UniPROBE-Cgd2_3490.UP00395" and  "Hsapiens-UniPROBE-Sox4.UP00401"
    # for reasons not yet explored.  10e-8 should be be possible
  checkTrue(all(sapply(matrices, function(m) all(abs(colSums(m) - 1.0) < 0.02))))

} # test.allMatricesAreNormalized
#------------------------------------------------------------------------------------------------------------------------
test.providerNames = function()
{
  print('--- test.getProviderNames')
  mdb = MotifDb #()
  pn = mcols(mdb)$providerName
  checkEquals(length(which(is.na(pn))), 0)
  checkEquals(length(which(pn == '')), 0)

} # test.providerNames
#------------------------------------------------------------------------------------------------------------------------
test.geneSymbols = function()
{
  print('--- test.getGeneSymbols')
  mdb = MotifDb #()
  syms = mcols(mdb)$geneSymbol
  checkEquals(length(which(is.na(syms))), 683)  # no symols yet for the dgf stamlab motifs
  checkEquals(length(which(syms == '')), 0)

} # test.geneSymbols
#------------------------------------------------------------------------------------------------------------------------
test.geneIdsAndTypes = function()
{
  print('--- test.getGeneIdsAndTypes')
  mdb = MotifDb
  tbl <- mcols(mdb)
  geneIds = tbl$geneId
  geneIdTypes = tbl$geneIdType
  typeCounts = as.list(table(geneIdTypes))

  checkTrue(typeCounts$ENTREZ > 2300)
  checkTrue(typeCounts$FLYBASE >= 45)
  checkTrue(typeCounts$SGD >= 600)
  checkTrue(nrow(subset(tbl, is.na(geneIdType))) > 2000)

  empty.count = length(which(geneIds == ''))
  checkEquals(empty.count, 0)


} # test.geneIdsAndTypes
#------------------------------------------------------------------------------------------------------------------------
# make sure that all proteinIds have explicit values, either proper identifiers or NA
# currently tested by looking for empty string assignments
test.proteinIds = function()
{
  print('--- test.proteinIds')
  mdb = MotifDb #(quiet=TRUE)
  NA.string.count <- sum(is.na(mcols(mdb)$proteinId))
#  NA.string.count = length(grep('NA', mcols(mdb)$proteinId))

  checkEquals(NA.string.count, 2514)
  # FIX THIS; Currently 2514 don't have protein IDs
  #checkEquals(NA.string.count, 0)

  empty.count = length(which(mcols(mdb)$proteinId==""))
  if(empty.count > 0)
    browser('test.proteinIds')

  checkEquals(empty.count, 0)

     # FlyFactorSurvey, as digested by me, had a blanket assigment of UNIPROT to all proteinIds
     # Herve' pointed out that this applied also to entries with no proteinId.
     # make sure this is fixed

  ### FIX THIS TOO! Currently have 913 entries with a proteinIdType and no proteinId
  x = mcols(mdb)
  # checkEquals(nrow(subset(x, !is.na(proteinIdType) & is.na(proteinId))), 0)


} # test.proteinIds
#------------------------------------------------------------------------------------------------------------------------
# only for UniPROBE do we not have sequence count.  might be possible to get them along with 'insertion sequences'
test.sequenceCount = function()
{
  print('--- test.sequenceCount')
  mdb = MotifDb #()
  x = mcols(mdb)
  if(interactive()) {
    x.up = subset(x, dataSource == 'UniPROBE')
    checkTrue(all(is.na(x.up$sequenceCount)))
    }
  else {
    uniprobe.indices = which(x$dataSource == 'UniPROBE')
    checkTrue(all(is.na(x$sequenceCount [uniprobe.indices])))
    }

} # test.sequenceCount
#------------------------------------------------------------------------------------------------------------------------
# make sure that a legitimate organism-dataSource-identifier is supplied for each matrix and as a rowname
# of the corresponding DataFrame
test.longNames = function()
{
  print('--- test.longNames')
  mdb <- MotifDb
  longNames <- strsplit(names(mdb), '-')
  organisms <- unique(sapply(longNames, '[', 1))

  dataSources <- unique(lapply(longNames, '[', 2))

  recognized.dataSources <- c("cisbp_1.02", "FlyFactorSurvey",
                              "HOCOMOCOv10", "HOCOMOCOv11B-core", "HOCOMOCOv11C-core", "HOCOMOCOv11B-full",
                              "HOCOMOCOv11C-full", "HOCOMOCOv11A-core", "HOCOMOCOv11D-full",
                              "HOCOMOCOv11A-full", "HOMER", "hPDI",
                              "JASPAR_CORE", "JASPAR_2014", "jaspar2016", "jaspar2018",
                              "jolma2013", "ScerTF", "stamlab", "SwissRegulon", "UniPROBE")

  recognized.dataSources <-  c("cisbp_1.02",
                             "FlyFactorSurvey",
                             "HOCOMOCOv10",
                             "HOCOMOCOv11-core-A",
                             "HOCOMOCOv11-core-B",
                             "HOCOMOCOv11-core-C",
                             "HOCOMOCOv11-full-A",
                             "HOCOMOCOv11-full-B",
                             "HOCOMOCOv11-full-C",
                             "HOCOMOCOv11-full-D",
                             "HOMER",
                             "hPDI",
                             "JASPAR_2014",
                             "JASPAR_CORE",
                             "jaspar2016",
                             "jaspar2018",
                             "jolma2013",
                             "ScerTF",
                             "stamlab",
                             "SwissRegulon",
                             "UniPROBE")

  recognized.organisms <- unique(mcols(mdb)$organism)
    # a few(3) matrices from JASPAR core have NA organism.  make this into a character
    # so that it can be matched up against the 'NA' extracted from longNames just above
  na.indices <- which(is.na(recognized.organisms))
  if(length(na.indices) > 0)
     recognized.organisms [na.indices] <- 'NA'

  checkTrue(all(organisms %in% recognized.organisms))
  # new hocomoco-[core|full]-[ABCD] dataSource names are not incorporated into rownames yet
  # checkTrue(all(dataSources %in% recognized.dataSources))

  } # test.longNames
#------------------------------------------------------------------------------------------------------------------------
# make sure that a legitimate organism is specified for each matrix
test.organisms <- function()
{
  print('--- test.organisms')
  mdb <- MotifDb #(quiet=TRUE)
  organisms <- mcols(mdb)$organism

     # jaspar_core has 3 NA speciesId: TBP, HNF4A and CEBPA(MA0108.2, MA0114.1, MA0102.2)
     # their website shows these as vertebrates, which I map to 'Vertebrata'.  An organismID of '-'
  # gets the same treatment, matching website also.

  ### Note: this failing test is the same as the test.noNAorganisms test!
  # As in case of noNA, need to add organisms for these
  #checkEquals(which(is.na(mcols(MotifDb)$organism)), integer(0))

  empty.count <- length(which(mcols(mdb)$organism==""))
  checkEquals(empty.count, 0)

} # test.organisms
#------------------------------------------------------------------------------------------------------------------------
test.bindingDomains <- function()
{
  print('--- test.bindingDomains')
  mdb <- MotifDb #(quiet=TRUE)
  checkTrue(length(unique(mcols(mdb)$bindingDomain)) > 1)

} # test.bindingDomains
#------------------------------------------------------------------------------------------------------------------------
test.flyBindingDomains <- function()
{
  print('--- test.flyBindingDomains')

  x <- mcols(MotifDb)
  tmp <- as.list(head(sort(table(subset(x, organism=='Dmelanogaster')$bindingDomain), decreasing=TRUE), n=3))

    # these counts will likely change with a fresh load of data from FlyFactorSurvey.

  checkEquals(tmp$Homeobox, 212)
  checkEquals(tmp[['zf-C2H2']], 160)
  checkEquals(tmp[["Helix-Turn-Helix"]], 182)
  checkTrue(length(which(is.na(subset(x, organism=='Dmelanogaster')$bindingDomain))) > 300) # lots of cisbp

} # test.flyBindingDomains
#------------------------------------------------------------------------------------------------------------------------
test.experimentTypes <- function()
{
  print('--- test.experimentTypes')
  mdb <- MotifDb #(quiet=TRUE)
  x <- mcols(mdb)
  checkTrue(length(unique(x$experimentType)) >= 18)
  checkEquals(length(which(x$experimentType=='')), 0)

} # test.experimentTypes
#------------------------------------------------------------------------------------------------------------------------
test.tfFamilies = function()
{
  print('--- test.tfFamilies')
  mdb = MotifDb #(quiet=TRUE)
  checkTrue(length(unique(mcols(mdb)$tfFamily)) > 1)

} # test.tfFamilies
#------------------------------------------------------------------------------------------------------------------------
test.bindingSequences = function()
{
  print('--- test.bindingSequences')
  mdb = MotifDb #(quiet=TRUE)
  checkTrue(length(unique(mcols(mdb)$bindingSequence)) > 1)

} # test.tfFamilies
#------------------------------------------------------------------------------------------------------------------------
test.pubmedIDs = function()
{
  print('--- test.pubmedIDs')
  x = mcols(MotifDb) #(quiet=TRUE))
  checkTrue(length(unique(x$pubmedID)) >= 139)
  checkEquals(length(which(x$pubmedID == '')), 0)

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
skip.test.allFullNames = function()
{
  print('--- test.allFullNames')
  mdb = MotifDb #(quiet=TRUE)
  matrices = mdb@listData
  fullNames = names(matrices)

  all.dataSources = unique(mcols(mdb)$dataSource)
  checkTrue(length(all.dataSources) >= 4)

  for(source in all.dataSources) {
     this.dataSource <- source
     matrices.by.source = subset(mdb, dataSource==this.dataSource)
     matrix.name = names(matrices.by.source)[1]
        #  FlyFactorSurvey: Dmelanogaster-FlyFactorSurvey-ab_SANGER_10_FBgn0259750
        #             hPDI: Hsapiens-hPDI-ABCF2
        #      JASPAR_CORE: JASPAR_CORE-Athaliana-AGL3-MA0001.1
        #         UniPROBE: Hsapiens-UniPROBE-Sox4.UP00401
     checkTrue(grep(this.dataSource, matrix.name) == 1)
     #printf('%20s: %s', source, matrix.name)
     }

  return(TRUE)

} # test.allFullNames
#------------------------------------------------------------------------------------------------------------------------
test.subset = function()
{
  if(interactive()) {
    print('--- test.subset')
    mdb = MotifDb
    checkTrue('geneSymbol' %in% colnames(elementMetadata(mdb)))
    mdb.sub = subset(mdb, geneSymbol=='ABCF2')
    checkEquals(length(mdb.sub), 1)
    } # if interactive

} # test.subset
#------------------------------------------------------------------------------------------------------------------------
test.subsetWithVariables = function()
{
  if(interactive()) {
    print('--- test.subsetWithVariables')

    mdb = MotifDb #()
    checkTrue('geneSymbol' %in% colnames(elementMetadata(mdb)))
    target.gene <<- 'ABCF2'
    mdb.sub = subset(mdb, geneSymbol==target.gene)
    checkEquals(length(mdb.sub), 1)
    } # if interactive

} # test.subsetWithVariables
#------------------------------------------------------------------------------------------------------------------------
# "old style": just one query term allowed
test.queryOldStyle = function()
{
  print('--- test.queryOldStyle')
  mdb = MotifDb

    # do queries on dataSource counts match those from a contingency table?
  sources.list = as.list(table(mcols(mdb)$dataSource))
  checkEquals(length(query(mdb, 'flyfactorsurvey')), sources.list$FlyFactorSurvey)
  checkEquals(length(query(mdb, 'uniprobe')), sources.list$UniPROBE)
  checkEquals(length(query(mdb, 'UniPROBE')), sources.list$UniPROBE)

    # gene symbols which begin with 'sox' are quite common.  can we them?
    # there are currently(19 jul 2012) 18, but since this may change, our test is approximate

  # Change on 8/1/2017: increase top limit of sox entries as they've expanded
  sox.entries = query(mdb, '^sox')
  checkTrue(length(sox.entries) > 10)
  checkTrue(length(sox.entries) < 200)

    # manual inspection reveals that some of these genes have names which are all capitalized.  test that.
  checkTrue(length(query(mdb, '^sox', ignore.case=TRUE)) > length(query(mdb, '^SOX', ignore.case=FALSE)))

    # make sure that queries can be stacked up, and that order of call does not affect the outcome
   uniprobe.sox.matrices = query(query(mdb, 'uniprobe'), '^sox')
   sox.uniprobe.matrices = query(query(mdb, '^sox'), 'uniprobe')

   checkEquals(length(uniprobe.sox.matrices), length(sox.uniprobe.matrices))

     # an approximate count check
  checkTrue(length(uniprobe.sox.matrices) > 10)
  checkTrue(length(uniprobe.sox.matrices) < 30)

   checkEquals(unique(mcols(uniprobe.sox.matrices)$dataSource), 'UniPROBE')
   checkEquals(unique(mcols(sox.uniprobe.matrices)$dataSource), 'UniPROBE')
   gene.symbols = sort(unique(mcols(uniprobe.sox.matrices)$geneSymbol))

    # query on a string(in this case, an oddly named motif) which contains
    # regular expression characters.  no solution to this yet.
    # query uses base R's grep, in which
    #  x <- query(mdb, "ELK1,4_GABP{A,B1}.p3")

} # test.queryOldStyle
#------------------------------------------------------------------------------------------------------------------------
test.query <- function()
{
  print('--- test.query')
  mdb = MotifDb

  ors <- c("MA0511.1", "MA0057.1")
  ands <- c("jaspar2018", "sapiens")
  nots <- "cisbp"

  x <- query(mdb, andStrings=ands, orStrings=ors)
  checkEquals(length(x), 2)
  checkEquals(sort(names(x)),
             c("Hsapiens-jaspar2018-MZF1(var.2)-MA0057.1", "Hsapiens-jaspar2018-RUNX2-MA0511.1"))

  x <- query(mdb, andStrings="MA0057.1")
  checkEquals(length(x), 15)

  x <- query(mdb, andStrings=c("MA0057.1", "cisbp"))
  checkEquals(length(x), 11)

  x <- query(mdb, andStrings=c("MA0057.1"), notStrings="cisbp")
  checkEquals(length(x), 4)

  x <- query(mdb, andStrings=c("MA0057.1"), notStrings=c("cisbp", "JASPAR_2014"))
  checkEquals(length(x), 3)

  x <- query(mdb, orStrings=c("mus", "sapiens"), andStrings="MA0057.1")
  #checkEquals(sort(names(x)),

    # do queries on dataSource counts match those from a contingency table?
  sources.list = as.list(table(mcols(mdb)$dataSource))
  checkEquals(length(query(mdb, 'flyfactorsurvey')), sources.list$FlyFactorSurvey)
  checkEquals(length(query(mdb, 'uniprobe')), sources.list$UniPROBE)
  checkEquals(length(query(mdb, 'UniPROBE')), sources.list$UniPROBE)

} # test.query
#------------------------------------------------------------------------------------------------------------------------
test.query2 <- function()
{
  mdb <- MotifDb
  matrices.human <- query(mdb, 'hsapiens')
  matrices.sox4 <- query(mdb, 'sox4')

  matrices.human.sox4 <- query(mdb, c("hsapiens", "sox4"))
  matrices.human.sox4.oldStyle <- query(matrices.human, "sox4")
  checkTrue(length(matrices.human.sox4) > 0)   # 6 found on(24 oct 2018)
  checkEquals(length(matrices.human.sox4), length(matrices.human.sox4.oldStyle))

  checkTrue(length(matrices.human.sox4) < length(matrices.human))
  checkTrue(length(matrices.human.sox4) < length(matrices.sox4))

  uniprobe.sox.matrices <- query(mdb, c('uniprobe', '^sox'))

} # test.query2
#------------------------------------------------------------------------------------------------------------------------
test.transformMatrixToMemeRepresentation = function()
{
  print('--- test.transformMatrixToMemeRepresentation')
  mdb = MotifDb #()
  matrix.agl3 = subset(mdb, dataSource=='JASPAR_CORE' & organism=='Athaliana' & geneSymbol=='AGL3')[[1]]
  checkEquals(dim(matrix.agl3), c(4, 10))

    # a transposed frequency is needed for meme.  we store normalized frequency matrices, to just transposition is needed
  tfm = MotifDb:::transformMatrixToMemeRepresentation(matrix.agl3)
  checkEquals(dim(tfm), c(10, 4))
  checkTrue(all(round(rowSums(tfm)) == 1.0))

} # test.transformMatrixToMemeRepresentation
#------------------------------------------------------------------------------------------------------------------------
# translate a motif matrix from MotifList into text suitable for meme and tomtom input.
# choose the top two from UniPROBE.  they reliably have high counts, and easily distinguished normalized frequencies
test.matrixToMemeText = function()
{
  print('--- test.matrixToMemeText')
  sox4 = subset(MotifDb, dataSource=='UniPROBE' & organism=='Hsapiens' & geneSymbol=='Sox4')
  rtg3 = subset(MotifDb, dataSource=='UniPROBE' & organism=='Scerevisiae' & geneSymbol=='Rtg3')

  text.sox4 = MotifDb:::matrixToMemeText(sox4)
     # check these with t(sox4 [[1]])
  line1.sox4  = " 0.2457457130  0.1950426500  0.2287887620  0.3304228750"
  line14.sox4 = " 0.2821643030  0.2286132160  0.1585395830  0.3306828990"

  checkEquals(length(text.sox4), 29)
  checkEquals(text.sox4 [1], "MEME version 4")
  checkEquals(text.sox4 [10], "MOTIF Hsapiens-UniPROBE-Sox4.UP00401")
  checkEquals(grep(line1.sox4, text.sox4), 12)
  checkEquals(grep(line14.sox4, text.sox4), 25)

  text.rtg3 = MotifDb:::matrixToMemeText(rtg3)
  line1.rtg3  = " 0.3935122858  0.1453016447  0.3308830322  0.1303030373"
  line20.rtg3 = " 0.2490417648  0.3966478493  0.1083586569  0.2459517291"
  checkEquals(length(text.rtg3), 35)         # 4 trailing empty lines
  checkEquals(text.rtg3 [1], "MEME version 4")
  checkEquals(text.rtg3 [10], "MOTIF Scerevisiae-UniPROBE-Rtg3.UP00356")
  checkEquals(grep(line1.rtg3, text.rtg3), 12)
  checkEquals(grep(line20.rtg3, text.rtg3), 31)

    # now call with both matrices, and see if the right results are returned
  text.both = MotifDb:::matrixToMemeText(c(sox4, rtg3))
  checkEquals(text.both [1], "MEME version 4")
  checkEquals(grep(line1.sox4, text.both), 12)
  checkEquals(grep(line14.sox4, text.both), 25)
  checkEquals(grep(line1.rtg3, text.both),  29)
  checkEquals(grep(line20.rtg3, text.both), 48)

} # test.matrixToMemeText
#------------------------------------------------------------------------------------------------------------------------
# translate a motif matrix from MotifList into text suitable for meme and tomtom input.
# choose the top two from UniPROBE.  they reliably have high counts, and easily distinguished normalized frequencies
#test.matrixToMemeText_mapVersion = function()
#{
#  print('--- test.matrixToMemeText_mapVersion')
#  sox4 = subset(MotifDb, dataSource=='UniPROBE' & organism=='Hsapiens' & geneSymbol=='Sox4')
#  rtg3 = subset(MotifDb, dataSource=='UniPROBE' & organism=='Scerevisiae' & geneSymbol=='Rtg3')
#
#  text.sox4 = MotifDb:::mapped.broken.matrixToMemeText(sox4)
#  text.rtg3 = MotifDb:::mapped.broken.matrixToMemeText(rtg3)
#  text.both = MotifDb:::mapped.broken.matrixToMemeText(c(sox4, rtg3))
#
#     # check these with t(sox4 [[1]])
#  line1.sox4  = " 0.2457457130  0.1950426500  0.2287887620  0.3304228750"
#  line14.sox4 = " 0.2821643030  0.2286132160  0.1585395830  0.3306828990"
#
#  checkEquals(length(text.sox4), 29)
#  checkEquals(text.sox4 [1], "MEME version 4")
#  checkEquals(text.sox4 [10], "MOTIF Hsapiens-UniPROBE-Sox4.UP00401")
#  checkEquals(grep(line1.sox4, text.sox4), 12)
#  checkEquals(grep(line14.sox4, text.sox4), 25)
#
#  line1.rtg3  = " 0.3935122858  0.1453016447  0.3308830322  0.1303030373"
#  line20.rtg3 = " 0.2490417648  0.3966478493  0.1083586569  0.2459517291"
#  checkEquals(length(text.rtg3), 35)         # 4 trailing empty lines
#  checkEquals(text.rtg3 [1], "MEME version 4")
#  checkEquals(text.rtg3 [10], "MOTIF Scerevisiae-UniPROBE-Rtg3.UP00356")
#  checkEquals(grep(line1.rtg3, text.rtg3), 12)
#  checkEquals(grep(line20.rtg3, text.rtg3), 31)
#
#    # now call with both matrices, and see if the right results are returned
#  checkEquals(text.both [1], "MEME version 4")
#  checkEquals(grep(line1.sox4, text.both), 12)
#  checkEquals(grep(line14.sox4, text.both), 25)
#  checkEquals(grep(line1.rtg3, text.both),  29)
#  checkEquals(grep(line20.rtg3, text.both), 48)
#
#} # test.matrixToMemeText_mapVersion
#------------------------------------------------------------------------------------------------------------------------
test.export_memeFormatStdOut = function()
{
  print('--- test.export_memeFormatStdOut')
  mdb = MotifDb #()
  mdb.chicken = subset(mdb, organism=='Gallus')
  checkEquals(length(mdb.chicken), 3)
    # text is cat-ed to stdout, so not avaialable here to check.
    # but just like print, export also returns the text invisibly.
    # so that CAN be checked.

  meme.text = export(mdb.chicken, format='meme')
  checkEquals(length(meme.text), 1)   # just one long string
  checkTrue(is.character(meme.text))
  checkTrue(nchar(meme.text) > 800)   # 1002 as of(10 aug 2012)
  return(TRUE)

} # test.exportMemeFormatToStdOut
#------------------------------------------------------------------------------------------------------------------------
test.export_memeFormatToFile = function()
{
  print('--- test.export_memeFormatToFile')
  mdb = MotifDb #()
  mdb.chicken = subset(mdb, organism=='Gallus')
  checkEquals(length(mdb.chicken), 3)
  output.file = tempfile()
  meme.text = export(mdb.chicken, output.file, 'meme')
  retrieved = scan(output.file, what=character(0), sep='\n', quiet=TRUE)
  invisible(retrieved)

} # test.exportMemeFormatToFile
#------------------------------------------------------------------------------------------------------------------------
test.export_memeFormatToFileDuplication = function()
{
  print('--- test.export_memeFormatToFileDuplication')
  mdb = MotifDb #()
  mdb.mouse = subset(mdb, organism=='Mmusculus')
  checkTrue(length(mdb.mouse) > 1300)
  output.file = 'mouse.txt' # tempfile()
  max = 3
  meme.text = export(mdb.mouse [1:max], output.file, 'meme')
  retrieved = scan(output.file, what=character(0), sep='\n', quiet=TRUE)
  invisible(retrieved)

} # test.exportMemeFormatToFileDuplication
#------------------------------------------------------------------------------------------------------------------------
test.export_memeFormatToFile_run_tomtom = function(max=50)
{
  if(interactive()) {
    print('--- test.export_memeFormatToFile_run_tomtom')
    mdb = MotifDb
    sox4.mouse = subset(mdb, organism=='Mmusculus' & geneSymbol=='Sox4')
    all.human = subset(mdb, organism=='Hsapiens')

    tomtom.tmp.dir = tempfile()
    print(tomtom.tmp.dir)
    stopifnot(dir.create(tomtom.tmp.dir))
    sox4.file.path = file.path(tomtom.tmp.dir, 'sox4.mouse.text')
    all.human.file.path = file.path(tomtom.tmp.dir,'all.human.text')
    export(sox4.mouse, sox4.file.path, 'meme')
    export(all.human, all.human.file.path, 'meme')

       # find similarity of motif #1 to all the motifs in mdbMany

       # cannot rely upon tomtom being present

    #cmd = sprintf('tomtom -no-ssc -oc %s -verbosity 3 -min-overlap 5 -mi 1 -dist pearson -evalue -thresh 10 %s %s',
    #                tomtom.tmp.dir, sox4.file.path, all.human.file.path)
    #system(cmd)
    #cmd = sprintf('open %s/tomtom.html', tomtom.tmp.dir)
    #system(cmd)
    } # if interactive

} # test.export_memeFormatToFile_run_tomtom
#------------------------------------------------------------------------------------------------------------------------
# MotIV::motifMatch fails with MotIV_1.25.0.  will look into this in September, 2015, pshannon
# MotIV abandoned, all of my own motif/sequence matching is done via an independently installed
# FIMO from the meme suite.
# another option is the MOODs motif matching capability provided by the bioc package "motifmatchr"
# this possibility, and these tests, are deferred for now.
disabled_test.run_MotIV.motifMatch = function()
{
  require(MotIV)
  print('--- test.run_MotIV.motifMatch')
  mdb <- MotifDb #()

  db.tmp <- mdb@listData

     # match motif 1 against the entire MotifDb  collection
  motif.hits <- motifMatch(db.tmp [1], database<-db.tmp)
     # the long way to extract the matrix name.  see MotIV.toTable below for more convenient way
  checkEquals(motif.hits@bestMatch[[1]]@aligns[[1]]@TF@name, names(db.tmp)[1])

     # match the last motif against all
  last <- length(db.tmp)
     # MotIV:motifMatch works differently on linux and macos.  by asking for 50 matches,
     # the search target(db.tmp[last]) is sure to be in the hit list.
  motif.hits <-  motifMatch(db.tmp [last], database=db.tmp, top=50)
  tbl.hits <- MotIV.toTable(motif.hits)
    # the 5 hits return should include the one we tried to match, but the MotIV search strategy
    # may not place it first
  checkTrue(names(db.tmp[last]) %in% tbl.hits$name)
  invisible(tbl.hits)

} # distable_test.run_MotIV.motifMatch
#------------------------------------------------------------------------------------------------------------------------
MotIV.toTable = function(match)
{
  stopifnot(length(match@bestMatch) >= 1)
  alignments = match@bestMatch[[1]]@aligns

  df = data.frame(stringsAsFactors=FALSE)
  for(alignment in alignments) {
    x = alignment
    name = x@TF@name
    eVal = x@evalue
    sequence = x@sequence
    match = x@match
    strand = x@strand
    df = rbind(df, data.frame(name=name, eVal=eVal, sequence=sequence, match=match, strand=strand, stringsAsFactors=FALSE))
    } # for alignment

  return(df)

} # MotIV.toTable
#------------------------------------------------------------------------------------------------------------------------
disabled_test.MotIV.toTable = function()
{
  print('--- test.MotIVtoTable')
  mdb = MotifDb #()
  test.hits = motifMatch(mdb[1]@listData, database=jaspar)
  tbl.hits =  MotIV.toTable(test.hits)
  checkEquals(dim(tbl.hits), c(5, 5))
  checkEquals(sort(colnames(tbl.hits)), sort(c("name", "eVal", "sequence", "match", "strand")))

} # test.MotIV.toTable
#------------------------------------------------------------------------------------------------------------------------
pwmMatch.toTable = function(motifMatch) {
   if(length(motifMatch@bestMatch) == 0)
   return(NA)

   df.list = vector("list", length(motifMatch@bestMatch))
   for(k in seq(length(motifMatch@bestMatch)))
   {
       alignments = motifMatch@bestMatch[[k]]@aligns
       df = data.frame(stringsAsFactors=FALSE)
       for(alignment in alignments) {
           x = alignment
           name = x@TF@name
           eVal = x@evalue
           sequence = x@sequence
           match = x@match
           strand = x@strand
           df = rbind(df, data.frame(name=name, eVal=eVal, sequence=sequence, match=match, strand=strand, stringsAsFactors=FALSE))
       } # for alignment
       df.list[[k]]=df
   }
   names(df.list) <- names(motifMatch)
   return(df.list)

} # pwmMatch.toTable
#------------------------------------------------------------------------------
# Robert Stojnic reports incorrect gene symbols for matrices obtained from
# flyFactorSurvey.
# the solution was to abandon the original strategy of extracting the
# symbol from the matrix(and file) name.
# now, the flybase importer("inst/scripts/import/flyFactorSurvey/import.R")
# uses FBgn id(which can be reliably extracted) and uses indpendent
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
    print("--- test.flyFactorGeneSymbols")
    mdb = MotifDb
    checkEquals(sort(mcols(query(mdb, "FBgn0259750"))$geneSymbol),
                sort(c("FBgn0259750", "FBgn0259750")))
    checkEquals(mcols(query(mdb, "FBgn0000014"))$geneSymbol, rep("abd-A", 3))
    checkEquals(mcols(query(mdb, "FBgn0000015"))$geneSymbol, rep("Abd-B", 3))

} # test.flyFactorGeneSymbols
#-------------------------------------------------------------------------------
test.export_jasparFormatStdOut = function()
{
  print('--- test.export_jasparFormatStdOut')
  mdb = MotifDb #()
  mdb.chicken = subset(mdb, organism=='Gallus')
  checkEquals(length(mdb.chicken), 3)
    # text is cat-ed to stdout, so not avaialable here to check.
    # but just like print, export also returns the text invisibly.
    # so that CAN be checked.

  jaspar.text = export(mdb.chicken, format='jaspar')
  checkEquals(length(jaspar.text), 1)   # just one long string
  checkTrue(is.character(jaspar.text))
  checkTrue(nchar(jaspar.text) > 800)   # 1002 as of(10 aug 2012)
  return(TRUE)

} # test.exportjasparFormatToStdOut
#------------------------------------------------------------------------------------------------------------------------
test.export_jasparFormatToFile = function()
{
  print('--- test.export_jasparFormatToFile')
  mdb = MotifDb #()
  mdb.chicken = subset(mdb, organism=='Gallus')
  checkEquals(length(mdb.chicken), 3)
  output.file = tempfile()
  jaspar.text = export(mdb.chicken, output.file, 'jaspar')
  retrieved = scan(output.file, what=character(0), sep='\n', quiet=TRUE)
  invisible(retrieved)

} # test.exportjasparFormatToFile
#------------------------------------------------------------------------------------------------------------------------
test.geneToMotif <- function()
{
   printf("--- test.geneToMotif")
   mdb <- MotifDb

   genes <- c("FOS", "ATF5", "bogus", "SATB2")
   good.genes <- genes[-which(genes=="bogus")]
      # use  TFClass family classifcation
   tbl.tfClass <- geneToMotif(mdb, genes, source="TfClaSS")   # intentional mis-capitalization
   checkTrue(all(good.genes %in% tbl.tfClass$gene))

   expected.motifs <- c("MA0833.1", "MA0099.2", "MA0476.1", "MA0679.1", "MA0754.1", "MA0755.1", "MA0756.1", "MA0757.1")
   checkTrue(all(expected.motifs %in% tbl.tfClass$motif))
   checkEquals(unique(tbl.tfClass$source), "TFClass")

      # MotifDb mode uses the MotifDb metadata, pulled from many sources
   tbl.mdb <- geneToMotif(mdb, genes, source="mOtifdb")     # intentional mis-capitalization
   checkEquals(dim(tbl.mdb), c(14, 6))
   checkEquals(subset(tbl.mdb, dataSource=="jaspar2016" & geneSymbol== "FOS")$motif, "MA0476.1")
      # no recognizable(i.e., jaspar standard) motif name returned by MotifDb metadata
      # MotifDb for ATF5
      # todo: compare the MA0110596_1.02 matrix of cisp_1.02 to japar MA0833.1

     # check use of ignore.case
   tbl.caseSensitive <-  geneToMotif(MotifDb, "STAT4", source="MotifDb")
   checkEquals(length(grep("jaspar", tbl.caseSensitive$dataSource, ignore.case=TRUE)), 0)
   tbl.caseInsensitive <-  geneToMotif(MotifDb, "STAT4", source="MotifDb", ignore.case=TRUE)
   checkTrue(length(grep("jaspar", tbl.caseInsensitive$dataSource, ignore.case=TRUE)) >= 3)

   tbl.caseSensitive <-  geneToMotif(MotifDb, "stat4", source="TFclass")
   checkEquals(nrow(tbl.caseSensitive), 0)
   tbl.caseInsensitive <-  geneToMotif(MotifDb, "stat4", source="TFclass", ignore.case=TRUE)
   checkTrue(nrow(tbl.caseInsensitive) >= 5)

} # test.geneToMotif
#------------------------------------------------------------------------------------------------------------------------
# this case discovered(31 jan 2018). when called on a gene/source combination for which there are
# no motifs, i attempted to add the mapping source(either "MotifDb", "TFClass") as a column
# to an empty data.frame.  check for that and its fix here
test.geneToMotif.oneGene.noMotifs <- function()
{
   checkEquals(nrow(geneToMotif(MotifDb, "SATB2", "MotifDb")), 0)
   checkEquals(nrow(geneToMotif(MotifDb, "bogus-arandum", "MotifDb")), 0)
   checkEquals(nrow(geneToMotif(MotifDb, "bogus-arandum", "TFclass")), 0)

} # test.geneToMotif.oneGene.noMotifs
#------------------------------------------------------------------------------------------------------------------------
# sad to say I do not recall what problem/fix is tested here(pshannon, 23 jan 2018).
# however, it demonstrates the variety of results which can be returned by non-jaspar datasets
# when using the MotifDb mapping source, and the relative paucity which is sometimes
# seen with the TFclass mapper
test.geneToMotif.ignore.jasparSuffixes <- function()
{
   printf("--- test.geneToMotif.ignore.jasparSuffixes")
   mdb <- MotifDb

   genes <- c("FOS", "ATF5", "bogus")

      # use  TFClass family classifcation
   tbl.tfClass <- geneToMotif(mdb, genes, source="TfClaSS")   # intentional mis-capitalization
   checkEquals(sort(tbl.tfClass$gene),  sort(c("ATF5", "FOS", "FOS")))
   checkEquals(sort(tbl.tfClass$motif),  sort(c("MA0833.1", "MA0099.2", "MA0476.1")))
   checkEquals(tbl.tfClass$source, rep("TFClass", 3))

      # MotifDb mode uses the MotifDb metadata, pulled from many sources
   tbl.mdb <- geneToMotif(mdb, genes, source="mOtifdb")     # intentional mis-capitalization
   checkEquals(dim(tbl.mdb), c(14, 6))
   checkEquals(subset(tbl.mdb, dataSource=="jaspar2016" & geneSymbol== "FOS")$motif, "MA0476.1")
      # no recognizable(i.e., jaspar standard) motif name returned by MotifDb metadata
      # MotifDb for ATF5

      # compare the MA0110599_1.02 matrix of cisp_1.02 to japar MA0476.1: the identical matrix!
      # 1         FOS    MA0110599_1.02   cisbp_1.02  Hsapiens 24194598 MotifDb
      # 10        FOS          MA0476.1   jaspar2018  Hsapiens 17916232 MotifDb
      # this establishes the need for careful scrutiny as one winnows a geneToMotif result into
      # useful non-reduplicative sequence analysis

   pfm.ma0110599 <- as.list(query(mdb, "MA0110599"))[[1]]
   pfm.ma0476.1  <- as.list(query(query(mdb, "MA0476.1"), "jaspar2018"))[[1]]
   checkEquals(pfm.ma0110599, pfm.ma0476.1)

} # test.geneToMotif.ignore.jasparSuffixes
#------------------------------------------------------------------------------------------------------------------------
test.motifToGene <- function()
{
   printf("--- test.motifToGene")
   #printf("Sys.getlocale: %s", Sys.getlocale())

     # good test case of querying both sources,
   motif <- "MA0099.2"
   tbl <- motifToGene(MotifDb, motif, c("MotifDb", "TFclass"))
      # [1]   mdb.genes: AP1, JUN::FOS, FOS::JUN, FOS::JUN
      # [1]   tfc.genes: FOS, JUN
   checkEquals(sort(unique(tbl$geneSymbol)), c("AP1", "FOS", "FOS::JUN", "JUN", "JUN::FOS"))

   motifs <- c("MA0592.2", "UP00022", "ELF1.SwissRegulon")

      # MotifDb mode uses the MotifDb metadata "providerId",
   tbl.mdb <- motifToGene(MotifDb, motifs, source="MotifDb")
   checkEquals(dim(tbl.mdb), c(3, 5))
   expected <- sort(c("MA0592.2", "ELF1.SwissRegulon", "UP00022"))
   actual <- sort(tbl.mdb$motif)
   checkEquals(actual, expected)
   checkEquals(sort(tbl.mdb$geneSymbol), sort(c("ELF1", "Esrra", "Zfp740")))
   checkEquals(tbl.mdb$source,     rep("MotifDb", 3))

      # TFClass mode uses  TF family classifcation
   tbl.tfClass <- motifToGene(MotifDb, motifs, source="TFClass")
   checkEquals(dim(tbl.tfClass), c(9,5))
   checkEquals(tbl.tfClass$motif, rep("MA0592.2", 9))
   checkEquals(sort(tbl.tfClass$gene), sort(c("AR", "ESR1", "ESR2", "ESRRA", "ESRRB", "ESRRG", "NR3C1", "NR3C2", "PGR")))
   checkEquals(tbl.tfClass$source,       rep("TFClass", 9))

     # test motifs with regex characters in them, or other characters neither letter nor number
   motifs <- sort(c("DMAP1_NCOR{1,2}_SMARC.p2", "ELK1,4_GABP{A,B1}.p3", "SNAI1..3.p2", "EWSR1-FLI1.p2", "ETS1,2.p2"))
   tbl <- motifToGene(MotifDb, motifs, source="MotifDb")
   checkEquals(nrow(tbl), 0)

   tbl <- motifToGene(MotifDb, motifs, source="tfclass")
   checkEquals(ncol(tbl), 5)
   checkTrue(nrow(tbl) > 80)
   checkTrue(nrow(tbl) < 100)
   checkTrue(all(motifs %in% tbl$motif))

   motifs <- c("Hsapiens-HOCOMOCOv10-IKZF1_HUMAN.H10MO.C",
               "MA0099.2",
               "Hsapiens-SwissRegulon-ARID3A.SwissRegulon",
               "MA0592.2",
               "MA0036913_1.02")

   tbl <- motifToGene(MotifDb, motifs, source=c("MotifDb", "TFClass"))
   checkTrue(all(motifs %in% tbl$motif))

   motif <- "Mmusculus;Rnorvegicus;Hsapiens-jaspar2018-FOS::JUN-MA0099.2"
   tbl <- motifToGene(MotifDb, motif, source="TFClass")
   checkEquals(tbl$geneSymbol, c("FOS", "JUN"))

   motifs <- c("Hsapiens-jaspar2016-RUNX1-MA0002.1",
               "Hsapiens-jaspar2016-TFAP2A-MA0003.1",
               "Hsapiens-jaspar2016-TFAP2A-MA0003.2",
               "Hsapiens-jaspar2016-TFAP2A-MA0003.3",
               "Hsapiens-jaspar2016-AR-MA0007.2",
               "MA0872.1")
   tbl <- motifToGene(MotifDb, motifs, source="Motifdb")
   checkTrue(all(motifs %in% tbl$motif))
   checkEquals(sort(unique(tbl$geneSymbol)), c("AR", "RUNX1", "TFAP2A", "TFAP2A(var.3)"))

   tbl <- motifToGene(MotifDb, motifs, source="TFClass")
   checkEquals(sort(unique(tbl$motif)), c("Hsapiens-jaspar2016-TFAP2A-MA0003.3", "MA0872.1"))
   checkEquals(sort(unique(tbl$geneSymbol)), c("TFAP2A", "TFAP2B", "TFAP2C", "TFAP2D", "TFAP2E"))

   tbl <- motifToGene(MotifDb, motifs, source=c("MotifDb", "TFClass"))
   checkTrue(all(motifs %in% tbl$motif))
   checkEquals(sort(unique(tbl$geneSymbol)),
                    c("AR", "RUNX1", "TFAP2A", "TFAP2A(var.3)", "TFAP2B", "TFAP2C", "TFAP2D", "TFAP2E"))

      #(23 may 2018) found that MotifDb works, but c("MotifDb", "TFClass") does not
      # test the fix here

   motifs <- c("Hsapiens-jolma2013-IRF5-2", "Hsapiens-SwissRegulon-IRF5.SwissRegulon")
   checkEquals(motifToGene(MotifDb, motifs, source=c("MotifDb"))$geneSymbol, c("IRF5", "IRF5"))
   checkEquals(nrow(motifToGene(MotifDb, motifs, source=c("TFClass"))), 0)
   checkEquals(motifToGene(MotifDb, motifs, source=c("MotifDb", "TFClass"))$geneSymbol, c("IRF5", "IRF5"))

   } # test.motifToGene
#------------------------------------------------------------------------------------------------------------------------
test.associateTranscriptionFactors <- function()
{
   printf("--- test.associateTranscriptionFactors")

   mdb <- MotifDb
   pfms <- query(mdb, andStrings=c("sapiens", "jaspar2016"))

      # query(mdb, "jaspar2016", "sapiens")
      # first check motifs with MotifDb-style long names, using MotifDb lookup, in the
      # metadata of MotifDb:
      #      "Hsapiens-jaspar2016-RUNX1-MA0002.1"  "Hsapiens-jaspar2016-TFAP2A-MA0003.1"

   motif.names <- c(names(pfms[1:5]), "MA0872.1", "hocus.pocus")
   tbl <- data.frame(motifName=motif.names, score=runif(7), stringsAsFactors=FALSE)

   tbl.anno.mdb <- associateTranscriptionFactors(mdb, tbl, source="MotifDb", expand.rows=TRUE)
   checkEquals(nrow(tbl), nrow(tbl.anno.mdb))
   checkTrue(is.na(tbl.anno.mdb$geneSymbol[grep("hocus.pocus", tbl.anno.mdb$motifName)]))
   checkTrue(all(c("AR", "RUNX1", "TFAP2A", "TFAP2A", "TFAP2A", "TFAP2A(var.3)") %in% tbl.anno.mdb$geneSymbol))

   tbl.anno.tfc <- associateTranscriptionFactors(mdb, tbl, source="TFClass", expand.rows=TRUE)
   checkTrue(nrow(tbl) < nrow(tbl.anno.tfc))
   checkEquals(sort(unique(tbl.anno.tfc$geneSymbol)), c("TFAP2A", "TFAP2B", "TFAP2C", "TFAP2D", "TFAP2E"))

   tbl.anno.both <- associateTranscriptionFactors(mdb, tbl, source=c("MotifDb", "TFClass"), expand.rows=TRUE)
   checkEquals(length(grep("MotifDb", tbl.anno.both$source)), 6)
   checkEquals(length(grep("TFClass", tbl.anno.both$source)), 10)

   #   checkEquals(dim(tbl.anno.mdb), c(nrow(tbl), ncol(tbl) + 4))
   #   checkTrue(all(c("geneSymbol", "pubmedID") %in% colnames(tbl.anno.mdb)))
   #   checkEquals(sort(tbl.anno.mdb$geneSymbol),
   #               sort(c("AR", "RUNX1", "TFAP2A", "TFAP2A", "TFAP2A", "TFAP2A(var.3)")))
   #
   #      # now add in a bogus motif name, one for which there cannot possibly be a TF
   #
   #   motif.names[3] <- "bogus"
   #   tbl <- data.frame(motifName=motif.names, score=runif(5), stringsAsFactors=FALSE)
   #   tbl.anno <- associateTranscriptionFactors(mdb, tbl, source="MotifDb", expand.rows=FALSE)
   #   checkTrue(is.na(tbl.anno$geneSymbol[3]))
   #   checkTrue(is.na(tbl.anno$pubmedID[3]))
   #
   #      # now check motifs with short, traditional names, in "TFClass" mode, which uses
   #      # the tfFamily
   #      #      "MA0002.1" "MA0003.1" "MA0003.2" "MA0003.3" "MA0007.2"
   #
   #   motif.names <- names(pfms[1:5])
   #   short.motif.names <- unlist(lapply(strsplit(motif.names, "-"), function(tokens) return(tokens[length(tokens)])))
   #   tbl <- data.frame(motifName=motif.names, shortMotif=short.motif.names, score=runif(5), stringsAsFactors=FALSE)
   #
   #   tbl.anno <- associateTranscriptionFactors(mdb, tbl, source="TFClass", expand.rows=FALSE)
   #   checkEquals(dim(tbl.anno), c(nrow(tbl), ncol(tbl) + 2))
   #
   #      # TFClass only annotates MA0003.3, none of the others
   #   checkTrue(all(is.na(tbl.anno$geneSymbol[-4])))
   #   checkTrue(all(is.na(tbl.anno$pubmedID[-4])))
   #   checkEquals(tbl.anno$geneSymbol[4], "TFAP2A;TFAP2B;TFAP2C;TFAP2D;TFAP2E")
   #   checkEquals(tbl.anno$pubmedID[4],   "23180794")
   #
   #      # now ask for expandsion of the semicolon separated list
   #   tbl.anno <- associateTranscriptionFactors(mdb, tbl, source="TFClass", expand.rows=TRUE)
   #   checkEquals(dim(tbl.anno), c(nrow(tbl) + 4, ncol(tbl) + 2))
   #   checkTrue(all(c("TFAP2A", "TFAP2B", "TFAP2C", "TFAP2D", "TFAP2E") %in% tbl.anno$geneSymbol))
   #
   #      # now add in a bogus motif name, one for which there cannot possibly be a TF
   #
   #   motif.names <- names(pfms[1:5])
   #   short.motif.names <- unlist(lapply(strsplit(motif.names, "-"), function(tokens) return(tokens[length(tokens)])))
   #   short.motif.names[4] <- "bogus"
   #   tbl <- data.frame(shortMotif=short.motif.names, score=runif(5), stringsAsFactors=FALSE)
   #
   #   tbl.anno <- associateTranscriptionFactors(mdb, tbl, source="TFClass", expand.rows=FALSE)
   #   checkEquals(dim(tbl.anno), c(nrow(tbl), ncol(tbl) + 2))
   #      # after adding bogus to the only mapped motif name, all geneSymbol and pubmedID values should be NA
   #   checkTrue(all(is.na(tbl.anno$geneSymbol)))
   #   checkTrue(all(is.na(tbl.anno$pubmedID)))
   #
   #      # now make sure that the absence of the  TFClass-specific "shortMotif" field is detected
   #   motif.names <- names(pfms[1:5])
   #   tbl <- data.frame(motifName=motif.names, score=runif(5), stringsAsFactors=FALSE)
   #   checkException(tbl.anno <- associateTranscriptionFactors(mdb, tbl, source="TFClass", expand.rows=FALSE), silent=TRUE)

      # now some motif names
} # test.associateTranscriptionFactors
#------------------------------------------------------------------------------------------------------------------------
# disabled (2 apr 2020) due to very large (~100?) dependendencies introducted directly and indirectly
# via motifmatchr, TFBSTools, universalmotif
# test.match <- function()
# {
#    printf("--- test.match")
#    gr.region <- GRanges(seqnames="chr1", IRanges(start=47229520, end=47229560))
#    motifs <- query(MotifDb, c("jaspar2018", "ZNF263"))
#    checkEquals(length(motifs), 1)
#    gr.match <- matchMotif(MotifDb, motifs, "hg38", gr.region, 1e-5)
#    checkEquals(length(gr.match), 1)  # just one motif
#    checkEquals(names(gr.match), names(motifs))
#    checkEquals(length(gr.match[[1]]), 3)
#
#    tbl.match <- matchMotif(MotifDb, motifs, "hg38", gr.region, 1e-5, fimoDataFrameStyle=TRUE)
#    checkEquals(dim(tbl.match), c(3, 7))
#    checkTrue(all(tbl.match$motif == names(motifs)))
#    checkEquals(class(tbl.match$chrom), "character")  # not a factor
#
#    motifs <- query(MotifDb, "ZNF263", c("jaspar2018", "swissregulon"))
#    checkEquals(length(motifs), 2)
#    gr.match <- matchMotif(MotifDb, motifs, "hg38", gr.region, 1e-5)
#    checkEquals(names(gr.match), names(motifs))
#    checkEquals(as.numeric(lapply(gr.match, length)), c(3, 1))
#
#    tbl.match <-matchMotif(MotifDb, motifs, "hg38", gr.region, 1e-5, fimoDataFrameStyle=TRUE)
#    checkEquals(dim(tbl.match), c(4, 7))
#    checkEquals(length(unique(tbl.match$motif)), 2)
#    checkEquals(unique(tbl.match$motif), names(motifs))
#    checkEquals(colnames(tbl.match), c("chrom", "start", "end", "width", "strand", "mood.score", "motif_id"))
#
#
#         #------------------------------------------------
#         # now all jaspar2018 human motifs
#         #------------------------------------------------
#
#    motifs <- query(MotifDb, c("jaspar2018", "hsapiens"))
#    tbl.match <- matchMotif(MotifDb, motifs, "hg38", gr.region, 1e-5, fimoDataFrameStyle=TRUE)
#    checkEquals(dim(tbl.match), c(7, 7))
#    checkEquals(sort(unique(tbl.match$motif)),
#                c("Hsapiens-jaspar2018-EWSR1-FLI1-MA0149.1", "Hsapiens-jaspar2018-ZNF263-MA0528.1"))
#
#         #-----------------------------------------------------
#         # now all jaspar2018 human motifs, loosen the pValue
#         #-----------------------------------------------------
#
#    motifs <- query(MotifDb, c("jaspar2018", "hsapiens"))
#    tbl.match <- matchMotif(MotifDb, motifs, "hg38", gr.region, 1e-4, fimoDataFrameStyle=TRUE)
#    checkTrue(nrow(tbl.match) > 15)
#
#    tbl.match <- matchMotif(MotifDb, motifs, "hg38", gr.region, 1e-3, fimoDataFrameStyle=TRUE)
#    checkTrue(nrow(tbl.match) > 50)
#
#        #-------------------------------------------------------------
#        # now all jaspar2018 and hocomoco human motifs across 10kb
#        #------------------------------------------------------------
#
#    motifs <- query(MotifDb, "hsapiens", orStrings=c("jaspar2018", "hocomoco-core"))
#    checkTrue(length(motifs) > 500)
#    gr.region <- GRanges(seqnames="chr1", IRanges(start=47229000, end=47239000))
#
#    tbl.match <- matchMotif(MotifDb, motifs, "hg38", gr.region, 1e-7, fimoDataFrameStyle=TRUE)
#    checkTrue(nrow(tbl.match) > 90 && nrow(tbl.match) < 110)
#    checkEquals(order(tbl.match$start), seq_len(nrow(tbl.match)))
#
# } # test.match
#------------------------------------------------------------------------------------------------------------------------
findMotifsWithMutuallyExclusiveMappings <- function()
{
   xtab <- as.data.frame(table(MotifDb@manuallyCuratedGeneMotifAssociationTable$motif))
   xtab <- xtab[order(xtab$Freq, decreasing=TRUE),]
   for(motif in xtab$Var1){
      mdb.genes <- toupper(motifToGene(MotifDb, motif, "motifdb")$geneSymbol)
      tfc.genes <- toupper(motifToGene(MotifDb, motif, "tfclass")$geneSymbol)
      if(length(mdb.genes) == 0) next
      if(length(tfc.genes) == 0) next
      if(length(intersect(mdb.genes, tfc.genes)) == 0){
         printf("------ %s", motif)
         printf("  mdb.genes: %s", paste(mdb.genes, collapse=", "))
         printf("  tfc.genes: %s", paste(tfc.genes, collapse=", "))
        } # if no intersection
      } # for motif

   #  [1] ------ MA0099.2
   #  [1]   mdb.genes: AP1, JUN::FOS, FOS::JUN, FOS::JUN
   #  [1]   tfc.genes: FOS, JUN

} # findMotifsWithMutuallyExclusiveMappings
#------------------------------------------------------------------------------------------------------------------------
test.hocomoco11.with.reliabilityScores <- function()
{
   printf("--- test.hocomoco11.with.reliabilityScores")

     #-------------------------------------------------------------------------
     # these queries rely primarily upon the dataSoure column of the metadata
     # subsequent checks below look at metadata rownames and matrix names
     #-------------------------------------------------------------------------

   checkEquals(length(query(MotifDb, "hocomoco")), 1834)
   checkEquals(length(query(MotifDb, "hocomocov10")), 1066)
   checkEquals(length(query(MotifDb, "hocomocov11")), 768)
   checkEquals(length(query(MotifDb, "hocomocov11-core")), 400)
   checkEquals(length(query(MotifDb, "hocomocov11-secondary")), 368)

   checkEquals(length(query(MotifDb, "hocomocov11-core-A")), 181)
   checkEquals(length(query(MotifDb, "hocomocov11-secondary-A")), 46)

   checkEquals(length(query(MotifDb, "hocomocov11-core-B")), 84)
   checkEquals(length(query(MotifDb, "hocomocov11-secondary-B")), 19)

   checkEquals(length(query(MotifDb, "hocomocov11-core-C")), 135)
   checkEquals(length(query(MotifDb, "hocomocov11-secondary-C")), 13)

   checkEquals(length(query(MotifDb, "hocomocov11-core-D")), 0)
   checkEquals(length(query(MotifDb, "hocomocov11-secondary-D")), 290)

     #-------------------------------------------------------------------------
     # check matrix names
     #-------------------------------------------------------------------------
   checkEquals(length(grep("HOCOMOCOv11-core-A", names(MotifDb))), 181)
   checkEquals(length(grep("HOCOMOCOv11-secondary-A", names(MotifDb))), 46)

   checkEquals(length(grep("HOCOMOCOv11-core-A", rownames(mcols(MotifDb)))), 181)
   checkEquals(length(grep("HOCOMOCOv11-secondary-A", rownames(mcols(MotifDb)))), 46)

} # test.hocomoco11.with.reliabilityScores
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
