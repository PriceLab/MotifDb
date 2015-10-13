# uniprobe/test.R
#------------------------------------------------------------------------------------------------------------------------
library(RUnit)
source("import.R")
#------------------------------------------------------------------------------------------------------------------------
run.tests = function (dataDir)
{
  stopifnot(!missing(dataDir))
  dataDir <- file.path(dataDir, "uniprobe")
  stopifnot(file.exists(dataDir))  
  
  txxa <- test.createMatrixNameUniqifier ()
  txx1 <- test.extractNativeNames (dataDir)
  txx2 <- test.createPublicationRefTable ()
  txxb <- test.standardizeSpeciesNames ()
  txx4 <- test.createGeneRefTable(dataDir)  # read mysql database to get, uniprobe-style, pwm-specific metadata
  txx3 <- test.uniprotToStandardID ()
  txx5 <- test.getAllStandardIDs (dataDir)   # build a 'stdID' column to add to tbl.geneRef
  txx5a <- test.parsePWMfromText (dataDir)

  txx6 <- test.extractPWMfromFile (dataDir)
  txx7 <- test.translateFileNameToGeneName ()
  txx8 <- test.readAndParse (dataDir)
  txx9 <- test.createMetadata (dataDir)
  txxb <- test.createMetadata.sox4 (dataDir)
  txxb <- test.createMetadata.cbf1 (dataDir)
  txxba <- test.createMetadata.fixTrailingSpaceInBindingDomain (dataDir)
  txxc <- test.renameMatrices ()
  txxd <- test.emptyStringProteinId (dataDir)
  
} # run.tests
#------------------------------------------------------------------------------------------------------------------------
# what tables?  what columns in each?
exploreDatabase = function ()
{
  if (!exists ('db'))
    db <<- dbConnect (MySQL (), dbname='uniprobe')
  tables = dbListTables (db)
  excluders = grep ('mers_', tables)
  if (length (excluders) > 0)
    tables = tables [-excluders]
  for (table in tables) {
    fields = dbListFields (db, table)
    printf ('--- %s: %d', table, length (fields))
    printf ('  %s', list.to.string (fields)) 
    } # for table


} # exploreDatabase
#------------------------------------------------------------------------------------------------------------------------
test.parsePWMfromText = function (dataDir)
{
  print ('--- test.parsePWMfromText')

  sample.file = file.path(dataDir, 'All_PWMs/Cell08/Phox2a_3947.1_pwm.txt')
  lines.of.text = scan (sample.file, what=character(0), sep='\n', skip=2, quiet=TRUE)
  yy <<- lines.of.text
  checkEquals (grep ('^A:', lines.of.text), 1)
  pwm.phox2a = parsePWMfromText (lines.of.text)
  checkTrue (is.matrix (pwm.phox2a))
  checkEquals (dim (pwm.phox2a), c (4, 16))
  xx <<- pwm.phox2a
  for (c in 1:ncol (pwm.phox2a)) {
    checkEqualsNumeric (sum (pwm.phox2a [, c]), 1.0, tol=0.01)
    } # for c
  
  invisible (pwm.phox2a)

} # test.parsePWMfromText
#------------------------------------------------------------------------------------------------------------------------
test.extractPWMfromFile = function (dataDir)
{
  print ('--- test.extractPWMfromFile')

    # a well-behaved file, with a single unaccompanied tab separating columns
  file = file.path(dataDir, 'All_PWMs/SCI09/Arid3a_pwm_primary.txt')
  checkTrue(file.exists(file))
  
  arid2a.matrix <<- extractPWMfromFile (file)
  checkEquals (rownames (arid2a.matrix), c ('A','C', 'G', 'T'))
  checkTrue (all (as.integer (round (as.double (colSums (arid2a.matrix)))) == 1))

    # a file with a confusion of white spaces
    # "A:\t  0.247560633305795  0.264091385393575\t\t 0.204731944407109 ....
    # depends on perl=T whitespace+ regex in parsePWMfromText

  file = file.path(dataDir, 'All_PWMs/GD09/Nsy-7.pwm')
  nsy7.matrix <<- extractPWMfromFile (file)
  checkEquals (rownames (nsy7.matrix), c ('A','C', 'G', 'T'))
  checkTrue (all (as.integer (round (as.double (colSums (nsy7.matrix)))) == 1))

} # test.extractPWMfromFile
#------------------------------------------------------------------------------------------------------------------------
test.createPublicationRefTable = function ()
{
  print ('--- test.createPublicationRefTable')
  tbl.ref = createPublicationRefTable ()
  checkEquals (dim (tbl.ref), c (10, 6))
  checkEquals (colnames (tbl.ref), c ('folder', 'author', 'pmid', 'organism', 'count', 'title'))

     # make sure neither author nor pmid repeat
  #checkEquals (length (unique (tbl.ref$author)), nrow (tbl.ref))
  #checkEquals (length (unique (tbl.ref$pmid)), nrow (tbl.ref))

  invisible (tbl.ref)

} # test.createPublicationRefTable
#------------------------------------------------------------------------------------------------------------------------
test.identifyFiles = function ()
{
  print ('--- test.identifyFiles')

    # just one legit pwm file in GD09:

  files.found = identifyFiles ('All_PWMs/GD09')
  checkEquals (length (files.found), 1)


  files.found = identifyFiles ('.')

  checkTrue (length (files.found) > 375)

  invisible (files.found)


} # test.identifyFiles
#------------------------------------------------------------------------------------------------------------------------
test.readAndParse = function (dataDir)
{
  print ('--- test.readAndParse')
  all.files = identifyFiles (file.path(dataDir, 'All_PWMs'))
  test.file = grep ('GD09', all.files, v=TRUE)
  checkEquals (length (test.file), 1)
  mtx = readAndParse (test.file)
  checkEquals (length (mtx), 1)
  checkTrue (is.matrix (mtx [[1]]))

  random.sample <- c (98, 99, 251)
  test.files.3 = all.files [random.sample]
  mtx3 = readAndParse (test.files.3)
  checkEquals (length (mtx3), 3)
  expected.matrix.names <- c("PNAS08/PF14_0633.pwm", "PNAS08/PFF0200c.pwm",
                             "Cell08/Pou4f3_2791.1_pwm.txt")
  checkEquals(names(mtx3), expected.matrix.names)

  invisible (mtx3)

} # test.readAndParse
#------------------------------------------------------------------------------------------------------------------------
test.translateFileNameToGeneName = function ()
{
  print ('--- test.translateFileNameToGeneName')
  checkEquals (translateFileNameToGeneName ('./All_PWMs/GD09/Nsy-7.pwm'), 'Nsy-7')
  checkEquals (translateFileNameToGeneName ('Nsy-7.pwm'), 'Nsy-7')
  checkEquals (translateFileNameToGeneName ('Nsy-7'), 'Nsy-7')

} # test.translateFileNameToGeneName
#------------------------------------------------------------------------------------------------------------------------
test.createMetadata = function (dataDir)
{
  print ('--- test.createMetadata')
  all.files = identifyFiles (file.path(dataDir,'All_PWMs'))

  sample.files = all.files [c (12, 36, 259, 269, 309)]   #   sample (1:length (all.files), 5)]
    # add in all entries from the NBT paper, which are a mix of organisms
    # Scerevisiae;Hsapiens;Mmusculus.  want to make sure that organism is drawn from tbl.geneRef, rather than
    # the old approach, which learned matrix/gene/organism relationships from tbl.pubRef
  sample.files = c (sample.files, grep ('NBT', all.files, v=T))

  matrices.10 = readAndParse (sample.files)
  tbl.pubRef = createPublicationRefTable ()
  tbl.geneRef = createGeneRefTable (dataDir)
  tbl.md = createMetadata (matrices.10, tbl.pubRef, tbl.geneRef)

  checkEquals (dim (tbl.md), c (10, 15))
  checkEquals (tbl.md$geneId, c ("YPR104C", "YBR089C-A", "194738", "21410", "54131", "YJR060W", "179485", "5451", "YNL216W", "13653"))

  checkEquals (tbl.md$geneIdType, c ("SGD", "SGD", "ENTREZ", "ENTREZ", "ENTREZ", "SGD", "ENTREZ", "ENTREZ", "SGD", "ENTREZ"))

  checkEquals (tbl.md$organism, c ("Scerevisiae", "Scerevisiae", "Mmusculus", "Mmusculus", "Mmusculus", "Scerevisiae", 
                                   "Celegans", "Hsapiens", "Scerevisiae", "Mmusculus"))
  checkEquals (tbl.md$geneSymbol, c ("Fhl1", "Nhp6b", "Rhox11", "Tcf2", "Irf3", "Cbf1", "Ceh-22", "Oct-1", "Rap1", "Zif268"))
  
  checkEquals (tbl.md$proteinId,     c ("P39521", "P11633", "Q810N8", "P27889", "P70671", "P17106", "P41936", "P14859", "P11938", "P08046"))
  checkEquals (tbl.md$proteinIdType, c ("UNIPROT", "UNIPROT", "UNIPROT", "UNIPROT", "UNIPROT", "UNIPROT", "UNIPROT", "UNIPROT", "UNIPROT", "UNIPROT"))
  checkEquals (tbl.md$bindingDomain, c ("Fork-head", "HMG_box", "Homeo", "Homeo", "IRF", "HLH", "Homeo", "POU", "Myb", "ZnF_C2H2"))
  checkEquals (tbl.md$bindingSequence, c (NA_character_, NA_character_,
                                         'PRKAYRFTPGQLWELQAVFVENQYPDALKRKELAGLLNVDEQKIKDWFNNKRAKYRK',
                                         'RRNRFKWGPASQQILYQAYDRQKNPSKEEREALVEECNRAECLQRGVSPSKAHGLGSNLVTEVRVYNWFANRRKEEAF',
                                         NA_character_,
                                         NA_character_,
                                         NA_character_,
                                         NA_character_,
                                         NA_character_,
                                         NA_character_))


  checkEquals (length (which (duplicated (tbl.md$providerName))), 0)

  full.names = rownames (tbl.md)
  checkEquals (length (grep ('.UP0', full.names)), nrow (tbl.md))
  checkEquals (length (grep ('-UniPROBE', full.names)), nrow (tbl.md))

     # should be 10 full.names, each of 3 tokens: "Scerevisiae-UniPROBE-Fhl1.UP00303"
  tokens = strsplit (full.names, '-')
  checkEquals (length (tokens), nrow (tbl.md))
  all (sapply (tokens, function (sub.tokens) checkEquals (length (sub.tokens), 3)))

  invisible (tbl.md)

} # test.createMetadata
#------------------------------------------------------------------------------------------------------------------------
# uniprobe reports yeast Cbf1 from two studies. one is from NBT06 which is a multiple species pwm collection.  this
# will be handled properly in the future (TODO) but not yet.  the more common UniPROBE practive of one-organism/one-paper
# is used throughout to assign organism to pwm.  
# createMetadata can return too many rows in some cases, mostly (always) 
#
test.createMetadata.cbf1 = function (dataDir)
{
  print ('--- test.createMetadata.cbf1')
  all.files = identifyFiles (file.path(dataDir, 'All_PWMs'))
  cbf1.files = all.files [grep ('cbf1', all.files, ignore.case=T)]
  stopifnot (length (cbf1.files) == 2)   # one in GR09, one in NBT06

  matrices.cbf1 = readAndParse (cbf1.files)
  tbl.pubRef = createPublicationRefTable ()
  tbl.geneRef = createGeneRefTable (dataDir)
  tbl.md = createMetadata (matrices.cbf1, tbl.pubRef, tbl.geneRef)
  checkEquals (dim (tbl.md), c (2, 15))
  checkEquals (tbl.md$geneSymbol, c ('Cbf1', 'Cbf1'))
  checkEquals (tbl.md$providerName, c ('GR09/Cbf1.pwm', 'NBT06/Cbf1.pwm'))

  invisible (tbl.md)

} # test.createMetadata.cbf1
#------------------------------------------------------------------------------------------------------------------------
# some of the metadata for matrices in Cell08, obtained through sql in createGeneRefTable,  ended up with a
# trailing blank ('Homeo ') in the bindingDomain field.  make sure this is fixed.
# currently the fix is to define and use a string 'trim' method in createMetaData ()
test.createMetadata.fixTrailingSpaceInBindingDomain = function (dataDir)
{
  print ('--- test.createMetadata.fixTrailingSpaceInBindingDomain')
  all.files = identifyFiles (file.path(dataDir,'All_PWMs'))
  sample.trailingSpaceCulpritFiles = head (grep ('Irx6', all.files))
  stopifnot (length (sample.trailingSpaceCulpritFiles) == 1) 
  filenames = all.files [sample.trailingSpaceCulpritFiles]
  matrices.tmp = readAndParse (filenames)
  tbl.pubRef = createPublicationRefTable ()
  tbl.geneRef = createGeneRefTable (dataDir)
    # the fix is in here.  
  tbl.md = createMetadata (matrices.tmp, tbl.pubRef, tbl.geneRef)
  checkEquals (dim (tbl.md), c (1, 15))
  checkEquals (tbl.md$bindingDomain, "Homeo")

} # test.createMetadata.fixTrailingSpaceInBindingDomain
#------------------------------------------------------------------------------------------------------------------------
# uniprobe reports sox4 in both mouse and human.  these different organisms must be recognized for this to work.
# it used to fail when organism was variously (e.g.) 'Homo sapiens' and 'Hsapiens'.
# the function standardizeSpeciesNames was added to fix this, converting everything in tbl.geneRef$species
#
test.createMetadata.sox4 = function (dataDir)
{
  print ('--- test.createMetadata.sox4')
  all.files = identifyFiles (file.path(dataDir, 'All_PWMs'))
  sox4.files = all.files [grep ('sox4', all.files, ignore.case=T)]
  stopifnot (length (sox4.files) == 2)   # one in CR09, one in SCI09

  matrices.sox4 = readAndParse (sox4.files)
  tbl.pubRef = createPublicationRefTable ()
  tbl.geneRef = createGeneRefTable (dataDir)
  tbl.md = createMetadata (matrices.sox4, tbl.pubRef, tbl.geneRef)
  checkEquals (dim (tbl.md), c (2, 15))

} # test.createMetadata.sox4
#------------------------------------------------------------------------------------------------------------------------
test.extractNativeNames = function (dataDir)
{
  print ('--- test.extractNativeNames')

     # check a random sample, extract called one raw name at at time

  checkEquals (extractNativeNames ("Cell08/Lhx8_2247.2_pwm.txt"), 'Lhx8')
  checkEquals (extractNativeNames ("Cell08/Phox2a_3947.1_pwm.txt"), 'Phox2a')

  checkEquals (extractNativeNames ("GR09/Srd1.pwm"), 'Srd1')
  checkEquals (extractNativeNames ("Cell08/Pou3f3_3235.2_pwm.txt"), 'Pou3f3')
  checkEquals (extractNativeNames ("Cell08/Msx3_3206.1_pwm.txt"), 'Msx3')

    # now do the sample again in one call

  names.5 = extractNativeNames (c ("Cell08/Lhx8_2247.2_pwm.txt",
				   "Cell08/Phox2a_3947.1_pwm.txt",
				   "GR09/Srd1.pwm",
				   "Cell08/Pou3f3_3235.2_pwm.txt",
				   "Cell08/Msx3_3206.1_pwm.txt"))
  checkEquals (names.5, c ("Lhx8", "Phox2a", "Srd1", "Pou3f3", "Msx3"))

  all.files = identifyFiles (file.path(dataDir,'All_PWMs'))
  matrices.all = readAndParse (all.files)
  all.raw.names = names (matrices.all)
  cooked.names.all = c ()

  for (raw.name in all.raw.names) {
    cooked.name = extractNativeNames (raw.name)
    cooked.names.all = c (cooked.names.all, cooked.name)
    #printf ('%30s: %10s', raw.name, cooked.name)
    } # for raw.name
  
  cooked.names.all.v2 = extractNativeNames (all.raw.names)
  checkEquals (cooked.names.all, cooked.names.all.v2)

  invisible (cooked.names.all)

} # test.extractNativeNames
#------------------------------------------------------------------------------------------------------------------------
test.uniprotToStandardID = function ()
{
  print ('--- test.uniprotToStandardID')
  checkEquals (uniprotToStandardID ('mouse', 'Q14A64')$std, '11298')
  checkEquals (uniprotToStandardID ('human', 'P14859')$std, '5451')
  checkEquals (uniprotToStandardID ('yeast', 'P17106')$std, "YJR060W")
  checkEquals (uniprotToStandardID ('worm',  'P41936')$std, "179485")

  yeast.uids = c ("P22149", "Q04052", "P40467", "P22035", "P17106")
  checkEquals (uniprotToStandardID ('yeast', yeast.uids)$std,
               c ("YGL071W", "YDR421W", "YIL130W", "YKR099W", "YJR060W"))

  checkTrue (is.na (uniprotToStandardID ('mouse', 'Q8BI45')$std))
  mouse.uids = c ("O70137", "O35137", "Q62431", "Q8BI45", "O35085", "O35885")
  result = uniprotToStandardID ('mouse', mouse.uids)
  checkEquals (result$std, c ("11694", "11695", "13496", NA, "11878", "17173"))

  worm.uniprots =  c ("P41936", "Q4PIU8")
  result = uniprotToStandardID ('worm', worm.uniprots)
  checkEquals (result$uniprot, worm.uniprots)
  checkEquals (result$std, c ('179485', NA))

     # expose and solve a problem with duplicates
     # "P02831"  "P97436" "Q9QZ28" "Q66JY7" "P08046" 
  mouse.uids = c ("P02831", "P02831", '', '')
  result = uniprotToStandardID ('mouse', mouse.uids)
  checkEquals (result$uniprot, mouse.uids)
  checkEquals (result$std, c ('15400', '15400', NA, NA))

     # mouse Zic[123] maps to geneID 22771-3.  but 22772 is not found for Zic2.  why?
  checkEquals (uniprotToStandardID ('mouse', 'P46684')$std, '22771')
  #checkEquals (uniprotToStandardID ('mouse', 'Q62520')$std, '22772')    TODO!
  checkEquals (uniprotToStandardID ('mouse', 'Q62521')$std, '22773')

} # test.uniprotToStandardID
#------------------------------------------------------------------------------------------------------------------------
test.createGeneRefTable = function (dataDir)
{ 
  print ('--- test.createGeneRefTable')
  tbl.geneRef = createGeneRefTable (dataDir)  
  checkEquals (length (grep (' ', tbl.geneRef$species)), 0)   # no spaces in the species names
  checkEquals (ncol (tbl.geneRef), 13)
  checkEquals (sort (colnames (tbl.geneRef)),
               c ("bindingSequence", "domain", "folder_name", "id", "ihop", "name", "pmid", "pub", "refseq_id", "species",  
                   "stdID", "type", "uniprot"))
  checkEquals (nrow (tbl.geneRef), 372)
  checkEquals (length (which (duplicated (tbl.geneRef [, 1:2]))), 0)
  checkEquals (length (which (!is.na (tbl.geneRef$bindingSequence))), 168)
  
  invisible (tbl.geneRef)

} # test.createGeneRefTable
#------------------------------------------------------------------------------------------------------------------------
test.getAllStandardIDs = function (dataDir)
{
  print ('--- test.getAllStandardIDs')
  if (!exists ('tbl.geneRef'))
     tbl.geneRef =  createGeneRefTable (dataDir)  

  stdID = getAllStandardIDs (tbl.geneRef)
  checkEquals (length (stdID), nrow (tbl.geneRef))
  tbl.g2 =  cbind (tbl.geneRef, stdID, stringsAsFactors=FALSE)
  checkEquals (nrow (tbl.g2), nrow (tbl.geneRef))

     # every row in the table for a yeast pwm should have stdID of NA or start with Y
  should.be.orf.names = subset (tbl.g2, species=='Scerevisiae')$stdID
  checkEquals ((length (which (is.na (should.be.orf.names))) +
                length (grep ('^Y', should.be.orf.names))),
               length (should.be.orf.names))
    # no other rows should have stdID starting with Y
  should.not.be.orf.names = subset (tbl.g2, species!='Scerevisiae')$stdID
  checkEquals (length (grep ('^Y', should.not.be.orf.names)), 0)

  checkTrue (length (which (is.na (tbl.g2$stdID))) < 70)

 	    # now handcheck a few from each species
  checkEquals (subset (tbl.g2, name=='Cbf1')$stdID, c ("YJR060W", "YJR060W"))
  checkEquals (subset (tbl.g2, name=='Yox1')$stdID, "YML027W")

  invisible (tbl.g2)

} # test.getAllStandardIDs
#------------------------------------------------------------------------------------------------------------------------
test.standardizeSpeciesNames = function ()
{
  print ('--- test.standardizeSpeciesNames')
  old.style = c ("Saccharomyces cerevisiae", "Mus musculus", "Caenorhabditis elegans", "Cryptosporidium parvum", 
                 "Vibrio harveyi", "Homo sapiens", "Plasmodium falciparum")
  new.style = standardizeSpeciesNames (old.style [1])
  checkEquals (new.style, 'Scerevisiae')

  new.style = standardizeSpeciesNames (old.style)
  checkEquals (new.style, c ("Scerevisiae", "Mmusculus", "Celegans", "Cparvum", "Vharveyi", "Hsapiens", "Pfalciparum"))

} # test.standardizeSpeciesNames
#------------------------------------------------------------------------------------------------------------------------
test.renameMatrices = function (matrices, tbl.md, tbl.anno)
{
  print ('--- test.renameMatrices')

} # test.renameMatrices
#------------------------------------------------------------------------------------------------------------------------
# Scerevisiae-UniPROBE-Gsm1.UP00301, GR09/Gsm1.pwm, UP00301, Gsm1, NA, NA, "", UNIPROT, NA,  Scerevisiae, Zn2Cys67
test.emptyStringProteinId = function (dataDir)
{
  print ('--- test.emptyStringProteinId')
  all.files = identifyFiles (file.path(dataDir, 'All_PWMs'))
  filename.gsm1 = grep ('Gsm1', all.files, v=T)
  mtx = readAndParse (filename.gsm1)
  tbl.pubRef = createPublicationRefTable ()
  tbl.geneRef = createGeneRefTable (dataDir)
  tbl.md = createMetadata (mtx, tbl.pubRef, tbl.geneRef)
  checkEquals (nrow (tbl.md), 1)
  md.list = as.list (tbl.md)
  checkEquals (md.list$proteinId, 'NP_012432')
  checkEquals (md.list$proteinIdType, 'REFSEQ')


} # test.emptyStringProteinId 
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
