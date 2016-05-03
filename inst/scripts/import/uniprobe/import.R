# uniprobe/import.R
#------------------------------------------------------------------------------------------------------------------------
#library (RMySQL)
library (RSQLite)
library (org.Hs.eg.db)
library (org.Mm.eg.db)
library (org.Sc.sgd.db)
library (org.Ce.eg.db)
library(tools)   # for md5sum
#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
run = function (dataDir)
{
  dataDir <- file.path(dataDir, "uniprobe")
  stopifnot(file.exists(dataDir))

  all.files = identifyFiles (file.path(dataDir, 'All_PWMs'))
  matrices = readAndParse (all.files)
  tbl.pubRef = createPublicationRefTable ()
  tbl.geneRef = createGeneRefTable (dataDir)
  tbl.md = createMetadata (matrices, tbl.pubRef, tbl.geneRef)
  stopifnot (length (matrices) == nrow (tbl.md))
  matrices = renameMatrices (matrices, tbl.md)

  serializedFile <- "uniprobe.RData"
  save (matrices, tbl.md, file=serializedFile)
  printf("saved %d matrices to %s", length(matrices), serializedFile)
  printf("next step: copy %s to <packageRoot>/MotifDb/inst/extdata, rebuild package", serializedFile)

} # run
#------------------------------------------------------------------------------------------------------------------------
createMatrixNameUniqifier = function (matrix)
{
  temporary.file <- tempfile ()
  write (as.character (matrix), file=temporary.file)
  md5sum.string <- as.character (md5sum (temporary.file))
  stopifnot (nchar (md5sum.string) == 32)
  md5.6chars = substr (md5sum.string, 29, 32)
  #unlink (temporary.file)

} # createMatrixNameUniqifier
#------------------------------------------------------------------------------------------------------------------------
parsePWMfromText = function (lines.of.text)
{
  z = lines.of.text
  stopifnot (sum (sapply (c ('A', 'C', 'T', 'G'), function (token) length (grep (token, lines.of.text)))) == 4)

    # determine the number of columns
  token.count.first.line = length (strsplit (lines.of.text [1], '\t')[[1]])
  column.count = token.count.first.line - 1  # subtract off the 'A:\t' token
  result = matrix (nrow=4, ncol=column.count, byrow=TRUE, dimnames=list (c ('A', 'C', 'G', 'T'), 1:column.count))

  for (line in lines.of.text) {
    tokens = strsplit (line, '\\s*[:\\|]') [[1]]
    nucleotide = tokens [1]
    numbers.raw = tokens [2]
    number.tokens = strsplit (numbers.raw, '\\s+', perl=T)[[1]]
    while (nchar (number.tokens [1]) == 0)
      number.tokens = number.tokens [-1]
    numbers = as.numeric (number.tokens)
    #printf ('adding %s: %s', nucleotide, list.to.string (numbers))
    result [nucleotide,] = numbers
    }

  return (result)

} # parsePWMfromText
#------------------------------------------------------------------------------------------------------------------------
extractPWMfromFile = function (filename)
{
  text = scan (filename, sep='\n', what=character (0), quiet=TRUE)
  matrix.start.lines = grep ('A\\s*[:\\|]', text)
  stopifnot (length (matrix.start.lines) == 1)
  #printf ('%50s: %s', filename, list.to.string (matrix.start.lines))
  start.line = matrix.start.lines [1]
  end.line = start.line + 3
  lines = text [start.line:end.line]
  pwm.matrix = parsePWMfromText (lines)
  return (pwm.matrix)

} # extractPWMfromFile
#------------------------------------------------------------------------------------------------------------------------
createPublicationRefTable = function ()
{
  options (stringsAsFactors=FALSE)
  tbl.ref = data.frame (folder=c('CR09', 'Cell08', 'Cell09', 'EMBO10', 'GD09', 'GR09', 'MMB08', 'NBT06', 'PNAS08', 'SCI09'))
  tbl.ref = cbind (tbl.ref, author=c('Scharer', 'Berger', 'Grove', 'Wei', 'Lesch', 'Zhu', 'Pompeani', 'Berger', 'De Silva', 'Badis'))
  tbl.ref = cbind (tbl.ref, pmid=c('19147588','18585359','19632181','20517297','19204119','19158363','18681939','16998473','18541913','19443739'))
  tbl.ref = cbind (tbl.ref, organism=c('Hsapiens', 'Mmusculus', 'Celegans', 'Mmusculus', 'Celegans', 'Scerevisiae', 'Vharveyi',
                                       'Scerevisiae;Hsapiens;Mmusculus', 'Apicomplexa', 'Mmusculus'))
  tbl.ref = cbind (tbl.ref, count=c(1, 168, 34, 25, 1, 89, 1, 5, 3, 104))
  titles = c ('Genome-wide promoter analysis of the SOX4 transcriptional network in prostate cancer cells',
              'Variation in homeodomain DNA binding revealed by high-resolution analysis of sequence preferences',
              'A Multiparameter Network Reveals Extensive Divergence between C. elegans bHLH Transcription Factors',
              'Genome-wide analysis of ETS-family DNA-binding in vitro and in vivo',
              'Transcriptional regulation and stabilization of left right neuronal identity in C. elegans',
              'High-resolution DNA-binding specificity analysis of yeast transcription factors',
              'The Vibrio harveyi master quorum-sensing regulator, LuxR...',
              'Compact, universal DNA microarrays...',
              'Specific DNA-binding by Apicomplexan AP2 transcription factors',
              'Diversity and complexity in DNA recognition by transcription factors')
  tbl.ref = cbind (tbl.ref, title=titles)

  tbl.ref

} # createPublicationRefTable
#------------------------------------------------------------------------------------------------------------------------
identifyFiles = function (filePath)
{
  stopifnot(file.exists(filePath))
  
  cmd = sprintf ('find %s -name "*pwm*"', filePath)

  files.raw = system (cmd, intern=TRUE)

      # all legit files end in ".pwm" or ".txt"
  files = c (grep (".pwm$", files.raw, value=T, ignore.case=T), 
             grep (".txt$", files.raw, value=T, ignore.case=T))

  reverseComplementMatrices = grep ('RC', files)
  if (length (reverseComplementMatrices) > 0)
    files = files [-reverseComplementMatrices]

      # don't know why these were excluded.  leave them in for now
  embo10.excluders = grep ('EMBO', files)
  if (length (embo10.excluders) > 0)
    files = files [-embo10.excluders]

  cell09.excluders = grep ('Cell09', files)   # include these once the sql tables are updated
  if (length (cell09.excluders) > 0)
    files = files [-cell09.excluders]
  secondary.excluders = grep ('secondary', files)
  if (length (secondary.excluders) > 0)
    files = files [-secondary.excluders]

  invisible (files)

} # identifyFiles
#------------------------------------------------------------------------------------------------------------------------
readAndParse = function (file.list)
{
  matrices = list ()

  for (file in file.list) {
    #printf ('read and parse %s', file)
    text = scan (file, sep='\n', what=character (0), quiet=TRUE)
    matrix.start.lines = grep ('A:', text)
    stopifnot (length (matrix.start.lines) == 1)
    start.line = matrix.start.lines [1]
    end.line = start.line + 3
    lines.of.text = text [start.line:end.line]
    pwm.matrix = parsePWMfromText (lines.of.text)
    name.tokens <- strsplit(file,"/")[[1]]
    token.count <- length(name.tokens)

    matrix.name <- paste(name.tokens[(token.count-1):token.count], collapse="/")
    matrices [[matrix.name]] = pwm.matrix
    }

  invisible (matrices)

} # readAndParse
#------------------------------------------------------------------------------------------------------------------------
# eg, ./All_PWMs/GD09/Nsy-7.pwm to simply 'Nsy-7'
#
translateFileNameToGeneName = function (long.name)
{
  if (length (grep ('/', long.name)) > 0) {
    tokens = strsplit (long.name, '/')[[1]]
    count = length (tokens)
    gene.name.raw = tokens [count]  # get the last one
    }
  else {
    gene.name.raw = long.name
    }

  gene.name = strsplit (gene.name.raw, '\\.')[[1]][1]   # remove any file suffix

  gene.name

} # test.translateFileNameToGeneName
#------------------------------------------------------------------------------------------------------------------------
# and update matrix names
createMetadata = function (matrices, tbl.pubRef, tbl.geneRef)
{
  options (stringsAsFactors=FALSE)
  dataSource = 'UniPROBE'
  trim = function (s) {sub (' +$', '', sub ('^ +', '', s))}
  tbl.md = data.frame ()
  
  removers = list (Cart1=110, Cutl1=115, Hoxa7=156, Irx3=185)
 
  for (m in 1: length (matrices)) {
    matrix.name = names (matrices) [m]
    #printf ('%d: %s', m, matrix.name)
    native.name.raw = gsub ('All_PWMs/', '', matrix.name)
    native.name = extractNativeNames (native.name.raw)

       # extractNativeNames needs to be rethought.  but for now, just hack in some fixes
    if (native.name == 'Cgd2') native.name='Cgd2_3490'   # inconsistency at uniprobe, accomodated here
    if (native.name == 'PF14') native.name='PF14_0633'
    if (native.name == 'Uncx4') native.name='Uncx4.1'

    if (!native.name %in% tbl.geneRef$name) {
      browser (text='native.name not in tbl.geneRef')
      }
    
    experiment.folder = strsplit (native.name.raw, '/')[[1]][1]
    up.id.number = subset (tbl.geneRef, name==native.name & folder_name==experiment.folder)$id
    # todo BUG here!
    full.uniprobe.id = sprintf ('UP%05d', up.id.number)
    if (length (full.uniprobe.id) != 1) {
      browser (text='full.uniprobe.id has wrong length')
      }

    geneId = subset (tbl.geneRef, name==native.name & folder_name==experiment.folder)$stdID
    if (is.na (geneId) | geneId == '')
       geneId = NA_character_

    bindingDomain = trim (subset (tbl.geneRef, name==native.name & folder_name==experiment.folder)$domain)
    organism = subset (tbl.geneRef, name==native.name & folder_name==experiment.folder)$species

      # a native name (a gene symbol) may have a dash, but for the long name, we want dashes to separate
      # the organism-dataSource-geneIdentifier matrix name.
      # so covert this here, but do not eclipse any dash-including native.names

    native.name.no.dashes = gsub ('-', '_', native.name)
    native.name.uniq = sprintf ('%s-%s-%s.%s', organism, dataSource, native.name.no.dashes, full.uniprobe.id)
    if (native.name.uniq %in% rownames (tbl.md)) {
      matrix = matrices [[m]]
      uniqifier = createMatrixNameUniqifier (matrix)
      #printf ('before, not unique: %s', native.name.uniq)
      native.name.uniq = paste (native.name.uniq, uniqifier, sep='.')
      #printf ('after, not unique: %s', native.name.uniq)
      }
      
    sequenceCount = NA_integer_
    bindingSequence = NA_character_
    tfFamily = NA_character_

    experimentType = 'protein binding microarray'
    pubmedID = subset (tbl.pubRef, folder==experiment.folder)$pmid
    #printf ('%12s: %20s', native.name, organism)

    if (is.na (geneId))
      geneIdType = NA
    else if (organism == 'Scerevisiae')
      geneIdType = 'SGD'
    else if (organism %in% c ('Mmusculus', 'Hsapiens', 'Celegans'))
      geneIdType = 'ENTREZ'
    else
      geneIdType = 'todo'

    proteinId = NA_character_
    proteinIdType = NA_character_

    proteinId.tmp = subset (tbl.geneRef, name==native.name & folder_name==experiment.folder)$uniprot
    if (!is.na (proteinId.tmp) & nchar (proteinId.tmp) > 1) {
      #printf ('getting uniprot id for %s: %s', native.name, proteinId.tmp)
      proteinId = proteinId.tmp
      proteinIdType = 'UNIPROT'
      }  # found good id in the uniprot column

    if (is.na (proteinId.tmp) | nchar (proteinId.tmp) == 0) {
      proteinId.tmp = subset (tbl.geneRef, name==native.name & folder_name==experiment.folder)$refseq_id
      if (!is.na (proteinId.tmp) & nchar (proteinId.tmp) > 1) {
        #printf ('getting refseq id for %s: %s', native.name, proteinId.tmp)
        proteinId = proteinId.tmp
        proteinIdType = 'REFSEQ'
        }
      } # need to look in the refseq column


    bindingSequence =  subset (tbl.geneRef, name==native.name & folder_name==experiment.folder)$bindingSequence
    #printf ('%12s: %40s', native.name, bindingSequence)

    new.row = list (providerName=native.name.raw,
                    providerId=full.uniprobe.id,
                    dataSource=dataSource,
                    geneSymbol=native.name,
                    geneId=geneId,
                    geneIdType=geneIdType,
                    proteinId=proteinId,
                    proteinIdType=proteinIdType,
                    organism=organism,
                    sequenceCount=NA,
                    bindingSequence=bindingSequence,
                    bindingDomain=bindingDomain,
                    tfFamily=tfFamily,
                    experimentType=experimentType,
                    pubmedID=pubmedID)
    if (native.name.uniq %in% rownames (tbl.md)) browser (text='dup row name')
    tbl.md = rbind (tbl.md, data.frame (new.row, stringsAsFactors=FALSE))
    if (length (native.name.uniq) != 1) browser (text='native.name.unique length != 1')
    rownames (tbl.md) [m] = native.name.uniq
    }

  invisible (tbl.md)

} # createMetadata
#------------------------------------------------------------------------------------------------------------------------
extractNativeNames = function (native.names.raw)
{
  name.count = length (native.names.raw)
  result = vector ('character', name.count)

  for (i in 1:name.count) {
    native.name.raw = native.names.raw [i]
    tokens = strsplit (native.name.raw, '/') [[1]]
    count = length (tokens)
    cooked.1 = native.name.raw
    if (count > 1)
      cooked.1 = tokens [length (tokens)]
    tokens = strsplit (cooked.1, '[_\\.]')[[1]]
    cooked.2 = tokens [1]
    result [i] = cooked.2
    }

  invisible (result)

} # extractNativeNames
#------------------------------------------------------------------------------------------------------------------------
# 'standard' is usually entrez gene ID.  for yeast it is orf. for fly ...
uniprotToStandardID = function (organism, uniprot.ids)
{
#  uniprot.ids = unique (uniprot.ids)

  if (!exists ('lst.yeast.uniprot')) {
    tbl.tmp = toTable (org.Sc.sgdUNIPROT)
    yeast.uniprot <<- tbl.tmp$systematic_name
    names (yeast.uniprot) <<- tbl.tmp$uniprot_id
    }

  if (!exists ('mouse.uniprot')) {
    tbl.tmp = toTable (org.Mm.egUNIPROT)
    mouse.uniprot <<- tbl.tmp$gene_id
    names (mouse.uniprot) <<- tbl.tmp$uniprot_id
    }

  if (!exists ('human.uniprot')) {
    tbl.tmp = toTable (org.Hs.egUNIPROT)
    human.uniprot <<- tbl.tmp$gene_id
    names (human.uniprot) <<- tbl.tmp$uniprot_id
    }

  if (!exists ('worm.uniprot')) {
    tbl.tmp = toTable (org.Ce.egUNIPROT)
    worm.uniprot <<- tbl.tmp$gene_id
    names (worm.uniprot) <<- tbl.tmp$uniprot_id
    }

  organism = unique (organism)
  stopifnot (length (unique (organism)) == 1)
  stopifnot (organism %in% (c ('human', 'mouse', 'yeast', 'worm')))
 
  if (organism == 'human') {
    ids = human.uniprot [uniprot.ids]
    }
  else if (organism == 'mouse') {
    ids = mouse.uniprot [uniprot.ids]
    }
  else if (organism == 'yeast') {
    ids = yeast.uniprot [uniprot.ids]
    }
  else if (organism == 'worm') {
    ids = worm.uniprot [uniprot.ids]
    }


    # embarrassing but true:  could not get this to work by creating a list directly
    # round-about solution:  make a 1-column data.frame, then construct a named list from that
  df = data.frame (uniprot=uniprot.ids, std=as.character (ids), stringsAsFactors=FALSE)
  #rownames (df) = uniprot.ids
  #result = df$stdID
  #names (result) = uniprot.ids
  return (df)

} # uniprotToStandardID
#------------------------------------------------------------------------------------------------------------------------
# see ~/v/snotes/log "* load uniprobe sql dump file (11 may 2012)" for info on the uniprobe mysql database 
# used here.  
createGeneRefTable = function (dataDir)
{
  if (!exists ('db')){
      dbFile <- file.path(dataDir, "uniprobe.sqlite")
      stopifnot(file.exists(dbFile))
      db <<- dbConnect (dbDriver("SQLite"), dbFile)
      }

  tbl.pubmed = data.frame (folder=c('CR09', 'Cell08', 'Cell09', 'EMBO10', 'GD09', 'GR09', 'MMB08', 'NBT06', 'PNAS08', 'SCI09'),
                 pmid=c('19147588', '18585359', '19632181', '20517297', '19204119', '19158363', '18681939', '16998473', '18541913', '19443739'),
                 stringsAsFactors=FALSE)
       #CR09	1	Scharer 	19147588	human	Genome-wide promoter analysis of the SOX4 transcriptional 
       #                                                        network in prostate cancer cells
       #Cell08	168	Berger  	18585359	mouse	Variation in homeodomain DNA binding revealed by high-resolution 
       #                                                        analysis of sequence preferences
       #Cell09	34	Grove   	19632181	worm	A Multiparameter Network Reveals Extensive Divergence between 
       #                                                        C. elegans bHLH Transcription Factors
       #EMBO10	25	Wei     	20517297	mouse	Genome-wide analysis of ETS-family DNA-binding in vitro and in vivo
       #GD09	1	Lesch   	19204119	worm	Transcriptional regulation and stabilization of left right neuronal 
       #                                                         identity in C. elegans
       #GR09	89	Zhu     	19158363	yeast	High-resolution DNA-binding specificity analysis of yeast transcription factors
       #MMB08	1	Pompeani	18681939	Vibrio	The Vibrio harveyi master quorum-sensing regulator, LuxR...
       #NBT06	5	Berger  	16998473	yeast;human;mouse	Compact, universal DNA microarrays...
       #PNAS08	3	De Silva	18541913	apicomplexa	Specific DNA-binding by Apicomplexan AP2 transcription factors
       #SCI09	104	Badis   	19443739	mouse	Diversity and complexity in DNA recognition by transcription factors


  tbl.gene = dbGetQuery (db, 'select * from gene_ids_public')
  colnames (tbl.gene) = c ('id', 'name', 'species', 'pub', 'type')
  tbl.genomic = dbGetQuery (db, 'select * from genomic_info') [, -c(2,3,8:11)]  # remove the longish 'description' field
  tbl.pub = dbGetQuery (db, 'select * from publication_ids') [,-3]  # remove 'full_ref' for legibility
  tbl = merge (tbl.gene, tbl.pub [, c (1,4)], by.x='pub', by.y='publication_id', all.x=TRUE)
  tbl = merge (tbl, tbl.pubmed, by.x='folder_name', by.y='folder', all.x=TRUE)
  tbl = merge (tbl, tbl.genomic, all.x=TRUE, by.x='name', by.y='gene_name')
  redundant.column = grep ('species.y', colnames (tbl))
  stopifnot (length (redundant.column) == 1)
  tbl = tbl [, -redundant.column]
  column.to.rename = grep ('species.x', colnames (tbl))
  stopifnot (length (column.to.rename) == 1)
  colnames (tbl) [column.to.rename] = 'species'
  

   # now add 'standard IDs' -- orf names for yeast, entrez geneIDs for everything else

  stdID = getAllStandardIDs (tbl)
  tbl = cbind (tbl, stdID)
  duplicates = which (duplicated (tbl [, c (1,2)]))
  if (length (duplicates) > 0)
    tbl = tbl [-duplicates,]
  
    # fix the organism (species) name, from (e.g.) "Homo sapiens" to "Hsapiens"
  tbl$species = standardizeSpeciesNames (tbl$species)
  invisible (tbl)

    # now add bindingSequence, from DBDs:  gene_name + seq, apparently the binding sequence when known, of 171 of 372 genes
  tbl.dbds = dbGetQuery (db, 'select * from DBDs')  
  dbds.list = tbl.dbds$seq
  names (dbds.list) = tbl.dbds$gene_name
  bindingSequence = rep (NA_character_, nrow (tbl))
  names (bindingSequence) = tbl$name
  shared.names = unique (intersect (tbl.dbds$gene_name, tbl$name))
  bindingSequence [shared.names] = dbds.list [shared.names]
  tbl = cbind (tbl, bindingSequence)
  invisible (tbl)
  
} # createGeneRefTable
#------------------------------------------------------------------------------------------------------------------------
# make successive species-specific calls to 'uniprotToStandardID', assembling a new column to be added to tbl.geneRef
getAllStandardIDs = function (tbl.geneRef)
{
  mouse.rows  = grep ('Mus musculus', tbl.geneRef$species)
  yeast.rows = grep ('Saccharomyces cerevisiae', tbl.geneRef$species)
  human.rows = grep ('Homo sapiens', tbl.geneRef$species)
  worm.rows = grep ('Caenorhabditis elegans', tbl.geneRef$species)

  stdID = rep (NA_character_, nrow (tbl.geneRef))

  mouse.uids = tbl.geneRef$uniprot [mouse.rows]
  yeast.uids = tbl.geneRef$uniprot [yeast.rows]
  human.uids = tbl.geneRef$uniprot [human.rows]
  worm.uids = tbl.geneRef$uniprot [worm.rows]

    # add mouse geneIDs
  tbl.mouse =  uniprotToStandardID ('mouse', mouse.uids)
  stopifnot (length (mouse.rows) == nrow (tbl.mouse))
  stdID [mouse.rows] = tbl.mouse$std
  
    # add yeast orfs
  tbl.yeast =  uniprotToStandardID ('yeast', yeast.uids)
  stopifnot (length (yeast.rows) == nrow (tbl.yeast))
  stdID [yeast.rows] = tbl.yeast$std
  
    # add human geneIDs
  tbl.human =  uniprotToStandardID ('human', human.uids)
  stopifnot (length (human.rows) == nrow (tbl.human))
  stdID [human.rows] = tbl.human$std

    # add worm geneIDs
  tbl.worm =  uniprotToStandardID ('worm', worm.uids)
  stopifnot (length (worm.rows) == nrow (tbl.worm))
  stdID [worm.rows] = tbl.worm$std

  invisible (stdID)

} # getAllStandardIDs
#------------------------------------------------------------------------------------------------------------------------
# change, e.g., "Homo sapiens" to "Hsapiens"
standardizeSpeciesNames = function (names)
{
   fix = function (name) {
     tokens = strsplit (name, ' ')[[1]]
     stopifnot (length (tokens) == 2)
     genus = tokens [1]
     species = tokens [2]
     return (paste (substr (genus, 1, 1), species, sep=''))
     } 

   fixed.names = as.character (sapply (names, fix))
   invisible (fixed.names)

} # standardizeSpeciesNames
#------------------------------------------------------------------------------------------------------------------------
renameMatrices = function (matrices, tbl.md)
{
  stopifnot (length (matrices) == nrow (tbl.md))
  names (matrices) = rownames (tbl.md)
  invisible (matrices)

} # renameMatrices
#------------------------------------------------------------------------------------------------------------------------
