# ScerTF/import.R
#------------------------------------------------------------------------------------------------------------------------
library(RUnit)
library(org.Sc.sgd.db)
library(tools)   # for md5sum
#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
run = function (dataDir)
{
  dataDir <- file.path(dataDir, "ScerTF")
  freshStart ()
  tbl.ref = createExperimentRefTable ()
  all.files = getMatrixFilenames (dataDir)
  matrices = readAndParse (file.path(dataDir, all.files))
  tbl.md = createMetadata (matrices, tbl.ref)
  matrices = renameMatrices (matrices, tbl.md)
  serializedFile <- 'ScerTF.RData'
  save (matrices, tbl.md, file=serializedFile)
  printf("saved %d matrices to %s", length(matrices), serializedFile)
  printf("next step: copy %s to <packageRoot>/MotifDb/inst/extdata, rebuild package", serializedFile)

  invisible (list (mtx=matrices, md=tbl.md))

} # run
#------------------------------------------------------------------------------------------------------------------------
freshStart = function ()
{
  output.files.easy.to.regenerate = grep ('RData$', dir (), v=T)
  if (length (output.files.easy.to.regenerate) > 0)
     unlink (output.files.easy.to.regenerate)

} # freshStart
#------------------------------------------------------------------------------------------------------------------------
createMatrixNameUniqifier = function (matrix)
{
  temporary.file <<- tempfile ()
  write (as.character (matrix), file=temporary.file)
  md5sum.string <<- as.character (md5sum (temporary.file))
  stopifnot (nchar (md5sum.string) == 32)
  md5.6chars = substr (md5sum.string, 29, 32)
  #unlink (temporary.file)

} # createMatrixNameUniqifier
#------------------------------------------------------------------------------------------------------------------------
createExperimentRefTable = function ()
{
  options (stringsAsFactors = FALSE)
  tbl.ref = data.frame (author=c('badis', 'foat', 'fordyce', 'harbison', 'macisaac', 'morozov', 'pachkov', 'spivak', 'zhao', 'zhu'))
  tbl.ref = cbind (tbl.ref, year=c(2008,2008, 2010, 2004, 2006, 2007, 2007, 2012, 2011, 2009))
  tbl.ref = cbind (tbl.ref, pmid=c('19111667', '17947326', '20802496', '15343339', '16522208', '17438293', '17130146', 
                                    '22140105', '21654662', '19158363'))
  tbl.ref = cbind (tbl.ref, organism=rep('Scerevisiae', nrow (tbl.ref)))

  titles = c ('A library of yeast transcription factor motifs reveals a widespread function for Rsc3 in targeting nucleosome exclusion at promoters',
    'TransfactomeDB: a resource for exploring the nucleotide sequence specificity and condition-specific regulatory activity of trans-acting factors.', 
               'De novo identification and biophysical characterization of transcription-factor binding sites with microfluidic affinity analysis.',
               'Transcriptional regulatory code of a eukaryotic genome.',
               'An improved map of conserved regulatory sites for Saccharomyces cerevisiae',
               'Connecting protein structure with predictions of regulatory sites.',
               'SwissRegulon: a database of genome-wide annotations of regulatory sites.', 
               'ScerTF: a comprehensive database of benchmarked position weight matrices for Saccharomyces species.',
               'Quantitative analysis demonstrates most transcription factors require only simple models of specificity.',
                'High-resolution DNA-binding specificity analysis of yeast transcription factors')

  tbl.ref = cbind (tbl.ref, titles)
 
  tbl.ref

} # createExperimentRefTable
#------------------------------------------------------------------------------------------------------------------------
parsePWMfromText = function (lines.of.text)
{
  if (sum (sapply (c ('A', 'C', 'T', 'G'), function (token) length (grep (token, lines.of.text)))) != 4) {
    print (lines.of.text)
    return (NA)
    }     

  for (line in lines.of.text) {
    tokens = strsplit (line, '\\s*[:\\|]') [[1]]
    nucleotide = tokens [1]
    numbers.raw = tokens [2]
    number.tokens = strsplit (numbers.raw, '\\s+', perl=T)[[1]]
    while (nchar (number.tokens [1]) == 0)
      number.tokens = number.tokens [-1]
    numbers = as.numeric (number.tokens)
    if (!exists ('result'))
      result = matrix (nrow=4, ncol=length (numbers), byrow=TRUE, dimnames=list (c ('A', 'C', 'G', 'T'), 1:length(numbers)))
    result [nucleotide,] = numbers
    }

  return (result)

} # parsePWMfromText
#------------------------------------------------------------------------------------------------------------------------
# we expect exactly one matrix per file
readAndParse = function (filenames)
{
    # this script is in the directory with all of the PWM files.   exclude it

  files.to.remove = c ('go.R', 'matrices.ScerTF.RData', 'tbl.md.ScerTF.RData')
  for (file.to.remove in files.to.remove) { 
    removers = grep (file.to.remove, filenames)
    if (length (removers) > 0)
      filenames = filenames [-removers]
    } # for file.to.remove

  matrices = list ()

  for (filename in filenames) {
    lines.of.text = scan (filename, what=character(0), sep='\n', quiet=TRUE)
    mtx = parsePWMfromText (lines.of.text)
    mtx.normalized = apply (mtx, 2, function (colvector) colvector / sum (colvector))
    matrices [[basename(filename)]] = mtx.normalized
    }
 
   invisible (matrices)

} # readAndParse
#-----------------------------------------------------------------------------------------------------------------------
toOrf = function (gene.names, quiet=FALSE)
{
  orfs = mget (gene.names, org.Sc.sgdCOMMON2ORF, ifnotfound=NA)
  failures = which (is.na (orfs))
      # make sure that all of the failures are simply gene.names which are already orfs

    # some genes -- unpopular ones? :} -- have a classic orf name as gene symbol.  these do not apper in sgdCOMMON2ORF
    # don't protest these, rather, just notify the caller of *other* symbols which failed to map
    # this failure situation can be recognized if there are more failures (NAs) than orf names (Y.*) found in the input
    # gene.names

  if (length (grep ('^Y', names (failures))) != length (failures)) {
    if (!quiet) {
       printf ('error.  could not find orf for "%s"', list.to.string (gene.names [failures]))
       } # !quiet
    } 
 
  orfs [as.integer (failures)] = names (failures)
  for (gene.name in names (orfs)) {
     orfs [[gene.name]] = orfs [[gene.name]][1]
     } # for gene.name

  unlist (unname (orfs))

} # toOrf
#------------------------------------------------------------------------------------------------------------------------
toUniprot = function (orfs, quiet=FALSE)
{
  uniprots = mget (orfs, org.Sc.sgdUNIPROT, ifnotfound=NA)
  failures = which (is.na (uniprots))

  for (orf in names (uniprots)) {
     uniprots [[orf]] = uniprots [[orf]][1]
     } # for orf

  unlist (unname (uniprots))

} # toUniprot
#------------------------------------------------------------------------------------------------------------------------
# and rename the matrices to match the rownames of the tbl.md created here
createMetadata = function (matrices, tbl.ref)
{
  options ('stringsAsFactors'=FALSE)
  dataSource = 'ScerTF'
  organism = 'Scerevisiae'

  #tbl.ref = createExperimentRefTable ()
  native.names.raw = names (matrices)
  native.names.reversed = as.character (sapply (native.names.raw, 
                                 function (s) {tokens = strsplit (s, '\\.')[[1]]; return (paste (tokens [2], tokens[1], sep='-'))}))

  native.names = as.character (sapply (names (matrices), function (name) strsplit (name, '\\.') [[1]][2]))
  gene.symbols = native.names

  md5sum.suffix.uniqifiers = as.character (sapply (matrices, function (matrix) createMatrixNameUniqifier (matrix)))
  native.names.uniqified = paste (native.names, md5sum.suffix.uniqifiers, sep='-')
  full.names = paste (organism, dataSource, native.names.reversed, sep='-')

  tbl.md = data.frame (providerName=native.names.raw)
  rownames (tbl.md) = full.names
  names (matrices) = full.names
  tbl.md = cbind (tbl.md, providerId=gene.symbols)
  tbl.md = cbind (tbl.md, dataSource = rep ('ScerTF', nrow (tbl.md)))
  tbl.md = cbind (tbl.md, geneSymbol=gene.symbols)

  orfs = as.character (toOrf (native.names))
  tbl.md = cbind (tbl.md, geneId=orfs)
  tbl.md = cbind (tbl.md, geneIdType=rep('SGD', nrow (tbl.md)))

  uniprots <<- toUniprot (orfs)
  tbl.md = cbind (tbl.md, proteinId=uniprots)
  protein.id.types = rep('UNIPROT', nrow (tbl.md))
  NA.protein.ids = which (is.na (uniprots))
  if (length (NA.protein.ids) > 0)
    protein.id.types [NA.protein.ids] = NA
  
  tbl.md = cbind (tbl.md, proteinIdType=protein.id.types)
    
  tbl.md = cbind (tbl.md, organism = rep ('Scerevisiae', nrow (tbl.md)))

    # all matrices have colSums of 100, or very close, suggesting that these are not real counts but are normalized
  tbl.md = cbind (tbl.md, sequenceCount = rep (NA_integer_, nrow (tbl.md)))
  tbl.md = cbind (tbl.md, bindingSequence=rep(NA_character_, nrow (tbl.md)))
  tbl.md = cbind (tbl.md, bindingDomain=rep(NA_character_, nrow (tbl.md)))
  

  tbl.md = cbind (tbl.md, tfFamily=rep(NA_character_, nrow (tbl.md)))
  tbl.md = cbind (tbl.md, experimentType=rep(NA_character_, nrow (tbl.md)))
  authors = as.character (sapply (tbl.md$providerName, function (provider.name) strsplit (provider.name, '\\.')[[1]][1]))
  pmid = as.character (sapply (authors, function (athr) subset (tbl.ref, author==athr)$pmid))
  tbl.md = cbind (tbl.md, pubmedID=pmid)
  
  tbl.md

} # createMetadata
#-----------------------------------------------------------------------------------------------------------------------
getMatrixFilenames = function (dataDir)
{
  all.files = list.files (dataDir)
  files.to.exclude = c ('go.R', 'RCS', 'rdata')
  for (file in files.to.exclude) {
    removers = grep (file, all.files, ignore.case=T)
    if (length (removers) > 0)
      all.files = all.files [-removers]
    }

  invisible (all.files)

} # getMatrixFilenames
#-----------------------------------------------------------------------------------------------------------------------
renameMatrices = function (matrices, tbl.md)
{
  stopifnot (length (matrices) == nrow (tbl.md))
  names (matrices) = rownames (tbl.md)
  invisible (matrices)

} # renameMatrices
#------------------------------------------------------------------------------------------------------------------------
