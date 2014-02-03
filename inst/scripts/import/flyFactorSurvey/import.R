# flyFactorSurvey/import.R
#------------------------------------------------------------------------------------------------------------------------
library(org.Dm.eg.db)
library(biomaRt)
#------------------------------------------------------------------------------------------------------------------------
bindingDomainXrefSourceFile <- function() {"TFfile2b.tsv"}
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
run = function (dataDir)
{
  dataDir <- file.path(dataDir, "flyFactorSurvey")
  filenames = getMatrixFilenames (dataDir)

     # the rootDir has only matrix files, with one exception:
     # a file which describes the protein binding domains
     # remove that from the matrix file list
  
  removers <- grep (bindingDomainXrefSourceFile(), filenames)
  if(length(removers) > 0)
      filenames <- filenames[-removers]
  
  full.filenames = file.path(dataDir, filenames)

  list.BD = createXrefBindingDomain (dataDir)
  xref = createXref ()  # maps flybase ids to uniprot
  tbl.ref = createExperimentRefTable ()
  matrices = readAndParse (full.filenames)
  tbl.md = createMetadata (matrices, tbl.ref, xref, list.BD)
     # now normalize, not in readAndParse, because we need the counts to go into tbl.md$sequenceCount
  matrices = normalizeMatrices (matrices)
  matrices = renameMatrices (matrices, tbl.md)
  serializedFile <- "flyFactorSurvey.RData"
  save (matrices, tbl.md, file=serializedFile)
  printf("saved %d matrices to %s", length(matrices), serializedFile)
  printf("next step: copy %s to <packageRoot>/MotifDb/inst/extdata, and rebuild package", serializedFile)

} # run
#------------------------------------------------------------------------------------------------------------------------
# rec'd extra xref info from michael brodsky on (18 jul 2012)
# protein binding domains are of immediate interest.
createXrefBindingDomain = function (dataDir)
{
  f1 = file.path(dataDir, bindingDomainXrefSourceFile())
  tbl.raw = read.table (f1, sep='\t', header=TRUE, as.is=TRUE, comment.char='')
  list.BD = tbl.raw$DNAbindingDomain_Primary
  names (list.BD) = tbl.raw$FlybaseID
  invisible (list.BD)

} # createXrefBindingDomain
#------------------------------------------------------------------------------------------------------------------------
# a named list, with uniprot id as content, flybase gene names (FBgn) as names
createXref = function ()
{
  tbl.egUNIPROT = toTable (org.Dm.egUNIPROT)
  tbl.egFLYBASE = toTable (org.Dm.egFLYBASE)
  tbl.xref = merge (tbl.egFLYBASE, tbl.egUNIPROT)
  xref = tbl.xref [, 3]
  names (xref) = tbl.xref [, 2]
  invisible (xref)

} # createXrefp
#------------------------------------------------------------------------------------------------------------------------
getMatrixFilenames = function (dataDir)
{
  all.files = list.files (dataDir)
  files.to.exclude = c ('go.R', 'RCS', 'rdata')
  for (file in files.to.exclude) {
    removers = grep (file, all.files, ignore.case=T)
    if (length (removers) > 0)
      all.files = all.files [-removers]
    } # for

  invisible (sort(all.files))

} # getMatrixFilenames
#-----------------------------------------------------------------------------------------------------------------------
createExperimentRefTable = function ()
{
  experiment.names = c ('SANGER', 'SOLEXA', 'NAR', 'Cell', 'FlyReg', 'NBT')
  pubmedID = c (NA_character_, NA_character_, '1858360', '18332042', NA_character_, NA_character_)
  experimentType = c ('bacterial 1-hybrid, SANGER sequencing',
                      'bacterial 1-hybrid, SOLEXA sequencing',
                      'bacterial 1-hybrid, SANGER sequencing',
                      'bacterial 1-hybrid, SANGER sequencing',
                      'bacterial 1-hybrid, SANGER sequencing',
                      'bacterial 1-hybrid, SANGER sequencing')

  tbl.ref = data.frame (pmid=pubmedID, experimentType=experimentType, stringsAsFactors=FALSE)
  rownames (tbl.ref) = experiment.names
  invisible (tbl.ref)

} # createExperimentRefTable
#------------------------------------------------------------------------------------------------------------------------
parsePWMfromText = function (lines.of.text)
{
  stopifnot (sum (sapply (c ('A', 'C', 'T', 'G'), function (token) length (grep (token, lines.of.text)))) == 4)

  for (line in lines.of.text) {
    tokens = strsplit (line, '\\s*[:\\|]') [[1]]
    nucleotide = tokens [1]
    numbers.raw = tokens [2]
    number.tokens = strsplit (numbers.raw, '\\s+', perl=T)[[1]]
    while (nchar (number.tokens [1]) == 0)
      number.tokens = number.tokens [-1]
    numbers = as.numeric (number.tokens)
    if (!exists ('pwm.result'))
      pwm.result = matrix (nrow=4, ncol=length (numbers), byrow=TRUE, dimnames=list (c ('A', 'C', 'G', 'T'), 1:length(numbers)))
    pwm.result [nucleotide,] = numbers
    }

  pwm.result

} # parsePWMfromText
#------------------------------------------------------------------------------------------------------------------------
# we expect exactly one matrix per file
readAndParse = function (filenames)
{
  matrices = list ()

  for (filename in filenames) {
      lines.of.text = scan (filename, what=character(0), sep='\n', comment='#', quiet=TRUE)
      mtx = parsePWMfromText (lines.of.text)
      matrices [[basename(filename)]] = mtx
       }
 
   invisible (matrices)

} # readAndParse
#-----------------------------------------------------------------------------------------------------------------------
createMetadata = function (matrices, tbl.ref, xref, list.BD)
{
  options ('stringsAsFactors'=FALSE)
  dataSource = 'FlyFactorSurvey'
  organism = 'Dmelanogaster'
  stopifnot (length (grep ('\\.pfm$', names (matrices))) == length (matrices))
  native.names.raw = gsub ('\\.pfm$', '', names (matrices))
  #geneSymbols = sapply (strsplit (native.names.raw, '_'), function (tokens) return (tokens [1]))
  flybase.ids = sapply (strsplit (native.names.raw, '_'), function (tokens) return (tokens [length (tokens)]))
  tbl.ids <- fbgnToIDs(flybase.ids)
     # we may have duplicate FBgns in the otherwise unique matrix names,
     # creating a tbl.ids with fewer rows than we have matrices.
     # expand tbl.ids to align with the matrix nmes
  tbl.ids <- tbl.ids[flybase.ids,]
  stopifnot(nrow(tbl.ids) == length(matrices))
     #  next up: call fbgnToIDs(flybase.ids), assign geneSymbols
     #  and other variables from the returned data.frame
  stopifnot (length (grep ('^FBgn', flybase.ids)) ==  length (flybase.ids))
  stopifnot (all (nchar (flybase.ids) == 11))
  native.names.cooked = gsub ('-', '.', native.names.raw)   # so that the dash reliably separates the three parts of full.names
  full.names = paste (organism, dataSource, native.names.cooked, sep='-')

  tbl.md = data.frame (providerName=native.names.raw)
  rownames (tbl.md) = full.names
  tbl.md = cbind (tbl.md, providerId=flybase.ids)
  tbl.md = cbind (tbl.md, dataSource = rep (dataSource, nrow (tbl.md)))

  tbl.md = cbind (tbl.md, geneSymbol=tbl.ids$symbol)
  tbl.md = cbind (tbl.md, geneId=tbl.ids$gene_id)
  #geneIdType <- rep("ENTREZ", nrow(tbl.md))
     # not all FBgn ids can be mapped (by fbgnToIds) to entrez genes; for those,
     # the gene_id field is the FBgn id
  #geneIdIsFBgn <- grep("^FBgn", tbl.ids$gene_id)
  #if(length(geneIdIsFBgn) > 0)
  #    geneIdType [geneIdIsFBgn] <- "FLYBASE"
  tbl.md = cbind (tbl.md, geneIdType=tbl.ids$geneIdType)

  uniprot.ids = as.character (xref [flybase.ids])
  tbl.md = cbind (tbl.md, proteinId=uniprot.ids)
  protein.id.types = rep ('UNIPROT', nrow (tbl.md))
  NA.protein.indices = which (is.na (tbl.md$proteinId))
  if (length (NA.protein.indices) > 0)
    protein.id.types [NA.protein.indices] = NA
  tbl.md = cbind (tbl.md, proteinIdType=protein.id.types)
  #tbl.md = cbind (tbl.md, proteinIdType=rep ('UNIPROT', nrow (tbl.md)))

  tbl.md = cbind (tbl.md, organism = rep (organism, nrow (tbl.md)))

  counts = as.integer (sapply (matrices, function (mtx) as.integer (mean (round (colSums (mtx))))))
  tbl.md = cbind (tbl.md, sequenceCount=counts)

  tbl.md = cbind (tbl.md, bindingSequence=rep (NA_character_, nrow (tbl.md)))
  bindingDomains = rep (NA, nrow (tbl.md))
  indices = match (tbl.md$providerId, names (list.BD))
  bindingDomain = as.character (list.BD)[indices]
  bindingDomain [bindingDomain == 'NA'] = NA
  bindingDomain [bindingDomain == ''] = NA
  tbl.md = cbind (tbl.md, bindingDomain)

  experiment.types = rep (NA_character_, nrow (tbl.md))
  from.Cell = grep ('_Cell_', native.names.raw)
  from.NAR  = grep ('_NAR_', native.names.raw)
  from.SOLEXA = grep ('_SOLEXA_', native.names.raw)
  from.SANGER = grep ('_SANGER_', native.names.raw)
  from.FlyReg = grep ('_FlyReg_', native.names.raw) 

  experiment.types [from.SOLEXA] = 'bacterial 1-hybrid, SOLEXA sequencing'
  experiment.types [from.SANGER] = 'bacterial 1-hybrid, SANGER sequencing'
  experiment.types [from.Cell] = 'bacterial 1-hybrid, SANGER sequencing'
  experiment.types [from.NAR] = 'bacterial 1-hybrid, SANGER sequencing'
  experiment.types [from.FlyReg] = 'DNASE I footprinting'

  tbl.md = cbind (tbl.md, tfFamily=rep(NA_character_, nrow (tbl.md)))

  tbl.md = cbind (tbl.md, experimentType=experiment.types)

  pubmedIDs = rep (NA_character_, nrow (tbl.md))
  pubmedIDs [from.Cell] = '18332042'
  pubmedIDs [from.NAR] = '1858360'

  tbl.md = cbind (tbl.md, pubmedID=pubmedIDs)

  invisible (tbl.md)

} # createMetadata
#------------------------------------------------------------------------------------------------------------------------
normalizeMatrices = function (matrices)
{
  mtx.normalized = sapply (matrices, function (mtx) apply (mtx, 2, function (colvector) colvector / sum (colvector)))
  invisible (mtx.normalized)

} # normalizeMatrices
#------------------------------------------------------------------------------------------------------------------------
renameMatrices = function (matrices, tbl.md)
{
  stopifnot (length (matrices) == nrow (tbl.md))
  names (matrices) = rownames (tbl.md)
  invisible (matrices)

} # renameMatrices
#-------------------------------------------------------------------------------
fbgnToIDs <- function(fbgns, useInputForMissingValues=TRUE)
{
   fbgns <- unique(fbgns)
   if(!exists("tbl.names")) {
       tbl.sym <- toTable(org.Dm.egSYMBOL)
       tbl.fbgn <- toTable(org.Dm.egFLYBASE)
       tbl.names <<- merge(tbl.fbgn, tbl.sym)[, c("flybase_id", "symbol", "gene_id")]
       flyMart <<- useMart(biomart="ensembl", dataset="dmelanogaster_gene_ensembl")
       }
   tbl.bm <- getBM(attributes=c("flybase_gene_id", "entrezgene"), filters=c("flybase_gene_id"),
                   values=fbgns, mart=flyMart)
   fbgns.notfound <- setdiff(fbgns, tbl.bm$flybase_gene_id)
   if(length(fbgns.notfound) > 0){
       new.rows <- data.frame(flybase_gene_id=fbgns.notfound,
                              entrezgene=NA,
                              stringsAsFactors=FALSE)
       tbl.bm <- rbind(tbl.bm, new.rows)
       }

   tbl.bm <- unique(tbl.bm)

      # some fbgns could be mapped to multiple entrez genes.  remove all but the
      # one which occurs first
   
   duplicates <- which(duplicated(tbl.bm$flybase_gene_id))
   if(length(duplicates) > 0)
       tbl.bm <- tbl.bm[-duplicates,]
   
   stopifnot(sort(fbgns) == sort(tbl.bm$flybase_gene_id))

       # order the biomart results so that it matches the order of the incoming fgbns
   tbl.bm <- tbl.bm[match(fbgns, tbl.bm$flybase_gene_id),]
   hits <- match(tbl.bm$entrezgene, tbl.names$gene_id)
   result <- tbl.names[hits,]
   rownames(result) <- fbgns

       # not all FBgn ids have an associated entrez geneID.
       # use these orphan FBgn ids as their own "geneID"
   
   failed.lookups <- which(is.na(result$flybase_id))
   failed.count <- length(failed.lookups)
   if(useInputForMissingValues & failed.count > 0) {
      failed.names <- rownames(result)[failed.lookups]
        # an n x 3 matrix, all entries in each row the same
      names.matrix <- matrix(rep(failed.names, 3), nrow=failed.count, byrow=FALSE)
      result[failed.lookups,] <- names.matrix
      }

   geneIdType <- rep("ENTREZ", nrow(result))
   geneIdType [failed.lookups] <- "FLYBASE"
   result <- cbind(result, geneIdType, stringsAsFactors=FALSE)
   result


} # fbgnToIDs
#------------------------------------------------------------------------------------------------------------------------
