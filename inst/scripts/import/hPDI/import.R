# hPDI/import.R
#------------------------------------------------------------------------------------------------------------------------
library (org.Hs.eg.db)
#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
run = function (dataDir)
{
  dataDir <- file.path(dataDir, "hPDI")

  filenames = getMatrixFilenames (dataDir)
  matrices = readMatrices (filenames)
  tbl.anno = readAnnotation (dataDir)

     # make sure that all of the matrix names are geneSymbols represented in tbl.anno

  stopifnot (length (intersect (names (matrices), tbl.anno$geneSymbol)) == length (names (matrices)))

  tbl.md = createMetadata (matrices, tbl.anno)

     # also required:  the order of tbl.anno entries exactly follow that of the matrices

  stopifnot (names (matrices) == tbl.md$providerName)
  matrices = renameMatrices (normalizeMatrices (matrices), tbl.md)
  serializedFile <- "hPDI.RData"
  save (matrices, tbl.md, file=serializedFile)
  printf("saved %d matrices to %s", length(matrices), serializedFile)
  printf("next step:  copy %s to <packageRoot>/MotifDb/inst/extdata, rebuild package", serializedFile)

  invisible (list (matrices=matrices, md=tbl.md))

} # run
#------------------------------------------------------------------------------------------------------------------------
getMatrixFilenames = function (dataDir)
{
  files = list.files (file.path(dataDir, "pwm"))
  invisible (file.path(dataDir, "pwm", files))

} # getMatrixFilenames
#------------------------------------------------------------------------------------------------------------------------
readMatrices = function (filenames)
{
  max = length (filenames)

  matrices = list ()

  for (i in 1:max) {
    filename = filenames [i]
    #gene.name = strsplit (filename, '\\.')[[1]][1]
    gene.name = sub(".output", "", basename(filename))
    mtx = as.matrix (read.table (filename, header=FALSE, sep='\t', row.names=1))
    colnames (mtx) = 1:ncol (mtx)
    matrices [[gene.name]] = mtx
    } # readMatrices

  invisible (matrices)

} # readMatrices
#------------------------------------------------------------------------------------------------------------------------
extractProteinIDs = function (raw.list)
{
     # first go after all the explicit uniprot ids, map them to refseq NP
  uniprot.rows = grep ('Source:Uniprot', raw.list)
  uniprot.ids.raw = raw.list [uniprot.rows]
  uniprot.names = sapply (strsplit ((uniprot.ids.raw), ';Acc:'), function (tokens) tokens [2])
  uniprot.names = gsub ("\\].*$", '', uniprot.names)
  mapped.refseq.names = uniprotToRefSeq (uniprot.names)

  protein = rep (NA_character_, length (raw.list))
  protein [uniprot.rows] = as.vector (mapped.refseq.names)

    # now grab the rows with explicit RefSeq NP idneitifers
  refseq.rows  = grep ('Source:RefSeq', raw.list)
  refseq.ids.raw = raw.list [refseq.rows]
  refseq.names = sapply (strsplit ((refseq.ids.raw), ';Acc:'), function (tokens) tokens [2])
  refseq.names = gsub ("\\].*$", '', refseq.names)

  protein [refseq.rows] = refseq.names
  still.na.rows = which (is.na (protein))
  #mystery.rows = setdiff (1:nrow (tbl.anno), c (uniprot.rows, refseq.rows))
  invisible (protein)

} # extractProteinIDs
#------------------------------------------------------------------------------------------------------------------------
readRawAnnotation = function (dataDir)
{
  filename <- file.path(dataDir, "protein_annotation.txt")
  tbl.anno = read.table (filename, sep='\t', header=T, fill=T, stringsAsFactors=FALSE, quote="")

    # oddly, there are 37 empty columns in this data.frame.  locate the first of these, then delete all
  stopifnot ('X' %in% colnames (tbl.anno))
  first.empty.column = which (colnames (tbl.anno) == 'X')
  tbl.anno = tbl.anno [, 1:(first.empty.column-1)]
  stopifnot (length (colnames (tbl.anno)) == 16)

  #protein.ids = extractIdentifiers (tbl.anno$ensembl_description)
     # the goal is to create a new 'protein' column for the table, all refseq NP identifers, or NA
  #tbl.anno = cbind (tbl.anno, protein, stringsAsFactors=FALSE)
  invisible (tbl.anno)

} # readRawAnnotation
#------------------------------------------------------------------------------------------------------------------------
addProteinColumn = function (tbl.anno)
{
  target = "ensembl_description"
  stopifnot (target %in% colnames (tbl.anno))
  ensembl.desc = tbl.anno [, target]
  protein = extractProteinIDs (ensembl.desc)
  invisible (cbind (tbl.anno, protein, stringsAsFactors=FALSE))

} # addProteinColumn
#------------------------------------------------------------------------------------------------------------------------
# remove all but symbol, locusID, Pfam_domains and protein.  rename the first three
trimAnnoTable = function (tbl.anno)
{
  tbl.fixed = tbl.anno [, c ('symbol', 'locusID', 'Pfam_domains', 'protein')]
  colnames (tbl.fixed) = c ('geneSymbol', 'geneID', 'pfamDomain', 'protein')

  invisible (tbl.fixed)

} # trimAnnoTable
#------------------------------------------------------------------------------------------------------------------------
# read raw, then process and trim to get a tidy 4-column table
readAnnotation = function (dataDir)
{
  tbl.anno.raw = readRawAnnotation (dataDir)
  tbl.anno = addProteinColumn (tbl.anno.raw)
  tbl.anno = trimAnnoTable (tbl.anno)

  invisible (tbl.anno)

} # readAnnotation
#------------------------------------------------------------------------------------------------------------------------
createMetadata = function (matrices, tbl.anno)
{
  options ('stringsAsFactors'=FALSE)
  dataSource = 'hPDI'
  organism = 'Hsapiens'
  provider.names = names (matrices)
  full.names = paste (organism, dataSource, provider.names, sep='-')

  tbl.md = data.frame (providerName=provider.names)
  rownames (tbl.md) = full.names
  tbl.md = cbind (tbl.md, providerId=provider.names)

    # indices are crucial, retreiving values from tbl.anno columns in an order that corresponds to the
    # order of the provider.names obtained from the names of the matrices

  indices = match (provider.names, tbl.anno$geneSymbol)
  geneId = as.character (tbl.anno$geneID [indices])

  egSyms = as.character (mget (geneId, org.Hs.egSYMBOL, ifnotfound=NA))
  egSyms [egSyms=='NA'] == NA_character_

  tbl.md = cbind (tbl.md, dataSource = rep (dataSource, nrow (tbl.md)))

  tbl.md = cbind (tbl.md, geneSymbol=egSyms)
  tbl.md = cbind (tbl.md, geneId)
  tbl.md = cbind (tbl.md, geneIdType=rep('ENTREZ', nrow (tbl.md)))

   
    # find NP protein ids for each matrix.  these are provided by hPDI, and made into NPs consistently by me.
 

  protein.ids = tbl.anno$protein [indices]
  tbl.md = cbind (tbl.md, proteinId=tbl.anno$protein [indices])

  proteinIdType = rep (NA_character_, length (protein.ids))
  refseq.protein.ids = grep ('NP_', protein.ids)
  proteinIdType [refseq.protein.ids] = 'REFSEQ'
  tbl.md = cbind (tbl.md, proteinIdType)

  tbl.md = cbind (tbl.md, organism = rep (organism, nrow (tbl.md)))

  counts = as.integer (sapply (matrices, function (mtx) as.integer (mean (round (colSums (mtx))))))
  tbl.md = cbind (tbl.md, sequenceCount=counts)

  tbl.md = cbind (tbl.md, bindingSequence=rep(NA_character_, nrow (tbl.md)))

  bindingDomain = tbl.anno$pfamDomain [indices]
  tbl.md = cbind (tbl.md, bindingDomain)

  tbl.md = cbind (tbl.md, tfFamily=rep(NA_character_, nrow (tbl.md)))

  experiment.types = rep (NA_character_, nrow (tbl.md))
  tbl.md = cbind (tbl.md, experimentType=experiment.types)

  pubmedIDs = rep (NA_character_)

  tbl.md = cbind (tbl.md, pubmedID=pubmedIDs)

  invisible (tbl.md)

} # createMetadata
#------------------------------------------------------------------------------------------------------------------------
uniprotToRefSeq = function (ids)
{
  if (!exists ('tbl.prot')) {
    tbl.refseq = toTable (org.Hs.egREFSEQ)
    np.only = grep ('NP_', tbl.refseq$accession)
    stopifnot (length (np.only) > 1000)   # crude sanity check
    tbl.refseq = tbl.refseq [np.only,]
    tbl.uniprot = toTable (org.Hs.egUNIPROT)
    tbl.prot <<- merge (tbl.refseq, tbl.uniprot, all.x=TRUE)
    } # create tbl.protx

   indices = match (ids, tbl.prot$uniprot_id)
   result = tbl.prot$accession [indices]
   names (result) = ids
   return (result)

} # uniprotToRefSeq
#------------------------------------------------------------------------------------------------------------------------
normalizeMatrices = function (matrices)
{
  mtx.normalized = sapply (matrices, function (mtx) apply (mtx, 2, function (colvector) colvector / sum (colvector)))
  invisible (mtx.normalized)

} # normalizeMatrices
#------------------------------------------------------------------------------------------------------------------------
renameMatrices = function (matrices, tbl.md)
{
  stopifnot (names (matrices) == tbl.md$providerName)  # same length, same order
  names (matrices) = rownames (tbl.md)
  invisible (matrices)

} # renameMatrices
#------------------------------------------------------------------------------------------------------------------------
