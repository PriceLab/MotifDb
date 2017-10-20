# MotifDb/inst/scripes/import/cisbp/import.R
#------------------------------------------------------------------------------------------------------------------------
options (stringsAsFactors=FALSE)
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
run = function (dataDir)
{
  printf("cisbp matrix and metadata import");
  createMotifDbArchiveFile(dataDir, "cisbp.RData", count=NA)

} # run
#------------------------------------------------------------------------------------------------------------------------
# cisbp distributes metadata in a MySQL database with many tables.  i have cobbled together an sql query to
# extract a metadata data.frame from this database.  first, however, here are my notes on creating and filling
# that database.
# not that an alternative to all that fuss is to simply used this file, created by me and with all the
# information we currenly need:
#
#    file.path(dataDir, "cisbp", "cisbp-metadata-6559-matrices.Rdata")
#
# to create from scratch:
#   cd dataDir/cisbp/sqlTables
#   --- create accounts
#      bash> mysql -u root
#      mysql> create user 'pshannon'@'localhost';
#      mysql> grant all privileges on *.* to 'pshannon'@'localhost';
#      mysql> quit
#      bash> mysql -u pshannon
#      mysql> create database cisbp;
#      mysql> quit
#   --- load the data
#      bash> bash load.sh
#
createMetadataTableFromDatabase <- function(dataDir, user)
{

      # we want metadata only for those matrices found in the cisbp pwms download directory
      # so begin by getting those nams

   pwm.directory <- file.path(dataDir, "cisbp", "pwms")    
   matrix.names <- list.files(pwm.directory)
   matrix.names <- sub(".txt", "", matrix.names, fixed=TRUE)
   printf("matrix names: %d", length(matrix.names))
   s <- paste(matrix.names, collapse="','")
   formatted.matrix.names.as.group <- sprintf("'%s'", s)

   require(RMySQL)
   db <- dbConnect(MySQL(), user='pshannon', dbname='cisbp')
   matrix.selector <- sprintf("where ma.Motif_ID in (%s)", formatted.matrix.names.as.group)

   select <- "select ma.ma_id, ma.tf_id, ma.motif_id, ma.species, tf.TF_Name, ms.PMID, fa.Family_Name, pr.DBID"
   from <- "from motif_all as ma, tfs as tf, motifs as mo, motif_sources as ms, tf_families as fa, proteins as pr"
   where <- paste(matrix.selector,
                  "and ma.Motif_ID=mo.Motif_ID",
                  "and ma.TF_ID = tf.TF_ID",
                  "and ma.Evidence = 'D'",
                  "and mo.MSource_ID = ms.MSource_ID",
                  "and tf.Family_ID = fa.Family_ID",
                  "and pr.TF_ID = tf.TF_ID"
                  )

   query <- paste(select, from, where, sep=" ")
   tbl <- dbGetQuery(db, query)
   print(dim(tbl))
   print(dim(unique(tbl[, 1:8])))   # [1] 15694     8
   print(dim(unique(tbl[, 1:7])))   # [1]  6559    7
   dups <- which(duplicated(tbl[, 1:7]))
   if(length(dups) > 0)
      tbl <- tbl[-dups,]
   dim(tbl)
   filename <- file.path(dataDir, "cisbp", "cisbp-metadata-6559-matrices.Rdata")
   printf("saving metdata data.frame to %s", filename)
   save(tbl, file=filename)

} # createMetadataTableFromDatabase
#------------------------------------------------------------------------------------------------------------------------
readRawMatrices = function (dataDir)
{
    # our convention is that there is a shared "dataDir" visible to
    # the importer, and that within that directory there is one
    # subdirectory for each data source.
    # for this example importer, that directory will be <dataDir>/test
    # within which we will look for one small file "sample.pcm"
    
  filename <- file.path(dataDir, "cisbp", "sample.pcm")
  printf("checking for readable matrix file:")
  printf("     %s", filename)
  stopifnot(file.exists(filename))
  
  all.lines = scan (filename, what=character(0), sep='\n', quiet=TRUE)
  title.lines = grep ('^>', all.lines)
  title.line.count <<- length (title.lines)
  max = title.line.count - 1

  pwms = list ()
  
  for (i in 1:max) {
    start.line = title.lines [i]
    end.line = title.lines [i+1] - 1
    new.pwm = parsePwm (all.lines [start.line:end.line])
    pwms = c (pwms, list (new.pwm))
    } # for i

  
  invisible (pwms)

} # readRawMatrices
#------------------------------------------------------------------------------------------------------------------------
extractMatrices = function (pwm.list)
{
  matrices = sapply (pwm.list, function (element) element$matrix)
  matrix.names <- sapply (pwm.list, function (element) element$title)
  matrix.names <- sub("^>", "", matrix.names)
  names (matrices) <- matrix.names
  
  matrices

} # extractMatrices
#------------------------------------------------------------------------------------------------------------------------
# where incoming x is a list with elements like this
#
#   $ma_id =[1] "MA0116870_1.02"
#   $tf_id = [1] "T025973_1.02"
#   $motif_id = [1] "M0316_1.02"
#   $species = [1] "Mus_musculus"
#   $TF_Name = [1] "Nfil3"
#   $PMID = [1] "25215497"
#   $Family_Name = [1] "bZIP"
#   $DBID = [1] "ENSMUSP00000065363"
#
translateMetadataToMotifDbStandardForm <- function(x)
{
  checkTrue(all(c("tf_id", "motif_id", "species", "TF_Name", "PMID", "Family_Name", "DBID") %in% names(x)))
  standardizedOrganism <- standardizeOrganism(x$species)
  proteinIdType <- deduceProteinIdType(x$DBID, standardizedOrganism)
  
  std <- list(providerName=x$motif_id,
              providerId=x$ma_id,
              dataSource="cisbp_1.02",
              geneSymbol=x$TF_Name,
              geneId=NA,
              geneIdType=NA,
              proteinId=x$DBID,
              proteinIdType=proteinIdType,
              organism=standardizedOrganism,
              sequenceCount=NA,
              bindingSequence=NA,
              bindingDomain=NA,
              tfFamily=x$Family_Name,
              experimentType=NA,
              pubmedID=x$PMID)

  std  

} # translateMetadataToMotifDbStandardForm
#------------------------------------------------------------------------------------------------------------------------
deduceProteinIdType <- function(id, organism)
{
  if(substr(id, 1, 1) == "Y" & organism == "Scerevisiae")
     return("SGD")

  if(substr(id, 1, 2) == "EN")
     return("ENSEMBL")

  return(NA)
  
} # deduceProteinIdType
#------------------------------------------------------------------------------------------------------------------------
standardizeOrganism <- function(x)
{
   #printf("standardizeOrganism: %s", x);
   tokens <- strsplit(x, "_")[[1]]

   if(length(tokens) != 2){
      warning(sprintf("cisbp import could not standardize species name: '%s'", x))
      return(x)
      }

   result <- sprintf("%s%s", substr(tokens[1], 1, 1), tokens[2])

   result

} # standardizeOrganism
#------------------------------------------------------------------------------------------------------------------------
createMotifDbArchiveFile <- function(dataDir, RDataFileName, count=NA)
{
   if(is.na(count))
      file.names <- list.files(file.path(dataDir, "cisbp", "pwms"))
   else
      file.names <- head(list.files(file.path(dataDir, "cisbp", "pwms")), n=count)

   parseFunction <- function(filename){
      full.path <- file.path(dataDir, "cisbp", "pwms", filename)
      checkTrue(file.exists(full.path))
      text <- scan(full.path, sep="\n", quiet=TRUE, what=character(0)) 
      pwm <- parsePwm(title, text)
      pwm$matrix
      }

   matrices <- lapply(file.names, parseFunction)
   titles <- unlist(lapply(file.names, function(filename) sub(".txt", "", filename)))
   names(matrices) <- titles

   if(!exists("tbl")){
      metadata.file <- file.path(dataDir, "cisbp", "cisbp-metadata-6559-matrices.Rdata")
      load(metadata.file, env=.GlobalEnv)
      }

   count <- length(titles)
   tbl.md <- as.data.frame(list(providerName=vector("character", count),
                                providerId=vector("character", count),
                                dataSource=vector("character", count),
                                geneSymbol=vector("character", count),
                                geneId=vector("character", count),
                                geneIdType=vector("character", count),
                                proteinId=vector("character", count),
                                proteinIdType=vector("character", count),
                                organism=vector("character", count),
                                sequenceCount=vector("character", count),
                                bindingSequence=vector("character", count),
                                bindingDomain=vector("character", count),
                                tfFamily=vector("character", count),
                                experimentType=vector("character", count),
                                pubmedID=vector("character", count)))

   tbl.md.row = 0
   for(title in titles){
     tbl.md.row <- tbl.md.row + 1
     if(!title %in% tbl$motif_id){
        #printf("skipping matrix %d, no metadata for %s", tbl.md.row, title)
        next
        }
     row <- grep(title, tbl$motif_id)
     md <- as.list(tbl[row,])
     md.fixed <- translateMetadataToMotifDbStandardForm(md)
     tbl.md[tbl.md.row,] <- as.data.frame(md.fixed)
     } # for title

   empties <- which(nchar(tbl.md$providerName) == 0)
   if(length(empties) > 0){
      tbl.md <- tbl.md[-empties,]
      }
   rownames(tbl.md) <- paste(tbl.md$organism, tbl.md$dataSource, tbl.md$providerName, sep="-")
   matrices <- matrices[-empties]
   names(matrices) <- rownames(tbl.md)

   printf("saving %d matrices with metadata to %s", nrow(tbl.md), file.path(getwd(), RDataFileName))
   save(tbl.md, matrices, file=RDataFileName)

     # rebuild & install MotifDb, then watch as it loads:
     #  x <- MotifDb:::.MotifDb(TRUE, FALSE)
     #  mcols(query(x, "csipb"))
   
} # createMotifDbArchiveFile
#------------------------------------------------------------------------------------------------------------------------
## # files are named by Motif_ID, which also provides the database key used to create the metadata entries
## # for each motif's matrix, "M1093_1.02.txt" and "M1093_1.02"
## # cisbp at version 1.02 has 6559 matrices.  each of these is annotated to different organisms
## # producing maybe 70k metadata table entries
## # M1093
## createMetadataTable = function (dataDir, motifIDs)
## {
## 
##   select <- "select mb.Motif_ID, tf.TF_Name, Species, PMID"
##   from <- "from  motif_best as mb, tfs as tf, motifs as mo, motif_sources as ms"
##   where <- "where  mb.Motif_ID='M1903_1.02' and mb.Motif_ID=mo.Motif_ID and mb.TF_ID = tf.TF_ID and mo.MSource_ID=ms.MSource_ID"
## 
##   select <- "select mb.Motif_ID, tf.TF_Name, Species, PMID, Family_Name, pr.DBID"
##   from <- "from  motif_best as mb, tfs as tf, motifs as mo, motif_sources as ms, tf_families as fa, proteins as pr"
##   where <- paste(sprintf("where mb.Motif_ID='%s' and ", "M1903_1.02"),
##                  "mb.Motif_ID=mo.Motif_ID and",
##                  "mb.TF_ID = tf.TF_ID and",
##                  "mo.MSource_ID = ms.MSource_ID and",
##                  "tf.Family_ID = fa.Family_ID and",
##                  "pr.TF_ID = tf.TF_ID")
## 
##   tbl <- dbGetQuery(db, paste(select, from, where, sep=" "))
##      # multiple protein identifiers are often returned; we need just one
##   dups <- which(duplicated(tbl[, -(grep("DBID", colnames(tbl)))]))
##   if(length(dups) > 0)
##       tbl <- tbl[-dups,]
##   
##   filename <- file.path(dataDir, "cisbp", "md-raw.tsv")
##   printf("checking for readable metadata file:")
##   printf("   %s", filename)
##   stopifnot(file.exists(filename))
## 
##   tbl.raw <- read.table(filename, sep="\t", header=TRUE, as.is=TRUE)
##   tbl.md = data.frame ()
##   matrix.ids = names (matrices)
##   
##   for (matrix.id in matrix.ids) {
##     matrix <- matrices[[matrix.id]]
##     short.matrix.name <- sub("\\..*$", "", matrix.id)
##     stopifnot(length(grep(short.matrix.name, tbl.raw$ma.name)) == 1)
##     md <- as.list(subset(tbl.raw, ma.name==short.matrix.name))
##     dataSource <- "cisbp"
##     #browser()
##     organism <- ncbiTaxonimicCodeToLinnaean(md$ncbi.tax.code)
##     new.row = list (providerName=matrix.id,
##                     providerId=matrix.id,
##                     dataSource=dataSource,
##                     geneSymbol=md$gene.symbol,
##                     geneId=NA,
##                     geneIdType=NA,
##                     proteinId=md$uniprot,
##                     proteinIdType=guessProteinIdentifierType(md$uniprot),
##                     organism=ncbiTaxonimicCodeToLinnaean(md$ncbi.tax.code),
##                     sequenceCount=max(colSums(matrix)),
##                     bindingSequence=NA_character_,
##                     bindingDomain=NA,
##                     tfFamily=md$family,
##                     experimentType=md$type,
##                     pubmedID="24194598")
##     tbl.md = rbind (tbl.md, data.frame (new.row, stringsAsFactors=FALSE))
##     full.name = sprintf ('%s-%s-%s', organism, dataSource, matrix.id)
##     rownames (tbl.md) [nrow (tbl.md)] = full.name
##     } # for i
## 
##   invisible (tbl.md)
## 
## } # createMetadataTable
## #------------------------------------------------------------------------------------------------------------------------
## renameMatrices = function (matrices, tbl.md)
## {
##   stopifnot (length (matrices) == nrow (tbl.md))
##   names (matrices) = rownames (tbl.md)
##   invisible (matrices)
## 
## } # renameMatrices
## #------------------------------------------------------------------------------------------------------------------------
## # an empirical and not altogether trustworthy solution to identifying identifier types.
## guessProteinIdentifierType = function (moleculeName)
## {
##   if (nchar (moleculeName) == 0)
##     return (NA_character_)
##   if (is.na (moleculeName))
##     return (NA_character_) 
## 
##   first.char = substr (moleculeName, 1, 1)
## 
##   if (first.char == 'Y')
##     return ('SGD')
## 
##   if (first.char %in% c ('P', 'Q', 'O', 'A', 'C'))
##     return ('UNIPROT')
## 
##   if (length (grep ('^NP_', moleculeName)) == 1)
##     return ('REFSEQ')
## 
##    return (NA_character_)
## 
## } # guessProteinIdentifierType
## #------------------------------------------------------------------------------------------------------------------------
## normalizeMatrices = function (matrices)
## {
##   mtx.normalized = sapply (matrices,
##       function (mtx) apply (mtx, 2, function (colvector) colvector / sum (colvector)))
## 
##   invisible (mtx.normalized)
## 
## } # normalizeMatrices
## #------------------------------------------------------------------------------------------------------------------------
parsePwm = function (title, text)
{
   lines <- strsplit(text, '\t')
   stopifnot(lines[[1]] == c("Pos", "A", "C", "G", "T"))
   line.count <- length(lines)
   row.count <- length(lines) - 1
   col.count <- 4

     # our standard form is 4 rows (one per nucelotide) and n columns
     # cisbp matrices come in transposed: read them as-is, then transpose
     # to our format

  result <- matrix(nrow=row.count, ncol=col.count,
                   dimnames=list(as.character(1:row.count),
                                 c("A","C", "G", "T")))
                                 
  for(i in 1:row.count){
    row <- i + 1
    result [i,] = as.numeric (lines[[row]][2:5])
    } # for i

  result = t (result)
  return (list (title=title, matrix=result))

} # parsePwm
## #----------------------------------------------------------------------------------------------------
## ncbiTaxonimicCodeToLinnaean <- function(code)
## {
##   code <- as.character(code)
##   
##   lookup <- list("3702" = "Arabidopsis thaliana",
##                  "3888" = "Pisum sativum",
##                  "4094" = "Nicontia sp.",
##                  "4102" = "Petunia hybrida",
##                  "4151" = "Antirrhinum majus",
##                  "4513" = "Hordeum vulgare",
##                  "4565" = "Triticum aestivam",
##                  "4577" = "Zea mays",
##                  "4932" = "Saccaromyces cerevisiae",
##                  "6239" = "Caenorhabditis elegans",
##                  "7227" = "Drosophila melanogaster",
##                  "7729" = "Halocynthia roretzi",
##                  "7742" = "Vertebrata",
##                  "8022" = "Onchorhynchus mykiss",
##                  "8355" = "Xenopus laevis",
##                  "8364" = "Silurana tropicalis",
##                  "9031" = "Gallus gallus",
##                  "9606" = "Homo sapiens",
##                  "9913" = "Bos taurus",
##                  "9986" = "Oryctolagus cuniculus",
##                  "10090" = "Mus musculus",
##                  "10116" = "Rattus norvegicus",
##                  "10117" = "Rattus rattus")
## 
##   if (code %in% names(lookup))
##       return(lookup[[code]])
## 
##   NA
## 
## } # ncbiTaxonimicCodeToLinnaean
## #----------------------------------------------------------------------------------------------------
## 
