# Get metadata for JASPAR 2022
library(jsonlite)
library(RUnit)
library(BiocParallel)
#----------------------------------------------------------------------------------------------------
getIDs <- function()
{
    f <- "JASPAR2022_CORE_non-redundant_pfms_jaspar.txt"
    cmd <- sprintf("grep '>' %s", f)
    ids.raw <- system(cmd, intern=TRUE)
    length(ids.raw)  # 1956
    ids.cooked <- sub(">", "", ids.raw)
    tokens <- strsplit(ids.cooked, "\t")
    ids <- unlist(lapply(tokens, "[", 1))
    return(ids)

} # getIDs
#----------------------------------------------------------------------------------------------------
ncbiTaxonimicCodeToBiocLinnaean <- function(code)
{
  code <- as.character(code)

  lookup <- list("3702" = "Athaliana",
                 "3888" = "Psativum",
                 "4094" = "Nsp.",
                 "4102" = "Phybrida",
                 "4151" = "Amajus",
                 "4513" = "Hvulgare",
                 "4565" = "Taestivam",
                 "4577" = "Zmays",
                 "4932" = "Scerevisiae",
                 "6239" = "Celegans",
                 "7227" = "Dmelanogaster",
                 "7729" = "Hroretzi",
                 "7742" = "Vertebrata",
                 "8022" = "Omykiss",
                 "8355" = "Xlaevis",
                 "8364" = "Stropicalis",
                 "9031" = "Ggallus",
                 "9606" = "Hsapiens",
                 "9913" = "Btaurus",
                 "9986" = "Ocuniculus",
                 "10090" = "Mmusculus",
                 "10116" = "Rnorvegicus",
                 "10117" = "Rrattus")

  if (code %in% names(lookup))
    return(lookup[[code]])

  NA

} # ncbiTaxonimicCodeToLinnaean
#----------------------------------------------------------------------------------------------------
# Bring in a single matrix with its metadata
getMatrixMetaData <- function(matrix.id, verbose=FALSE)
{
    #browser()
    if(verbose) message(sprintf("%s", matrix.id))
    url <- sprintf("http://jaspar.genereg.net/api/v1/matrix/%s/", matrix.id)
    result <- fromJSON(url)
    result$pfm <- NULL # Get rid of PFM data

    # Split the species and convert the tax ID to the species we want
    result$tax_id <- result$species$tax_id
    result$species <- lapply(result$species$tax_id, ncbiTaxonimicCodeToBiocLinnaean)

    # Shorten the list
    desired.cols <- c("matrix_id", "name", "family", "species", "class", "uniprot_ids", "type", "pubmed_ids")
    short.result <- result[desired.cols]

    # Turn empty lists into "NA"
    idx <- sapply(short.result, function(x) length(x) == 0)
    short.result[idx] <- NA

    # Collapse things with multiple entries
    short.result <- lapply(short.result, paste0, collapse = ";")

    # Turn it into a row in a matrix and then a data.frame
    tbl <- as.data.frame(matrix(unlist(short.result), nrow = 1))
    names(tbl) <- desired.cols

    return(tbl)

} # getMatrixMetaData
#----------------------------------------------------------------------------------------------------
test_getMatrixMetaData <- function()
{
     # Test with 25 matrices
   test.motifs <- head(getIDs(), n=5)

   tbl <- getMatrixMetaData(test.motifs[1], verbose=TRUE)
   checkTrue(is.data.frame(tbl))
   checkEquals(colnames(tbl), c("matrix_id", "name", "family", "species", "class", "uniprot_ids",
                                "type", "pubmed_ids"))
      # make sure no factors
   checkTrue(all(unlist(lapply(tbl, class), use.names=FALSE) == "character"))
   checkEquals(as.character(tbl),
               c("MA0004.1", "Arnt", "PAS domain factors", "Mmusculus",
                 "Basic helix-loop-helix factors (bHLH)", "NA", "SELEX", "7592839"))


   system.time(tbls <- lapply(test.motifs, getMatrixMetaData)) # About 7.8 seconds
   tbl <- do.call(rbind, tbls)
   checkEquals(dim(tbl), c(5, 8))

     # Test with parallel process
   test.motifs <- head(getIDs(), n=10)
   register(MulticoreParam(workers = 3))
   system.time(tbls <- bplapply(test.motifs, getMatrixMetaData)) # About 2.2 seconds
   tbl <- do.call(rbind, tbls)
   checkEquals(dim(tbl), c(10, 8))

} # test_getMatrixMetaData
#----------------------------------------------------------------------------------------------------
getAll <- function()
{
   motif.ids <- getIDs()
   length(motif.ids)

   register(MulticoreParam(workers = 4))
   tbls <- bplapply(motif.ids, function(id) getMatrixMetaData(id, verbose=TRUE))
   #tbls <- lapply(motif.ids, function(id) getMatrixMetaData(id, verbose=TRUE))
   tbl.md <- do.call(rbind, tbls)
   save(tbl.md, file="metaData.RData")

} # getAll
#----------------------------------------------------------------------------------------------------
