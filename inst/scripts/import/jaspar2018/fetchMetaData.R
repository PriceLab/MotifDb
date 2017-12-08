# Get meta data for JASPAR 2018
library(jsonlite)
library(dplyr)

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

# Bring in a single matrix with its metadata
getMatrixMetaData <- function(matrix.id){
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
    df.row <- as_data_frame(matrix(unlist(short.result), nrow = 1)) 
    names(df.row) <- desired.cols

    return(df.row)
}

# Test with 25 matrices
my.motifs <- head(matrices$matrix_id, 100)
system.time(data.list <- lapply(my.motifs, getMatrixMetaData)) # About 7.8 seconds

# Test with parallel process
library(BiocParallel)
register(MulticoreParam(workers = 5))
system.time(data.list <- bplapply(my.motifs, getMatrixMetaData)) # About 2.2 seconds

# Do it with sleep and a couple small chunks
all.metadata <- data.frame()
chunksize <- 25
my.motifs <- matrices$matrix_id
for(i in 1:floor(length(my.motifs)/chunksize)){
  
  # Make the chunk
  chunk.start <- chunksize*i - (chunksize-1)
  chunk.end <- chunksize*i
  chunk <- my.motifs[chunk.start:chunk.end]
  
  # Run in parallel
  data.list <- bplapply(chunk, getMatrixMetaData)
  all.metadata <- bind_rows(all.metadata, data.list)
  
  # Sleep for 2 seconds
  Sys.sleep(2)
}

# Do the last chunk
last.index <- floor(length(my.motifs)/chunksize) * chunksize
if(last.index < length(my.motifs)){
  final.chunk <- my.motifs[(last.index + 1):length(my.motifs)]
  data.list <- bplapply(final.chunk, getMatrixMetaData)
  all.metadata <- bind_rows(all.metadata, data.list)
}
# Save the data frame
save(all.metadata, file = "./metaData.Rdata")
