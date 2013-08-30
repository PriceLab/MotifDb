# jolma2013/import.R
#------------------------------------------------------------------------------------------------------------------------
library (org.Hs.eg.db)
library (org.Mm.eg.db)
#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
run = function (dataDir)
{
   dataDir <- file.path(dataDir, "jolma2013")
   x <- readRawMatrices(dataDir)
   matrices <- x$matrices
   tbl.mAnno <- x$anno
   tbl.mAnno$symbol[tbl.mAnno$symbol=="Tp53"] <- "Trp53"
   tbl.mAnno$symbol[tbl.mAnno$symbol=="Tp73"] <- "Trp73"

   tbl.anno <- readRawAnnotationTable(dataDir)
   mouse.syms <- subset(tbl.anno, organism=="mouse")$symbol
   organism <- rep("Hsapiens", nrow(tbl.mAnno))
   mouse.indices <- which(names(matrices) %in% mouse.syms)
   organism[mouse.indices] <- "Mmusculus"
   tbl.mAnno <- cbind(tbl.mAnno, organism, stringsAsFactors=FALSE)

   tbl.md <- createMetadataTable(tbl.mAnno, matrices)
   rownames(tbl.md) <- tbl.md$providerName

   matrices <- normalizeMatrices(matrices)
   names(matrices) <- tbl.md$providerName

   serializedFile <- "jolma2013.RData"
   save (matrices, tbl.md, file=serializedFile)
   printf("saved %d matrices to %s", length(matrices), serializedFile)
   printf("next step: copy %s to <packageRoot>/MotifDb/inst/extdata, and rebuild package", serializedFile)


   browser("rr")

} # run
#-------------------------------------------------------------------------------
normalizeMatrices = function (matrices)
{
  mtx.normalized = sapply (matrices, function (mtx) apply (mtx, 2, function (colvector) colvector / sum (colvector)))
  invisible (mtx.normalized)

} # normalizeMatrices
#------------------------------------------------------------------------------------------------------------------------
createMetadataTable = function (tbl.mAnno, matrices)
{
    options (stringsAsFactors=FALSE)
    tbl.md = data.frame ()
    dataSource = 'Jolma2013'

    provider.names <- createProviderNames(tbl.mAnno)

    stopifnot(length(which(duplicated(provider.names))) == 0)

    stopifnot(nrow(tbl.mAnno) == length(matrices))

    size <- length(provider.names)
    tbl <- data.frame(list(providerName=vector("character", size),
                           providerId=vector("character", size),
                           dataSource=vector("character",size),
                           geneSymbol=vector("character",size),
                           geneId=vector("character",size),
                           geneIdType=vector("character",size),
                           proteinId=vector("character",size),
                           proteinIdType=vector("character",size),
                           organism=vector("character",size),
                           sequenceCount=vector("integer", size),
                           bindingSequence=vector("character",size),
                           bindingDomain=vector("character",size),
                           tfFamily=vector("character",size),
                           experimentType=vector("character",size),
                           pubmedID=vector("character",size)))

    #rownames(tbl) <- provider.names
                           
    geneIDs <- determineGeneIDs(tbl.mAnno$symbol)
    counts = as.integer (sapply (matrices, function (mtx) as.integer (mean (round (colSums (mtx))))))

    for (r in 1:length(provider.names)){
        tbl[r, "providerName"] <- provider.names[r]
        tbl[r, "providerId"]   <- tbl.mAnno$symbol[r]
        tbl[r, "dataSource"]   <- "jolma2013"
        tbl[r, "geneSymbol"]   <- tbl.mAnno$symbol[r]
        tbl[r, "geneId"]   <- geneIDs[r]
        tbl[r, "geneIdType"]   <- "ENTREZ"
        tbl[r, "proteinId"]   <- NA_character_
        tbl[r, "proteinIdType"]   <- NA_character_
        tbl[r, "organism"]   <-  tbl.mAnno$organism[r]
        tbl[r, "sequenceCount"]   <- counts[r]
        tbl[r, "bindingSequence"]   <- tbl.mAnno$sequence[r]
        tbl[r, "bindingDomain"]   <- NA_character_
        tbl[r, "tfFamily"]   <-  tbl.mAnno$family[r]
        tbl[r, "experimentType"]   <- "SELEX"
        tbl[r, "pubmedID"]   <- "23332764"
        } # for r

   invisible (tbl)

} # createMetadataTable
##------------------------------------------------------------------------------------------------------------------------
readRawMatrices = function (dataDir)
{
    data.file <- file.path(dataDir, "TableS3-PWM-models.tsv")
    tbl.raw <- read.table(data.file, skip=15, header=FALSE, as.is=TRUE, fill=TRUE, sep="\t")

     # the tsv file, exported from the excel spreadsheet, \r converted to \n,
     # has an unconventional structure, with 5 rows per pwm, including 1 header line
     # followed by by 4 rows of base counts, in which the first column is A,C,T or G
     # each header line has 10 columns, apparently these

     # symbol family clone         type ligand          sequence batch Seed multinomial cycle
     #  BCL6B   C2H2   DBD  TGCGGG20NGA     AC TGCTTTCTAGGAATTMM     2    4   monomeric   yes
     #   CTCF   C2H2  full TAGCGA20NGCT     AJ NGCGCCMYCTAGYGGTN     2    4   monomeric   yes
     #   EGR1   C2H2   DBD  TCTCTT20NGA      Y    NMCGCCCMCGCANN     2    2   monomeric   yes
     #   EGR1   C2H2  full TACTAT20NATC     AA    NACGCCCACGCANN     2    4   monomeric    no
     #   EGR2   C2H2   DBD  TCGGCC20NGA      W       MCGCCCACGCA     2    3   monomeric    no

     # though there are some left-over column titles.  don't know quite what to make of them
     #         site:
     #         type:
     #  comment:
     #  Matrix is one of the representative PWMs:
  
    
    matrix.start.rows <- seq(3, nrow(tbl.raw), 5)
    tbl.anno <- tbl.raw[matrix.start.rows, 1:10]
    colnames(tbl.anno) <- c("symbol", "family", "clone", "type", "ligand",
                            "sequence", "batch", "Seed", "multinomial",
                            "cycle")
    tbl.indices <- data.frame(start = matrix.start.rows + 1, end = matrix.start.rows + 4)
    matrix.list = apply(tbl.indices, 1, function(x) tbl.raw[x[1]:x[2],])
    names(matrix.list) = tbl.anno$symbol

    matrix.list = lapply(matrix.list, function(mtx) {
       tmp = mtx[,-1]   # drop the ACTG column
          # identify empty or NA columns
       columns.to.discard = apply(tmp, 2, function(x) all(is.na(x) | all(x == "")))
       tmp = tmp[, ! columns.to.discard]
       tmp2 = as.numeric(unlist(tmp))
       tmp = matrix(tmp2, nrow = 4, ncol = length(tmp2)/4)
       rownames(tmp) = mtx[,1]
       colnames(tmp) = 1:ncol(tmp)
       tmp
       })

    unique.ids <- apply(tbl.anno, 1, function(s) paste(s, collapse="."))
    #names(matrix.list) <- unique.ids

       # insert some naming repairs discovered by diego diez
    
    tbl.anno$symbol[tbl.anno$symbol=="Tp53"] <- "Trp53"
    tbl.anno$symbol[tbl.anno$symbol=="Tp73"] <- "Trp73"

    invisible (list(matrices=matrix.list, anno=tbl.anno))

} # readRawMatrices
#------------------------------------------------------------------------------------------------------------------------
readRawAnnotationTable <- function(dataDir)
{
    data.file <- file.path(dataDir, "Cell_2013_Jolma_annotations.tsv")
    tbl.raw <- read.table(data.file, header=TRUE, as.is=TRUE, fill=TRUE, sep="\t")
        # no need for the protein sequence
    tbl.raw <- tbl.raw[,-3]

       # insert some naming repairs discovered by diego diez
    
    tbl.raw$symbol[tbl.raw$symbol=="Egr1_E410D_FARSDERtoFARSDDR"] <- "Egr1"

    invisible(tbl.raw)
    
} # readRawAnnotationTable
#------------------------------------------------------------------------------------------------------------------------
determineGeneIDs <- function(symbols)
{
    geneIDs <- mget(symbols, org.Hs.egSYMBOL2EG, ifnotfound=NA)
    failures <- unique(names(which(is.na(geneIDs))))
    geneIDs <- geneIDs[which(!is.na(geneIDs))]

    alias.failures <- c()
    
    if(length(failures) > 0){ #
       aliasGeneIDs <- mget(failures, org.Hs.egALIAS2EG, ifnotfound=NA)
       alias.failures <- names(which(is.na(aliasGeneIDs)))
       alias.success <- aliasGeneIDs[which(!is.na(aliasGeneIDs))]
       if(length(alias.success) > 0)
          geneIDs <- c(geneIDs, alias.success)
       if(length(alias.failures) > 0){
          mouseGeneIDs <- mget(alias.failures, org.Mm.egSYMBOL2EG, ifnotfound=NA)
          mouse.failures <- names(which(is.na(mouseGeneIDs)))
          mouse.success <- mouseGeneIDs[which(!is.na(aliasGeneIDs))]
          if(length(mouse.success) > 0)
             geneIDs <- c(geneIDs, mouse.success)
          if(length(mouse.failures) > 0){
             names(mouse.failures) <- mouse.failures
             mouse.failures [mouse.failures] <- NA
             geneIDs <- c(geneIDs, mouse.failures)
             } # if mouse.failures
          } # if alias.failures
       }# if failures
       
    result <- as.character(sapply (symbols,
                     function(sym) {target <- sprintf("^%s$", toupper(sym));
                                    as.character(geneIDs[grep(target, names(geneIDs))])[1]}))
    names(result) <- symbols
    

    result

} # determineGeneIDs
#-------------------------------------------------------------------------------
createProviderNames <- function(tbl.mAnno)
{
    x <- as.character(apply(tbl.mAnno, 1,
                      function(row) paste(c(row[11], "jolma2013", row[c(1)]), collapse="-")))

    dups <- which(duplicated(x))
    # printf("incoming dups count: %d", length(dups))
    
    suffix <- 1
    
    while(length(dups) > 0) {
       suffix <- suffix + 1
       fixed.names <- as.character(sapply(x[dups],
                              function(dup) paste(c(dup, suffix), collapse="-")))
       x[dups] <- fixed.names
       dups <- which(duplicated(x))
       } # while dups

     x <- sub("2$", "2", x)
     x <- sub("2-3$", "3", x)
     x <- sub("2-3-4$", "4", x)
     x <- sub("2-3-4-5$", "5", x)
     x <- sub("2-3-4-5-6$", "6", x)
     x <- sub("2-3-4-5-6-7$", "7", x)
     x <- sub("2-3-4-5-6-7-8$", "8", x)
     x <- sub("2-3-4-5-6-7-8-9$", "9", x)

     x

} # createProviderNames
#-------------------------------------------------------------------------------
