# MotifDb/inst/scripes/import/demo/import.R
#------------------------------------------------------------------------------------------------------------------------
library(RUnit)
options (stringsAsFactors=FALSE)
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
run = function (dataDir)
{
    dataDir <- file.path(dataDir, "jaspar2022")
    stopifnot(file.exists(dataDir))
    rawMatrixList <- readRawMatrices(dataDir)
    checkTrue(length(rawMatrixList) == 1956)
    matrices <- extractMatrices(rawMatrixList)
    checkTrue(length(matrices) == 1956)
    incoming.metadata.file <- file.path(dataDir, "metaData.RData")
    checkTrue(file.exists(incoming.metadata.file))

    tbl.md <- createMetadataTable(dataDir, matrices, metadata.table = "./metaData.RData")
    matrices <- normalizeMatrices(matrices)
    matrices <- renameMatrices(matrices, tbl.md)

    serializedFile <- file.path(dataDir, "jaspar2022.RData")
    printf("writing %s to %s", "jaspar2022.RData", dataDir)

    save(matrices, tbl.md, file=serializedFile)
    printf("saved %d matrices to %s", length(matrices), serializedFile)
    printf("next step:  copy %s to <packageRoot>/MotifDb/inst/extdata, rebuild package", serializedFile)

} # run
#------------------------------------------------------------------------------------------------------------------------
readRawMatrices = function (dataDir)
{
    # our convention is that there is a shared "dataDir" visible to
    # the importer, and that within that directory there is one
    # subdirectory for each data source.
    # for this example importer, that directory will be <dataDir>/test
    # within which we will look for one small file "sample.pcm"


    filename <- file.path(dataDir, "JASPAR2022_CORE_non-redundant_pfms_jaspar.txt")
    printf("checking for readable matrix file:")
    printf("     %s", filename)
    stopifnot(file.exists(filename))

    all.lines = scan(filename, what=character(0), sep='\n', quiet=TRUE)
    title.lines = grep('^>', all.lines)
    title.line.count <<- length(title.lines)
    max = title.line.count -1

    # browser()
    # Matt's additions to fix brackets and delimiters
    all.lines <- gsub("[A,C,G,T]\\s+\\[\\s*","",all.lines)
    all.lines <- gsub("\\s*\\]","",all.lines)
    all.lines <- gsub("\\s+","\t",all.lines)

    pwms = list()

    for (i in 1:max) {
        # print(i)
        start.line = title.lines [i]
        end.line = title.lines [i+1] - 1
        new.pwm = parsePwm(all.lines [start.line:end.line])
        pwms = c(pwms, list(new.pwm))
    } # for i

    # Add the last motif, which got missed
    start.line <- title.lines[title.line.count]
    end.line <- start.line + 4
    new.pwm = parsePwm(all.lines [start.line:end.line])
    pwms = c(pwms, list(new.pwm))


    invisible(pwms)

} # readRawMatrices
#------------------------------------------------------------------------------------------------------------------------
extractMatrices = function (pwm.list)
{
    matrices = sapply(pwm.list, function(element) element$matrix)
    matrix.names <- sapply(pwm.list, function(element) element$title)
    matrix.names <- sub("^>", "", matrix.names)
    names(matrices) <- matrix.names

    matrices

} # extractMatrices
#------------------------------------------------------------------------------------------------------------------------
createMetadataTable = function(dataDir, matrices, metadata.table)
{
    filename <- file.path(dataDir, metadata.table)
    printf("checking for readable metadata file:")
    printf("   %s", filename)
    stopifnot(file.exists(filename))

    tbl.mdRaw <- get(load(filename)) # this is "tbl.mdRaw"

    tbl.md <- data.frame()

    matrix.ids = names(matrices)

    for(matrix.id in matrix.ids) {
        matrix <- matrices[[matrix.id]]
        # Simply subset by the full  matrix ID
        md <- as.list(subset(tbl.mdRaw, matrix_id == matrix.id))
        print(matrix.id)
        stopifnot(length(md$matrix_id) == 1)
        dataSource <- "jaspar2022"
        new.row = list(providerName=matrix.id,
                        providerId=matrix.id,
                        dataSource=dataSource,
                        geneSymbol=md$name,
                        geneId=NA,
                        geneIdType=NA,
                        proteinId=md$uniprot_ids,
                        proteinIdType=guessProteinIdentifierType(md$uniprot_ids),
                        organism=md$species,
                        sequenceCount=max(colSums(matrix)),
                        bindingSequence=NA_character_,
                        bindingDomain=NA,
                        tfFamily=md$family,
                        experimentType=md$type,
                        pubmedID= md$pubmed_ids)
        tbl.md = rbind(tbl.md, data.frame(new.row, stringsAsFactors=FALSE))
        full.name = sprintf('%s-%s-%s-%s',
                             md$species,
                             dataSource,
                             md$name,
                             matrix.id)
        rownames(tbl.md) [nrow(tbl.md)] = full.name
    } # for i

    invisible(tbl.md)

} # createMetadataTable
#------------------------------------------------------------------------------------------------------------------------
test_createMetadataTable <- function()
{

} # test_createMetadataTable
#------------------------------------------------------------------------------------------------------------------------
renameMatrices = function(matrices, tbl.md)
{
    stopifnot(length(matrices) == nrow(tbl.md))
    names(matrices) = rownames(tbl.md)
    invisible(matrices)

} # renameMatrices
#------------------------------------------------------------------------------------------------------------------------
# an empirical and not altogether trustworthy solution to identifying identifier types.
guessProteinIdentifierType = function(moleculeName)
{
    if(is.na(moleculeName))
        return(NA_character_)
    if(nchar(moleculeName) == 0)
        return(NA_character_)

    first.char = substr(moleculeName, 1, 1)

    if(first.char == 'Y')
        return('SGD')

    if(first.char %in% c('P', 'Q', 'O', 'A', 'C'))
        return('UNIPROT')

    if(length(grep('^NP_', moleculeName)) == 1)
        return('REFSEQ')

    return(NA_character_)

} # guessProteinIdentifierType
#------------------------------------------------------------------------------------------------------------------------
normalizeMatrices = function(matrices)
{
    mtx.normalized = sapply(matrices,
                             function(mtx) apply(mtx, 2, function(colvector) colvector / sum(colvector)))

    invisible(mtx.normalized)

} # normalizeMatrices
#------------------------------------------------------------------------------------------------------------------------
parsePwm = function(text)
{
    lines = strsplit(text, '\t')
    stopifnot(length(lines)==5) # title line, one line for each base
    title = lines [[1]][1]
    line.count = length(lines)

    # expect 4 rows, and a number of columns we can discern from
    # the incoming text.
    cols <- length(lines[[2]])
    result <- matrix(nrow=4, ncol=cols,
                      dimnames=list(c('A','C','G','T'),
                                    as.character(1:cols)))
    row = 1
    for(i in 2:line.count){
        result [row,] = as.numeric(lines[[i]])
        row = row + 1
    } # for i

    return(list(title=title, matrix=result))

} # parsePwm
#----------------------------------------------------------------------------------------------------
if(!interactive())
    run("..")
