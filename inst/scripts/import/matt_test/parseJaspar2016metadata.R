# go.R
#------------------------------------------------------------------------------------------------------------------------
options(stringsAsFactors=FALSE)
library (RUnit)
#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
run = function (levels)
{
  if (0 %in% levels) {
    tbl.matrix <- read.table("MATRIX.txt", sep="\t", header=FALSE, as.is=TRUE)
    colnames(tbl.matrix) <- c("id", "jaspar.class", "ma.name", "unknown", "gene.symbol")
    } # 0

  if (1 %in% levels) {  # a data.frame of 3 cols, multiple rows per id
    tbl.anno <- read.table("MATRIX_ANNOTATION.txt", sep="\t", header=FALSE, as.is=TRUE, quote="")
    colnames(tbl.anno) <- c("id", "attr", "value")
    } # 1

  if (2 %in% levels) {
    tbl.prot <- read.table("MATRIX_PROTEIN.txt", sep="\t", header=FALSE, as.is=TRUE, quote="")
    colnames(tbl.prot) <- c("id", "uniprot")
    } # 2

  if (3 %in% levels) {
    tbl.species <- read.table("MATRIX_SPECIES.txt", sep="\t", header=FALSE, as.is=TRUE, quote="")
    colnames(tbl.species) <- c("id", "ncbi.tax.code")
    } # 3

  if (4 %in% levels) {
    col.names <- unique(tbl.anno$attr)
    row.names <- as.character(unique(tbl.anno$id))
    nas <- rep(NA, length(row.names))
    tbl.a <- data.frame(nas, nas, nas, nas, nas, nas, nas, nas, nas, nas,
                        nas, nas, nas, nas, nas, nas, nas, nas, nas, nas, nas,
                        stringsAsFactors=FALSE)
    colnames(tbl.a) <- col.names
    rownames(tbl.a) <- row.names
    for(matrix.id in row.names){
       # browser()
       # print(matrix.id)
        tbl.sub <- subset(tbl.anno, id==matrix.id)
        # Pull out only unique attributes
        tbl.sub <- tbl.sub[!duplicated(tbl.sub$attr),]
        attrs <- tbl.sub$attr
        vals  <- tbl.sub$value
        tbl.a[matrix.id, attrs] <- vals
        }# for matrix.id
    tbl.a$id <- rownames(tbl.a)
    tbl.md <- tbl.a
    #Pull out only unique ids
  #  browser()
  #  tbl.md <- tbl.md[!duplicated(tbl.md$id),]
    } # 4

  if (5 %in% levels) {
      # "demo" source only uses  first 4 matrices in
      # http://jaspar.binf.ku.dk/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_all.txt
#    moi <- c("MA0004", "MA0006", "MA0008", "MA0009")
#    x <- subset(tbl.matrix, ma.name %in% moi)
#      tbl.md <- merge(merge(merge(x, tbl.prot, by="id"), tbl.species, by="id"), tbl.md, by="id")
      tbl.md <- merge(merge(merge(tbl.matrix,
                                  tbl.prot, by="id", all = TRUE),
                            tbl.species, by="id", all = TRUE),
                      tbl.md, by="id", all = TRUE)
      # This line might cut out some metadata
      tbl.md <- tbl.md[!duplicated(tbl.md$id),]
    printf("saving tbl.md, %d rows, %d cols", nrow(tbl.md), ncol(tbl.md))
    save(tbl.md, file="tbl.md.RData")
    write.table(tbl.md, quote=FALSE, file="md.tsv", sep="\t", row.names=FALSE)
    } # 5

  if (6 %in% levels) {
    } # 6

  if (7 %in% levels) {
    } # 7

  if (8 %in% levels) {
    } # 8

  if (9 %in% levels) {
    } # 9

  if (10 %in% levels) {
    } # 10

  if (11 %in% levels) {
    } # 11

  if (12 %in% levels) {
    } # 12

  if (13 %in% levels) {
    } # 13

  if (14 %in% levels) {
    } # 14

  if (15 %in% levels) {
    } # 15

  if (16 %in% levels) {
    } # 16

  if (17 %in% levels) {
    } # 17

  if (18 %in% levels) {
    } # 18

  if (19 %in% levels) {
    } # 19

  if (20 %in% levels) {
    } # 20


} # run
#------------------------------------------------------------------------------------------------------------------------
