library(MotifDb)
library(trena)   # only for MotifMatcher
library(igvR)
library(chipDB)

#------------------------------------------------------------------------------------------------------------------------
targetGene <- "GATA2"
tf <- "ZNF263"
motif <- query(MotifDb, c("sapiens", tf, "jaspar2018"))
stopifnot(length(motif) == 1)

igv <- igvR()
setGenome(igv, "hg38")
showGenomicRegion(igv, "GATA2")
getGenomicRegion(igv)

tbl.region <- with(getGenomicRegion(igv), data.frame(chrom=chrom, start=start, end=end, stringsAsFactors=FALSE))

matcher <- MultiMethodMotifMatcher("hg38", as.list(motif), tbl.region, "Biostrings matchPWM", .80)


tbl.enhancers <- get(load(system.file(package="TrenaProject", "extdata", "genomeAnnotation", "geneHancer.v4.7.allGenes.RData")))
tbl.targetGeneEnhancers <- subset(tbl.enhancers, geneSymbol==targetGene)
track <- DataFrameQuantitativeTrack("enhancers", tbl.targetGeneEnhancers[, c("chrom", "start", "end", "combinedScore")],
                                    autoscale=FALSE, min=0, max=50, color="brown")
displayTrack(igv, track)
showGenomicRegion(igv, "chr3:128,470,539-128,502,070")


motifMatcher <- MotifMatcher("hg38", as.list(motif))
tbl.out <- findMatchesByChromosomalRegion(motifMatcher, tbl.region, pwmMatchMinimumAsPercentage=80)[, c(2,3,4,7)]

track <- DataFrameQuantitativeTrack("moods-ZNF263", tbl.out, autoscale=TRUE, color="blue")
displayTrack(igv, track)

cdb <- chipDB(quiet=FALSE)
tbl.cdb <- with(getGenomicRegion(igv), getHits(cdb, chrom, start, end))
tbl.csHits <- subset(tbl.cdb, tf=="ZNF263")
track <- DataFrameAnnotationTrack(sprintf("%s-chip", tf), tbl.csHits[,c(1,2,3,5)], color="red")
displayTrack(igv, track)


