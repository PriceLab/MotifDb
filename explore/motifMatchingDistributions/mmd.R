library(MotifDb)
library(TrenaViz)   # only for MultiMethodMotifMatcher
library(igvR)
igv <- igvR()
igvR::setGenome(igv, "hg38")
showGenomicRegion(igv, "GATA2")

library(TrenaProjectErythropoiesis)
tp <- TrenaProjectErythropoiesis()
setTargetGene(tp, "GATA2")

tbl.enhancers <- getEnhancers(tp)[, c(1,2,3,5)]

track <- DataFrameQuantitativeTrack("GH", tbl.enhancers,
                                    autoscale=FALSE, min=0, max=50, color="brown")
displayTrack(igv, track)
tbl.bigRegion <- with(tbl.enhancers, data.frame(chrom=tbl.enhancers$chrom[1], start=min(start)-5000, end=max(end)+5000),
                      stringsAsFactors=FALSE)

bigRegion <- with(tbl.enhancers, sprintf("%s:%d-%d", chrom=tbl.enhancers$chrom[1], start=min(start)-5000, end=max(end)+5000))
showGenomicRegion(igv, bigRegion)

motif.tbx15 <- query(MotifDb, c("TBX15", "sapiens"), "jaspar2018")
m4.biostrings <- MultiMethodMotifMatcher("hg38", as.list(motif.tbx15), tbl.bigRegion, "Biostrings matchPWM", .80)
tbl.tbx15 <- matchMotifInSequence(m4.biostrings)
dim(tbl.tbx15)
track <- DataFrameQuantitativeTrack("matchPWM-tbx15", tbl.tbx15[, c(1,2,3,6)], autoscale=TRUE, color="blue")
displayTrack(igv, track)


motif.znf263 <- query(MotifDb, c("sapiens", "ZNF263", "jaspar2018"))

tbl.znf263 <- matchMotifInSequence(m4.biostrings)
dim(tbl.znf263)
tbl.znf263$chrom <- as.character(tbl.znf263$chrom)


tbl.tbx15$chrom <- as.character(tbl.tbx15$chrom)
track <- DataFrameQuantitativeTrack("jaspar2018-TBX15-MA0803.1", tbl.tbx15[, c(1,2,3,6)], autoscale=TRUE, color="blue")
displayTrack(igv, track)

m4.moods <- MultiMethodMotifMatcher("hg38", as.list(motif), tbl.bigRegion, "MOODS matchMotifs", 6)
tbl.tbx15.2 <-  matchMotifInSequence(m4.moods)
dim(tbl.tbx15.2)
track <- DataFrameQuantitativeTrack("moods jaspar2018-TBX15-MA0803.1", tbl.tbx15.2[, c(1,2,3,6)], autoscale=TRUE, color="green")
displayTrack(igv, track)


library(FimoClient)

FIMO_HOST <- "localhost"
FIMO_PORT <- 600161
if(!exists("fimoServerStarted")){
   fimoServerStarted <- TRUE
   export(motif, con="motif.meme", format="meme")
   cmd <- sprintf("make -f ~/github/fimoService/server/makefile", PORT=%d MOTIFS=motif.meme", FIMO_PORT)
   system(cmd)
   #meme.file <- system.file(package="FimoClient", "extdata", "human.jaspar2018.meme")
   #stopifnot(file.exists(meme.file))
   #cmd <- sprintf("make -f ~/github/fimoService/server/makefile PORT=%d MOTIFS=%s", FIMO_PORT, meme.file)
   #print(cmd)
   #system(cmd)
   printf("--- sleeping 5, making sure fimo server is awake")
   Sys.sleep(5)
   }


fc <- FimoClient(FIMO_HOST, FIMO_PORT, quiet=FALSE)
tbl.fimo <- requestMatchForRegions(fc, tbl.bigRegion, "hg38", pvalThreshold=0.00006)
track <- DataFrameQuantitativeTrack("fimo jaspar2018-TBX15-MA0803.1", tbl.fimo[, c(1,2,3,6)], autoscale=TRUE, color="darkred")
displayTrack(igv, track)

tbl.fimo.qScore <- tbl.fimo[, c(1,2,3,8)]
tbl.fimo.qScore$qValue <- -log10(tbl.fimo.qScore$qValue)
track <- DataFrameQuantitativeTrack("fimo qScore TBX15", tbl.fimo.qScore[, c(1,2,3,4)], autoscale=TRUE, color="orange")
displayTrack(igv, track)


library(chipDB)
cdb <- chipDB(quiet=FALSE)
tbl.cdb <- with(getGenomicRegion(igv), getHits(cdb, chrom, start, end))

tbl.tf <- as.data.frame(table(tbl.cdb$tf))
colnames(tbl.tf) <- c("tf", "count")
tbl.tf <- tbl.tf[order(tbl.tf$count, decreasing=TRUE),]
tfs <- unique(tbl.cdb$tf)
length(tfs)


genesInModelForMarjorie <- c("FOXM1", "RELA", "GATA1", "PAX6", "ATF2", "GTF2I", "ELK3", "SREBF2", "TBX15", "VDR", "PATZ1",
                             "SIN3A", "RFX5", "ZNF263", "BATF", "PLAGL1", "MYB")
length(genesInModelForMarjorie)  # [1] 17
subset(tbl.tf, tf %in% genesInModelForMarjorie)   # 10
 #         tf count
 # 218  SIN3A    44
 # 206   RELA    21
 # 90   GATA1    11
 # 291 ZNF263    10
 # 273    VDR     8
 # 85   FOXM1     5
 # 154    MYB     5
 # 208   RFX5     4
 # 239 SREBF2     2
 # 12    ATF2     1

track <- DataFrameAnnotationTrack("ZNF263", subset(tbl.cdb, tf=="ZNF263")[, c("chrom", "start", "end", "tissueOrCellType")])
displayTrack(igv, track)

motif.znf263 <- query(MotifDb, c("sapiens", "ZNF263", "jaspar2018"))

tbl.region <- with(getGenomicRegion(igv), data.frame(chrom=chrom, start=start, end=end, stringsAsFactors=FALSE))
m4.biostrings <- MultiMethodMotifMatcher("hg38", as.list(motif.znf263), tbl.region, "Biostrings matchPWM", .75)
tbl.znf263 <- matchMotifInSequence(m4.biostrings)
dim(tbl.znf263)
tbl.znf263$chrom <- as.character(tbl.znf263$chrom)
track <- DataFrameQuantitativeTrack("moods-ZNF263", tbl.znf263[, c(1,2,3,6)], autoscale=TRUE, color="blue")
displayTrack(igv, track)
