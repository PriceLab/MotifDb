# jolma2013/test.R
#------------------------------------------------------------------------------------------------------------------------
library (RUnit)
source("import.R")
#------------------------------------------------------------------------------------------------------------------------
# for pshannon development, dataDir <- "~/s/data/public/TFBS"
# for rhino deployment,     dataDir <- ""
run.tests = function (dataDir)
{
    dataDir <- file.path(dataDir, "jolma2013")
    
    x.rmat <<- test.readRawMatrices(dataDir)
    rmat <- x.rmat$matrices
    ranno <- x.rmat$anno
    
    x.anno <<- test.readRawAnnotationTable(dataDir)

    #test.readRaw(dataDir)

       # mapping from matrices to annotation is done by symbol
       # make sure that diego's fixes are in place, and that
       # these correspond.  note, however, that some are duplicated
    
    checkTrue(all(x.rmat$anno$symbol %in% x.anno$symbol))
    checkTrue(all(x.anno$symbol %in% x.rmat$anno$symbol))

    test.createProviderNames(dataDir)
    
} # run.tests
#------------------------------------------------------------------------------------------------------------------------
test.readRaw = function (dataDir)
{
    print ("--- test.readRaw")
    x = readRaw (dataDir)

} # test.readRaw
#------------------------------------------------------------------------------------------------------------------------
test.readRawAnnotationTable = function (dataDir)
{
    print("--- test.readRawAnnotationTable")
    tbl.annoRaw <- readRawAnnotationTable(dataDir)
    checkEquals(nrow(tbl.annoRaw), 538)
    checkEquals(colnames(tbl.annoRaw),
                c("symbol", "ensembl", "family", "organism", "clone_type", "method", "clone_source"))
    checkEquals(as.list(table(tbl.annoRaw$organism)), list(human=454, mouse=84))
  
    print(noquote("TRUE"))
    invisible(tbl.annoRaw)

} # test.readRawMatrices
#------------------------------------------------------------------------------------------------------------------------
test.readRawMatrices = function (dataDir)
{
    print ("--- test.readRawMatrices")
    x = readRawMatrices (dataDir)
    checkEquals(names(x), c("matrices", "anno"))
    checkEquals(length(x$matrices), 843)
    checkEquals(nrow(x$anno), 843)
    #checkEquals(names(x$matrices)[1],
    #            "BCL6B.C2H2.DBD.TGCGGG20NGA.AC.TGCTTTCTAGGAATTMM.2.4.monomeric.")
  
    print(noquote("TRUE"))
    invisible(x)

} # test.readRawMatrices
#------------------------------------------------------------------------------------------------------------------------
test.readRawAnnotationTable = function (dataDir)
{
    print("--- test.readRawAnnotationTable")
    tbl.annoRaw <- readRawAnnotationTable(dataDir)
    checkEquals(nrow(tbl.annoRaw), 538)
    checkEquals(colnames(tbl.annoRaw),
                c("symbol", "ensembl", "family", "organism", "clone_type", "method", "clone_source"))
    checkEquals(as.list(table(tbl.annoRaw$organism)), list(human=454, mouse=84))
  
    print(noquote("TRUE"))
    invisible(tbl.annoRaw)

} # test.readRawMatrices
#------------------------------------------------------------------------------------------------------------------------
test.convertRawMatricesToStandard = function (matrices)
{
    print("--- test.convertRawMatricesToStandard")
    convertRawMatricesToStandard(matrices[1:2])
    browser("t.cRMTS")
    
} # test.convertRawMatricesToStandard 
#------------------------------------------------------------------------------------------------------------------------
test.determineGeneIDs <- function()
{
     print("--- test.determineGeneIDs")
     aliases.good <- c("ZNF238", "ZNF238", "ZNF306", "POU5F1P1",  "BHLHB2")
     symbols.good <- c("BCL6B", "Ctcf", "EGR1", "EGR1",  "EGR2", "EGR2", "SOX1")
     symbols.unmappable <- c("HINFP1","Hinfp1","HINFP1","ZFP652","ZFP740",
                             "TCFAP2A","TCFAP2A","TCFAP2A","Msx3","MSX3","RHOX11")
     syms <- c(symbols.good, symbols.unmappable, aliases.good)
     geneIDs <- determineGeneIDs(syms)
     checkEquals(length(geneIDs), length(syms))
     checkEquals(names(geneIDs), syms)
     checkEquals(geneIDs[["SOX1"]], "6656")
     checkEquals(geneIDs[["BCL6B"]], "255877")
     checkEquals(geneIDs[["POU5F1P1"]], "5462")
     checkEquals(geneIDs[["BHLHB2"]], "8553")

     indices <- grep("HINFP1", names(geneIDs))
     checkTrue(all(is.na(geneIDs[indices])))
     indices <- grep("TCFAP2A", names(geneIDs))
     checkTrue(all(is.na(geneIDs[indices])))

     browser("mouse check")

     mouse.syms <- c("Elf5", "Elk3", "Spic", "Ascl2", "Atoh1",
                     "Bhlhb2", "Mlx", "Srebf1", "Tcf21", "Atf4", "Cebpb", "Creb3l2",
                     "Creb5", "Dbp", "Hlf", "Jdp2", "Mafb", "Alx1", "Alx4", "Arx",
                     "Barhl1", "Dlx1", "Dlx2", "En2", "Gbx1", "Gbx2", "HMBOX1",
                     "Hoxa11", "Hoxa2", "Hoxc10", "Hoxd13", "Hoxd3", "Hoxd9", "Irx3",
                     "Lhx4", "Lhx8", "Meis2", "Meis3", "Meox2", "Msx3", "Nkx3-1",
                     "Nkx6-1", "Otx1", "Pknox2", "Prrx2", "Rhox11", "Shox2", "Uncx",
                     "Vsx1", "Pou2f2", "Rfx2", "Rfx3", "Tcfap2a", "Egr1", "Egr3",
                     "Hic1", "Klf12", "Zfp652", "Zfp740", "Zic3", "Sox1", "Sox10",
                     "Sox11", "Sox17", "Sox3", "Tcf7", "Foxc1", "Foxg1", "Foxj3",
                     "Foxk1", "Trp53", "Trp73", "TBR1", "Ar", "Esrra", "Hnf4a",
                     "Nr2e1", "Nr2f6", "Rara", "Rarb", "Rarg", "Rxra", "Rxrb", "Vdr")
                 
     #indices <- grep("EGR1", names(geneIDs))
     #checkTrue(all(geneIDs[indices] == "1958"))


} # test.determineGeneIDs
#-------------------------------------------------------------------------------
test.createProviderNames <- function(dataDir)
{
   print("--- test.createProviderNames")
   x <- readRawMatrices(dataDir)
   tbl.mAnno <- x$anno
   provider.names <- createProviderNames(tbl.mAnno)

   checkEquals(length(which(duplicated(provider.names))), 0)
   
} # test.createProviderNames
#-------------------------------------------------------------------------------
