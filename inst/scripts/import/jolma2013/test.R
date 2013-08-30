# jolma2013/test.R
#------------------------------------------------------------------------------------------------------------------------
library (RUnit)
source("import.R")
#------------------------------------------------------------------------------------------------------------------------
run.tests = function (dataDir)
{
    #test.determineGeneIDs()
    dataDir <- file.path(dataDir, "jolma2013")

    x.rmat <<- test.readRawMatrices(dataDir)
    rmat <- x.rmat$matrices
    ranno <- x.rmat$anno
    
    x.anno <<- test.readRawAnnotationTable(dataDir)

    test.readRaw(dataDir)

       # mapping from matrices to annotation is done by symbol
       # make sure that diego's fixes are in place, and that
       # these correspond.  note, however, that some are duplicated
    
    checkTrue(all(x.rmat$anno$symbol %in% x.anno$symbol))
    checkTrue(all(x.anno$symbol %in% x.rmat$anno$symbol))

       # some duplicated genes in x.rmat$anno
       # head(sort(table(x.rmat$anno$symbol), decreasing=TRUE))
       #  SOX8    SOX9  BARHL2 CREB3L1   ESRRA   HNF4A 
       #     8       7       6       6       6       6
              
    #browser("rt")

 #  x.tbl.anno <<- test.createAnnotationTable (dataDir)
#  test.assignGeneId (dataDir)
#  x.tbl.md <<- test.createMetadataTable (x.tbl.anno, x.matrices)
#  x.matrices.renamed <<- test.renameMatrices (x.matrices, x.tbl.md, x.tbl.anno)
#  x.matrices.normalized <<- test.normalizeMatrices (x.matrices.renamed)
#
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
##------------------------------------------------------------------------------------------------------------------------
#test.createAnnotationTable = function (dataDir)
#{
#  print ('--- test.createAnnotationTable')
#  tbl.anno = createAnnotationTable (dataDir)
#  checkEquals (dim (tbl.anno), c (513, 13))
#  expected = c ("fullID", "id", "category", "mID", "version", "binder", "speciesID", "proteinID", "family", "tax", "class", "pubmed", "type")
#  checkEquals (colnames (tbl.anno), expected)
#
#  checkEquals (head (tbl.anno$fullID),  c ("MA0001.1", "MA0003.1", "MA0004.1", "MA0005.1", "MA0006.1", "MA0006.1"))
#  invisible (tbl.anno)
#
#} # test.createAnnotationTable
##------------------------------------------------------------------------------------------------------------------------
#test.createMetadataTable = function (tbl.anno, matrices)
#{
#  print ('--- test.createMetadataTable')
#   # try it first with just two matrices
#  tbl.md = createMetadataTable (tbl.anno, matrices [1:2])
#  checkEquals (dim (tbl.md), c (2, 15))
#  checkEquals (colnames (tbl.md), c ("providerName", "providerId", "dataSource", "geneSymbol", "geneId", "geneIdType", 
#                                     "proteinId", "proteinIdType", "organism", "sequenceCount", "bindingSequence",
#                                     "bindingDomain", "tfFamily", "experimentType", "pubmedID"))
#  checkEquals (tbl.md$proteinId, c ('P29383', 'P05549'))
#  checkEquals (tbl.md$proteinIdType, c ('UNIPROT', 'UNIPROT'))
#
#    # now use the whole table
#  tbl.md = createMetadataTable (tbl.anno, matrices)
#  checkEquals (dim (tbl.md), c (length (matrices), 15))
#    # test for proper conversion of speciesID = NA or '-' to Vertebrata
#  checkEquals (which (is.na (tbl.md$organism)), integer (0))
#  checkEquals (grep ('-', tbl.md$organism), integer (0))
#
#    # Mmusculus-JASPAR_CORE-NF-kappaB-MA0061.1 had 'NA' for proteinID, not <NA>. fixed?
#  checkEquals (grep ('NA', tbl.md$proteinId), integer (0))
#  invisible (tbl.md)
#
#} # test.createMetadataTable
##------------------------------------------------------------------------------------------------------------------------
#test.renameMatrices = function (matrices, tbl.md, tbl.anno)
#{
#    # try it with just the first two matrices
#  matrix.pair = matrices [1:2]
#  tbl.md = createMetadataTable (tbl.anno, matrix.pair)
#  checkEquals (dim (tbl.md), c (2, 15))
#  old.matrix.names = names (matrix.pair)
#  matrices.renamed = renameMatrices (matrix.pair, tbl.md)
#
#    # test:  the old name is an id, '9229'.  find, in tbl.anno, the fullID, 'MA0001.1'.  then make sure 'MA000.1' is
#    # in the new name of that same matrix
#
#  for (i in 1:length (matrix.pair)) {
#    fullID = subset (x.tbl.anno, id==old.matrix.names [i])$fullID
#    checkTrue (length (grep (fullID, names (matrices.renamed) [i])) == 1)
#    } # for i
#
#    # now try it for the whole set, with selective focused tests
#
#  tbl.md = createMetadataTable (tbl.anno, matrices)
#  checkEquals (nrow (tbl.md), length (matrices))
#  old.matrix.names = names (matrices)
#  matrices.renamed = renameMatrices (matrices, tbl.md)
# 
#  checkEquals (nrow (tbl.md), length (matrices.renamed))
#  checkEquals (length (grep ('-MA0', names (matrices.renamed))), length (matrices.renamed))
#
#  invisible (matrices.renamed)
#
#} # test.renameMatrices
##------------------------------------------------------------------------------------------------------------------------
#test.convertTaxonCode = function ()
#{
#  print ('--- test.convertTaxonCode')
#
#  checkEquals (convertTaxonCode ('9606'), 'Hsapiens')
#  checkEquals (convertTaxonCode (9606), 'Hsapiens')
#     # anomalous codes, which an examination of the jaspar website reveals as 'vertebrates'
#  checkEquals (convertTaxonCode (NA), 'Vertebrata')
#  checkEquals (convertTaxonCode ('NA'), 'Vertebrata')
#  checkEquals (convertTaxonCode (NA_character_), 'Vertebrata')
#  checkEquals (convertTaxonCode ('-'), 'Vertebrata')
#
#} # test.convertTaxonCode
##------------------------------------------------------------------------------------------------------------------------
#test.guessProteinIdentifierType = function (moleculeName)
#{
#  print ('--- test.guessProteinIdentifierType')
#  checkEquals (guessProteinIdentifierType ('P29383'), 'UNIPROT')
#
#  all.types = sapply (x.tbl.anno$proteinID, guessProteinIdentifierType)
#  checkTrue (length (which (is.na (all.types))) < 12)   # got most of them.
#
#} # test.guessProteinIdentifierType
##------------------------------------------------------------------------------------------------------------------------
#test.normalizeMatrices = function (matrices)
#{
#  print ('--- test.normalizeMatrices')
#
#  colsums = as.integer (sapply (matrices, function (mtx) as.integer (mean (round (colSums (mtx))))))
#  #checkTrue (all (colsums > 1))
#
#  matrices.norm = normalizeMatrices (matrices)
#
#  colsums = as.integer (sapply (matrices.norm, function (mtx) as.integer (mean (round (colSums (mtx))))))
#  checkTrue (all (colsums == 1))
#
#  invisible (matrices.norm)
#
#} # test.normalizeMatrices
##------------------------------------------------------------------------------------------------------------------------
#test.assignGeneId = function (dataDir, proteinId)
#{
#  print ('--- test.assignGeneId')
#  uniprot.ids = c ('Q9GRA5', 'P31314', 'AAC18941', 'O49397')
#  refseq.ids  = c ('NP_995315.1', 'NP_032840', 'NP_599022')
#  yeast.ids   = c ('YKL112W', 'YMR072W', 'YLR131C')
#
#  checkEquals (assignGeneId ('NP_995315.1'), list (geneId='4782', type='ENTREZ'))
#  checkEquals (assignGeneId ('NP_599022'),   list (geneId='6095', type='ENTREZ'))
#
#  checkEquals (assignGeneId ('P31314'),      list (geneId='3195', type='ENTREZ'))
#
#  checkEquals (assignGeneId ('YKL112W'),     list (geneId='YKL112W', type='SGD'))
#
#    # see how successful this is over all 513 proteinIds
#
#  tbl.anno = createAnnotationTable (dataDir)
#  mtx.geneId = as.data.frame (t (sapply (tbl.anno$proteinID, assignGeneId)))
#  tbl.types = as.data.frame (table (as.character (mtx.geneId$type), useNA='always'), stringsAsFactors=FALSE)
#  checkEquals (tbl.types$Var1, c ("ENTREZ", "SGD", NA))
#  checkEquals (tbl.types$Freq, c (142, 177, 194))
#
#} # test.assignGeneId
##------------------------------------------------------------------------------------------------------------------------
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
