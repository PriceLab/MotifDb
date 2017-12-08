library (MotifDb)
library (RUnit)
library (MotIV)
library (seqLogo)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
# Note: all matrix numbers were chosen randomly using "sample()"
runTests = function ()
{

} # runTests
#----------------------------------------------------------------------------------------------------
test_cisbp <- function(){

    # Load the cisbp matrix as "cisbp"
    load("./single_matrices/cisbp_matrix.Rdata")
    
    # Matrix name
    mtx.name <- "M0308_1.02"

    # Query for the same matrix
    queried <- query(MotifDb, mtx.name)[[1]]

    # Check that they're the same
    checkTrue(all(cisbp == queried))
    
} # test_cisbp
#----------------------------------------------------------------------------------------------------
test_FlyFactorSurvey <- function(){

    # Load the cisbp matrix as "cisbp"
    load("./single_matrices/FlyFactorSurvey.Rdata")
    
    # Matrix name
    mtx.name <- "M0308_1.02"

    # Query for the same matrix
    queried <- query(MotifDb, mtx.name)[[1]]

    # Check that they're the same
    checkTrue(all(cisbp == queried))
    
} # test_FlyFactorSurvey
#----------------------------------------------------------------------------------------------------
test_HOCOMOCOv10 <- function(){

    # Load the cisbp matrix as "cisbp"
    load("./single_matrices/HOCOMOCOv10.Rdata")
    
    # Matrix name
    mtx.name <- "M0308_1.02"

    # Query for the same matrix
    queried <- query(MotifDb, mtx.name)[[1]]

    # Check that they're the same
    checkTrue(all(cisbp == queried))
    
} # test_HOCOMOCOv10
#----------------------------------------------------------------------------------------------------
test_HOMER <- function(){

    # Load the cisbp matrix as "cisbp"
    load("./single_matrices/HOMER.Rdata")
    
    # Matrix name
    mtx.name <- "M0308_1.02"

    # Query for the same matrix
    queried <- query(MotifDb, mtx.name)[[1]]

    # Check that they're the same
    checkTrue(all(cisbp == queried))
    
} # test_HOMER
#----------------------------------------------------------------------------------------------------
test_hPDI <- function(){

    # Load the cisbp matrix as "cisbp"
    load("./single_matrices/hPDI.Rdata")
    
    # Matrix name
    mtx.name <- "M0308_1.02"

    # Query for the same matrix
    queried <- query(MotifDb, mtx.name)[[1]]

    # Check that they're the same
    checkTrue(all(cisbp == queried))
    
} # test_hPDI
#----------------------------------------------------------------------------------------------------
test_JASPAR_2014 <- function(){

    # Load the cisbp matrix as "cisbp"
    load("./single_matrices/JASPAR_2014.Rdata")
    
    # Matrix name
    mtx.name <- "M0308_1.02"

    # Query for the same matrix
    queried <- query(MotifDb, mtx.name)[[1]]

    # Check that they're the same
    checkTrue(all(cisbp == queried))
    
} # test_JASPAR_2014
#----------------------------------------------------------------------------------------------------
test_JASPAR_CORE<- function(){

    # Load the cisbp matrix as "cisbp"
    load("./single_matrices/JASPAR_CORE.Rdata")
    
    # Matrix name
    mtx.name <- "M0308_1.02"

    # Query for the same matrix
    queried <- query(MotifDb, mtx.name)[[1]]

    # Check that they're the same
    checkTrue(all(cisbp == queried))
    
} # test_JASPAR_CORE
#----------------------------------------------------------------------------------------------------
test_jaspar2016 <- function(){

    # Load the cisbp matrix as "cisbp"
    load("./single_matrices/jaspar2016.Rdata")
    
    # Matrix name
    mtx.name <- "M0308_1.02"

    # Query for the same matrix
    queried <- query(MotifDb, mtx.name)[[1]]

    # Check that they're the same
    checkTrue(all(cisbp == queried))
    
} # test_jaspar2016
#----------------------------------------------------------------------------------------------------
test_jolma2013 <- function(){

    # Load the cisbp matrix as "cisbp"
    load("./single_matrices/jolma2013.Rdata")
    
    # Matrix name
    mtx.name <- "M0308_1.02"

    # Query for the same matrix
    queried <- query(MotifDb, mtx.name)[[1]]

    # Check that they're the same
    checkTrue(all(cisbp == queried))
    
} # test_jolma2013
#----------------------------------------------------------------------------------------------------
test_ScerTF <- function(){

    # Load the cisbp matrix as "cisbp"
    load("./single_matrices/ScerTF.Rdata")
    
    # Matrix name
    mtx.name <- "M0308_1.02"

    # Query for the same matrix
    queried <- query(MotifDb, mtx.name)[[1]]

    # Check that they're the same
    checkTrue(all(cisbp == queried))
    
} # test_ScerTF
#----------------------------------------------------------------------------------------------------
test_stamlab <- function(){

    # Load the cisbp matrix as "cisbp"
    load("./single_matrices/stamlab.Rdata")
    
    # Matrix name
    mtx.name <- "M0308_1.02"

    # Query for the same matrix
    queried <- query(MotifDb, mtx.name)[[1]]

    # Check that they're the same
    checkTrue(all(cisbp == queried))
    
} # test_stamlab
#----------------------------------------------------------------------------------------------------
test_SwissRegulon <- function(){

    # Load the cisbp matrix as "cisbp"
    load("./single_matrices/SR_matrix.Rdata")
    
    # Matrix name
    mtx.name <- "M0308_1.02"

    # Query for the same matrix
    queried <- query(MotifDb, mtx.name)[[1]]

    # Check that they're the same
    checkTrue(all(cisbp == queried))
    
} # test_SwissRegulon
#----------------------------------------------------------------------------------------------------
test_UniPROBE <- function(){

    # Load the uniprobe matrix as "uniprobe"
    load("./single_matrices/UniPROBE_matrix.Rdata")
    
    # Matrix name
    mtx.name <- "UP00230"

    # Query for the same matrix
    queried <- query(MotifDb, mtx.name)[[1]]

    # Check that they're the same
    checkTrue(all(uniprobe == queried))
    
} # test_UniPROBE

