directories <- c("flyFactorSurvey", "hPDI", "jaspar", "ScerTF", "stamlab")#   "uniprobe"       
starting.directory <- getwd()
stopifnot(basename(starting.directory) == "import")

kDataDir <- "~/s/data/public/TFBS"
kDataDir <- "/shared/silo_researcher/Morgan_M/BioC/MotifDb"

for(directory in directories){
    print(noquote(sprintf("--- importing %s", directory)))
    setwd(file.path(starting.directory, directory))
    source("test.R")
    run.tests(kDataDir)
    source("import.R")
    run(kDataDir)
    }
    
