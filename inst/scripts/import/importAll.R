directories <- c("flyFactorSurvey", "hPDI", "jaspar", "ScerTF", "stamlab")#   "uniprobe"       
starting.directory <- getwd()
stopifnot(basename(starting.directory) == "import")

for(directory in directories){
    print(noquote(sprintf("--- importing %s", directory)))
    setwd(file.path(starting.directory, directory))
    source("test.R")
    run.tests()
    source("import.R")
    run()
    }
    
