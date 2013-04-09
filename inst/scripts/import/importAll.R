directories <- c("flyFactorSurvey", "hPDI", "jaspar", "ScerTF", "stamlab", "uniprobe")
starting.directory <- getwd()
stopifnot(basename(starting.directory) == "import")

#repoRoot <- "~/s/data/public/TFBS"
repoRoot <- "/shared/silo_researcher/Morgan_M/BioC/MotifDb"

for(directory in directories){
    print(noquote(sprintf("--- importing %s", directory)))
    setwd(file.path(starting.directory, directory))
    source("test.R")
    run.tests(repoRoot)
    source("import.R")
    run(repoRoot)
    }
    
