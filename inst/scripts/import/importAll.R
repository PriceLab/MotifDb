library(org.Ce.eg.db)
library(org.Dm.eg.db)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Sc.sgd.db)
# biocLite(c("org.Ce.eg.db", "org.Dm.eg.db", "org.Hs.eg.db", "org.Mm.eg.db", "org.Sc.sgd.db"))

directories <- c("flyFactorSurvey", "hPDI", "jaspar", "ScerTF", "stamlab",
                 "uniprobe", "jaspar2014", "jolma2013", "cisbp")
                 #"HOCOMOCO",

starting.directory <- getwd()
stopifnot(basename(starting.directory) == "import")

#repoRoot <- "~/s/data/public/TFBS"
#repoRoot <- "/shared/silo_researcher/Morgan_M/BioC/MotifDb"
repoRoot <- "/fh/fast/morgan_m/BioC/MotifDb-raw-data"

for(directory in directories){
    print(noquote(sprintf("--- importing %s", directory)))
    setwd(file.path(starting.directory, directory))
    #source("test.R")
    #run.tests(repoRoot)
    source("import.R")
    run(repoRoot)
    }
    
setwd(starting.directory)
