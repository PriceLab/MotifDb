MotifDb <- NULL
#-------------------------------------------------------------------------------
.MotifDb = function(loadAllSources=TRUE, quiet=TRUE)
{
  mdb = MotifList()

  if(loadAllSources) {
    data.path = system.file('extdata', package='MotifDb')
    data.files = dir(data.path, full.names=TRUE)
    for(data.file in data.files) {
       # define these to keep 'check' happy.  they are loaded by 'load'
      tbl.md = NA; matrices = NA;  
      variables = load(data.file)
      mdb = append(mdb, MotifList(matrices, tbl.md))
      if(!quiet) 
        message(noquote(sprintf('added %s(%d) matrices, length now: %d',
                 basename(data.file), length(matrices), length(mdb))))
    } # for data.file
  
    if(!quiet) {
      print(table(values(mdb)$dataSource))
      }
    } # if loadAllSources

   message('See system.file("LICENSE", package="MotifDb") for use restrictions.')
  
  return(mdb)

} # MotifDb
#-------------------------------------------------------------------------------
.onLoad <- function(libname, pkgname)
{
    MotifDb <<- .MotifDb(loadAllSources=TRUE, quiet=TRUE)
}
#-------------------------------------------------------------------------------

