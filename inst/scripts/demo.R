library (MotifDb)
library (MotIV)
library (seqLogo)
library (RUnit)
#-----------------------------------------------------------------------------------------
demo.MotIV = function ()
{
  print ('--- demo.MotIV')
  mdb = MotifDb (quiet=TRUE)
  db.tmp = mdb@listData

     # match motif 1 against the entire MotifDb  collection
  motif.hits = motifMatch (db.tmp [1], database=db.tmp)

     # the long way to extract the matrix name.  see MotIV.toTable below for more convenient way
  checkEquals (motif.hits@bestMatch[[1]]@aligns[[1]]@TF@name, names (db.tmp)[1])

     # jaspar, uniprobe and ScerTF contribute a total of 1035 matrices
  checkTrue (length (mdb) >= 1035)

  target = db.tmp [1035]
  printf ('target matrix: %s', names (target))
  motif.hits =  motifMatch (db.tmp [1035], database=db.tmp)
  tbl.hits = MotIV.toTable (motif.hits)
  checkEquals (names (db.tmp)[1035], tbl.hits [1, 'name'])   # first hit should be the target itself

  return (tbl.hits)

} # demo.MotIV
#------------------------------------------------------------------------------------------------------------------------
display.matrices= function ()
{
  mdb = MotifDb (quiet=TRUE)
  printf ('--- %s', 'ScerTF-Scerevisiae-PDR8-badis')
  print (subset (mdb, name=='ScerTF-Scerevisiae-PDR8-badis')[[1]])
  printf ('--- %s', 'JASPAR_CORE-Scerevisiae-YDR520C-MA0422.1')
  print (subset (mdb, name=='JASPAR_CORE-Scerevisiae-YDR520C-MA0422.1')[[1]])
  printf ('--- %s', 'ScerTF-Scerevisiae-PDR8-badis')
  print (subset (mdb, name=='ScerTF-Scerevisiae-PDR8-badis')[[1]])
  printf ('--- %s', 'ScerTF-Scerevisiae-YLR278C-badis')
  print (subset (mdb, name=='ScerTF-Scerevisiae-YLR278C-badis')[[1]])
  printf ('--- %s', 'JASPAR_CORE-Scerevisiae-HAP1-MA0312.1')
  print (subset (mdb, name=='JASPAR_CORE-Scerevisiae-HAP1-MA0312.1')[[1]])
  printf ('--- %s', 'ScerTF-Scerevisiae-YKL222C-zhu')
  print (subset (mdb, name=='ScerTF-Scerevisiae-YKL222C-zhu')[[1]])

} # display.logos
#------------------------------------------------------------------------------------------------------------------------
display.logos = function ()
{
  mdb = MotifDb (quiet=TRUE)
  par (mfrow=c (3,2))
  seqLogo (subset (mdb, name=='ScerTF-Scerevisiae-PDR8-badis')[[1]])
  #readline ()
  seqLogo (subset (mdb, name=='JASPAR_CORE-Scerevisiae-YDR520C-MA0422.1')[[1]])
  seqLogo (subset (mdb, name=='ScerTF-Scerevisiae-PDR8-badis')[[1]])
  seqLogo (subset (mdb, name=='ScerTF-Scerevisiae-YLR278C-badis')[[1]])
  seqLogo (subset (mdb, name=='JASPAR_CORE-Scerevisiae-HAP1-MA0312.1')[[1]])
  seqLogo (subset (mdb, name=='ScerTF-Scerevisiae-YKL222C-zhu')[[1]])

} # display.logos
#------------------------------------------------------------------------------------------------------------------------
MotIV.toTable = function (match)
{
  stopifnot (length (match@bestMatch) == 1)
  alignments = match@bestMatch[[1]]@aligns

  df = data.frame (stringsAsFactors=FALSE)
  for (alignment in alignments) {
    x = alignment
    name = x@TF@name
    eVal = x@evalue
    sequence = x@sequence
    match = x@match
    strand = x@strand
    df = rbind (df, data.frame (name=name, eVal=eVal, sequence=sequence, match=match, strand=strand, stringsAsFactors=FALSE))
    } # for alignment

  return (df)

} # MotIV.toTable 
#------------------------------------------------------------------------------------------------------------------------
test.MotIV.toTable = function ()
{
  print ('--- test.MotIVtoTable')
  mdb = MotifDb (quiet=TRUE)
  test.hits = motifMatch (mdb[1]@listData, database=jaspar)
  tbl.hits =  MotIV.toTable (test.hits)
  checkEquals (dim (tbl.hits), c (5, 5))
  checkEquals (colnames (tbl.hits), c ("name", "eVal", "sequence", "match", "strand"))

} # test.MotIV.toTable 
#------------------------------------------------------------------------------------------------------------------------
demo.tomtom = function ()
{
  mdb = MotifDb (quiet=TRUE)
  mdb.01 = mdb [1035]
  max = length (mdb)

  mdb.many = mdb [1:max]
  export (mdb.01, 'mdb1.text', 'meme')
  export (mdb.many, 'mdbMany.text', 'meme')
     # find similarity of motif #1 to all the motifs in mdbMany
  cmd = sprintf ('tomtom -no-ssc -oc ../tomtomTest -verbosity 3 -min-overlap 5 -mi 1 -dist pearson -evalue -thresh 10 %s %s',
                  'mdb1.text',  'mdbMany.text')
  system (cmd)
  cmd = 'open ../tomtomTest/tomtom.html'   # file.path, browseURL
  system (cmd)

} # demo.tomtom
#------------------------------------------------------------------------------------------------------------------------
