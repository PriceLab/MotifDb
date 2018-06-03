#------------------------------------------------------------------------------------------------------------------------
plotMotifs <- function(motifs)
{
   stopifnot(length(motifs) > 0)

  if(length(motifs) == 1){
    pcm <- new("pcm", mat=motifs[[1]], name=names(motifs))
    plot(pcm)
    }

  motifStack(lapply(names(motifs), function(mName) new("pfm", motifs[[mName]], name=mName)))

} # plotMotifs
#------------------------------------------------------------------------------------------------------------------------

