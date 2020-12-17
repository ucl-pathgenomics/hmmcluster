#' calculates heterozygosity for a position, in an alignment (ape package matrix)
#' @param i position
#' @param seq alignement as character matrix
#' @return vector os position and heterozygocuty
#' @export
#'
hzyg <- function(i, seq){
  # .6 secs per 1000
  # heterozygocity
  t = seq[,i]
  #n = length(t)
  t = as.data.frame(table(t))
  t = t[t$t != "n",] #remove n sequences
  t = t[t$t != "-",] #remove indel sequences, hahn book says.
  n = sum(t$Freq) # updated n only non-ambiguous bases
  t.nuc = sum((t$Freq / n)^2)
  t.hzyg = (n/(n-1))*(1-t.nuc)
  return(t.hzyg)
}

#' calculated positional heterozygosity along an alignment
#' @param alignment fasta file location
#' @return dataframe of postiion and heterozygosity
#' @export
#'

alignment_heterozygosity = function(alignment){

  mds = ape::read.dna(alignment,format = "fasta",as.matrix = T,as.character = T) #not add
  
  out = data.frame(pos = 1:ncol(mds), hzyg = 0)
  for(pos in 1:length(mds[1,])){ # for each pos
    
    #----------- heterozygocty
    t.hzyg = hzyg(pos, mds)
    out[pos,2] = t.hzyg
    
  }
  # NA values, returned from N or single sequanece insertions shuld return as 0
  out$hzyg[is.na(out$hzyg)] <- 0
  return(out)
}