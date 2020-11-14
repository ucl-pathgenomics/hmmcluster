#' Return the start and stop positions relative to my reference sequence
#'
#' @param alignment Multiple Sequence Alignment as dnabin
#' @param ref.file Fasta file containing Reference genome
#' @param ref.pattern Unique string to identift reference sequence in alignment labels
#' @return vector representing start and stop positions relative to reference
#' @export
#'
get_relative_ref_pos = function(alignment, ref.file = system.file("ref/NC_006273.2.fasta",package = "hmmcluster"), ref.pattern = "Merlin"){
  # takes alignment and reference merlin file
  # returns start, end location relative to merlin
  
  pattern = "[^ACTG-]"
  
  ref.seq = ape::read.dna(ref.file,format = "fasta", as.matrix = T)
  ref.seq.string = paste(as.character(ref.seq),collapse = "")
  
  alignemnt.ref.num = grep(x = labels(alignment), pattern = ref.pattern)
  alignment.ref.seq = alignment[alignemnt.ref.num,]
  
  # build string of first 20 characters from reference in alignment. ignoring indels. 30 felt fine.
  t = 1; t.start20 = c()
  while(length(t.start20) < 30){
    t.string = as.character(alignment.ref.seq[t])
    if(grepl(pattern, t.string)){
      t.start20 = c(t.start20,as.character(alignment.ref.seq[t]))
    }
    t = t+1
    #print(t.start20)
  }
  t.start20 = paste(t.start20, collapse = "")
  
  
  # build string of last 20 characters from reference in alignment. ignoring indels
  t = length(alignment[1,]); t.end20 = c() # from last value work back
  while(length(t.end20) < 30){
    t.string = as.character(alignment.ref.seq[t])
    if(grepl(pattern, t.string)){
      t.end20 = c(as.character(alignment.ref.seq[t]), t.end20)
    }
    t = t-1
    #print(t.end20)
  }
  t.end20 = paste(t.end20, collapse = "")
  
  
  ## locate the location of these strings in the reference sequence and output.
  t.start.merlinpos = stringr::str_locate(pattern = t.start20, string = ref.seq.string)[1] # grab pos of first base in start string
  t.end.merlinpos = stringr::str_locate(pattern = t.end20, string = ref.seq.string)[2] # grab pos of last base in end string
  
  return(c(t.start.merlinpos, t.end.merlinpos))
}



#' Safely, escape if reference sequence has been corrupted. Return the start and stop positions relative to my reference sequence
#'
#' @param alignment Multiple Sequence Alignment fasta file location
#' @param ref.file Fasta file containing Reference genome
#' @param ref.pattern Unique string to identift reference sequence in alignment labels
#' @return vector representing start and stop positions relative to reference
#' @export
#'
get_relative_ref_pos_safe <- function(alignment, ref.file = system.file("ref/NC_006273.2.fasta",package = "hmmcluster"), ref.pattern = "Merlin") {
  time_limit <- 2 # 2 seconds
  setTimeLimit(cpu = time_limit, elapsed = time_limit, transient = TRUE)
  on.exit({
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
  })
  
  tryCatch({
    # do some stuff
    out = get_relative_ref_pos(alignment, ref.file, ref.pattern)
  }, error = function(e) {
    if (grepl("reached elapsed time limit|reached CPU time limit", e$message)) {
      # we reached timeout, apply some alternative method or do something else
      out = NA
    } else {
      # error not related to timeout
      out = NA
      stop(e)
    }
  })
}

