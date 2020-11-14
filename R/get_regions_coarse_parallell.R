#' For a given sequence alignment, fasta format, file name. Return the start and stop positions relative to a reference sequence. Quick and Dirty method only.
#'
#' @param alignment Multiple Sequence Alignment fasta file location
#' @param ref.file Fasta file containing Reference genome
#' @param ref.pattern Unique string to identift reference sequence in alignment labels
#' @param run_width how many base pairs should each genome scan be
#' @param run_java command to run java executable, equired openjdk 13 or equivilent
#' @param run_model model selection either 'AIC' or 'BIC. Defaut is AIC
#' @return vector representing start and stop positions relative to reference
#' @export
#'
get_regions_coarse_parallel = function(alignment = system.file("extdata/cmv_msa.fasta",package = "hmmcluster"), ref.file = system.file("ref/NC_006273.2.fasta",package = "hmmcluster"), ref.pattern = "Merlin", run_width = 200, run_java = "/opt/jdk-13/bin/java", run_model = "AIC"){
  warning("this is unstable and unfinished!")
  stop()
  #-----------------
  # setup
  #-----------------
  # library(Biostrings)       # Provides DNAString, DNAStringSet, etc
  # library(GenomicRanges)    # Provides GRanges, etc
  # library(GenomicFeatures)  # provides Txdb
  run = list();  msa = list();  t = list();  gtreg = list()
  options(scipen=999) # disable scientific number format
  msa$file = alignment
  run$ref = ref.file
  run$ref.pattern = ref.pattern
  run$width = run_width
  run$flank = 0
  run$twidth = run$width + 2*run$flank
  run$java = run_java # "/opt/jdk-13/bin/java"
  if(run_model == "AIC"){
    run$jar = system.file("AIC.jar",package = "hmmcluster")
  }else if(run_model == "BIC"){
    run$jar = system.file("BIC.jar",package = "hmmcluster")
  }else{
    warning("obs: Please check your model selection is 'AIC' or 'BIC'!")
    stop()
  }
  
  # check requirements
  if(grepl(Sys.info()[1], "Linux")){
    if(nchar(system("which parallel")) >= 1){}else{
      warning("obs: parallelisation requires you have gnu parallell installed, please install it!")
      stop()
    }
  }else{
    warning("obs: This script must be run on linux!")
    stop()
  }
  
  # create dirs
  if(!dir.exists("out")){dir.create("out")}
  if(!dir.exists("out/1-raw")){dir.create("out/1-raw")}
  if(!dir.exists("out/2-gt")){dir.create("out/2-gt")}
  if(!dir.exists("out/3-clean")){dir.create("out/3-clean")}
  # clean dirs 
  suppressMessages(do.call(file.remove, list(list.files("out/1-raw", full.names = TRUE))))
  suppressMessages(do.call(file.remove, list(list.files("out/2-gt", full.names = TRUE))))
  suppressMessages(do.call(file.remove, list(list.files("out/3-clean", full.names = TRUE))))
  
  
  
  #-----------------
  # run
  #-----------------
  msa$seq <- ape::read.dna(msa$file, format = "fasta", as.matrix = T)
  msa$len = length(msa$seq[1,])
  msa$genomes = length(labels(msa$seq))
  run$iter = seq(run$flank, msa$len - (run$width + run$flank+1), run$width)
  
  t$tbases = run$twidth * msa$genomes # total bases
  t$previous = FALSE; gtreg$prev_end = 1
  t$loops = 0 ; gtreg$loops = 1
  
  #cat(paste(paste(c("clusters", "start_pos", "end_pos", "LL", "paramaters", "AIC", "AIC_relative"),collapse = "\t"), "\n", sep = ""))
  
  # parallellisation 1
  # process each fasta into .out file for reading in later
  print(paste("step 1 -" , date() ,"- generating a set of genome chunk fasta files"))
  for(i in 1:(length(run$iter) - 1)){
    t$start = (run$iter[i] - run$flank)
    t$end = (run$iter[i+1] + run$flank - 1)
    t$seq = msa$seq[,t$start:t$end] # alignment chunk
    #---------- ignore if it is mostly garbage!
    t$num_indel = stringr::str_count(paste(as.character(t$seq),collapse = ""), "-")
    t$frac_indel = t$num_indel / t$tbases
    if(t$frac_indel > 0.4){next}
    #---------- write fasta
    t$outfile_fasta = paste("out/1-raw/",run$iter[i],".fasta", sep = "")
    ape::write.FASTA(t$seq, file = t$outfile_fasta)
  }
  
  print(paste("step 2 -" , date() ,"- running hmmcluster on each genome chunk"))
  system(paste0('find out/1-raw/ -name "*.fasta" | parallel " ', run$java,' -jar ' , run$jar , ' {} > {}.out" '))
  
  
  # parallellisation 2
  # identify regions to run trimmed
  print(paste("step 3 -" , date() ,"- identifying coarse regions of genomic structure & generating fasta files"))
  print(paste("region", "start", "end", "ref_start", "ref_end"))
  t$previous = 1
  gtreg$prev_end = 1
  gtreg$number = 1
  files = gtools::mixedsort(list.files("out/1-raw", pattern = "*.out"))
  for(i in 1:length(files)){
    infile = files[i]
    infile = paste0("out/1-raw/", infile)
    t$text = readLines(infile)
    t$empty_lines = grep(pattern = "^$", t$text) # which lines have empty lines
    # key output line
    dat.out = utils::read.table(text = t$text[t$empty_lines[1] + 1], sep = "\t")
    t$current = dat.out[1, 1] # value = number of clusters
    if (t$previous == 1 && t$current == 1 ) {
      # do nothing
    }
    if (t$previous != 1 && t$current != 1 ) {
      #do nothing
    }
    if (t$previous == 1 && t$current != 1 ) {
      # if start of structural region
      gtreg$start = as.numeric(stringr::str_extract(infile, "[0-9]{3,9}"))
      gtreg$start = gtreg$start - (run$width / 2)
    }
    if( t$current == 1 & t$previous != 1 ){
      gtreg$end = as.numeric(stringr::str_extract(infile, "[0-9]{3,9}"))
      gtreg$end = gtreg$end + run$width / 2
      # create region fasta
      t$seq = msa$seq[, gtreg$start:gtreg$end] # alignment chunk
      
      # get reference position
      gtreg$ref_ss = get_relative_ref_pos_safe(t$seq, run$ref, run$ref.pattern)
      gtreg$ref_start = gtreg$ref_ss[1]
      gtreg$ref_end = gtreg$ref_ss[2]
      
      print(paste(gtreg$number, gtreg$start, gtreg$end, gtreg$ref_start, gtreg$ref_end, sep = ","))
      gtreg$number = gtreg$number + 1
    }
    t$previous = t$current
  }
  
  
  
}