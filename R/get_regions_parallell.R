# script to run hmmcluster, given a reference genome and MSA.
# automatically spit out the start and stop of genotypic regions.
# as a csv file.
# This runs in parallell on linux mchines for maximum SPEED!!

#' For a given sequence alignment, fasta format, file name. Return the start and stop positions relative to my reference sequence
#'
#' @param alignment Multiple Sequence Alignment fasta file location
#' @param ref.file Fasta file containing Reference genome
#' @param ref.pattern 
#' @return vector representing start and stop positions relative to reference
#' @export
#'
get_regions_parallel = function(alignment = "all_raw_best_msa_man3b.fasta", ref.file = "ref/NC_006273.2.fasta", ref.pattern = "Merlin"){
  #-----------------
  # setup
  #-----------------
  # library(Biostrings)       # Provides DNAString, DNAStringSet, etc
  # library(GenomicRanges)    # Provides GRanges, etc
  # library(GenomicFeatures)  # provides Txdb
  run = list()
  msa = list()
  t = list() # temp for each loop
  gtreg = list() # temp gt data
  options(scipen=999) # disable scientific number format
  
  # handle command line
  if(length(commandArgs(trailingOnly = TRUE)) > 0){
    args <- commandArgs(trailingOnly = TRUE)
    msa$file = args[1]
    run$ref = args[2]
  }else{ # run default for testing
    msa$file = alignment
    run$ref = ref.file
    run$ref.pattern = ref.pattern
  }
  
  #----------------------------- INSTRUCTIONS
  run$width = 200 # 300 # //todo make smaller
  run$flank = 0 # 100
  run$twidth = run$width + 2*run$flank
  run$java = "/opt/jdk-13/bin/java"
  run$jar = "AIC.jar"  #model selection
  run$hmmclust = T
  #----------------------------- end INSTRUCTIONS
  
  # dna characters regex
  pattern = "[^ACTG-]"
  
  # create dir
  if(!dir.exists("out")){dir.create("out")}
  if(!dir.exists("out/1-raw")){dir.create("out/1-raw")}
  if(!dir.exists("out/2-gt")){dir.create("out/2-gt")}
  if(!dir.exists("out/3-clean")){dir.create("out/3-clean")}
  
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
  
  cat(paste(paste(c("clusters", "start_pos", "end_pos", "LL", "paramaters", "AIC", "AIC_relative"),collapse = "\t"), "\n", sep = ""))
  
  # parallellisation 1
  # process each fasta into .out file for reading in later
  for(i in 1:length(run$iter - 1)){
  
    t$start = (run$iter[i] - run$flank)
    t$end = (run$iter[i+1] + run$flank - 1)
    t$seq = msa$seq[,t$start:t$end] # alignment chunk
    #---------- ignore if it is mostly garbage!
    t$num_indel = stringr::str_count(paste(as.character(t$seq),collapse = ""), "-")
    t$frac_indel = t$num_indel/ t$tbases
    if(t$frac_indel > 0.4){next}
    if(run$hmmclust == T){
      #---------- write fasta
      t$outfile_fasta = paste("out/1-raw/",run$iter[i],".fasta", sep = "")
      write.FASTA(t$seq, file = t$outfile_fasta)
    }
  }
  date()
  system(paste0('find out/1-raw/ -name "*.fasta" | parallel "', run$java,' -jar AIC.jar {} > {}.out"'))
  date()
  
  
  # parallellisation 2
  # identify regions to run trimmed
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
      t$outfile_fasta = paste("out/2-gt/", gtreg$start, ".fasta", sep = "")
      write.FASTA(t$seq, file = t$outfile_fasta)
      gtreg$number = gtreg$number + 1
      print(paste(gtreg$number, gtreg$start, gtreg$end))
    }
    t$previous = t$current
  }
  date()
  system(paste0('find out/2-gt/ -name "*.fasta" | parallel "', run$java,' -jar AIC.jar {} -trim > {}.out"'))
  date()
  
  
  
  # get final regions
  gtreg$prev_end = 1
  gtreg$number = 1
  files = gtools::mixedsort(list.files("out/2-gt", pattern = "*.out"))
  cat(paste("region", "start_orig", "start", "end", "clusters", "note", "ref_start", "ref_end", sep = ","),sep = "\n",file = "out/3-clean/gt-regions.csv",append = F)
  for(i in 1:length(files)){
    infile = files[i]
    infile = paste0("out/2-gt/", infile)
    t$text = readLines(infile)
    t$empty_lines = grep(pattern = "^$", t$text) # which lines have empty lines
    dat.out = utils::read.table(text = t$text[t$empty_lines[1] + 1], sep = "\t") # key output line
    gtreg$clusters = dat.out[1, 1] # value = number of clusters
    gtreg$start_orig = as.numeric(stringr::str_extract(infile, "[0-9]{3,9}"))
    if(gtreg$clusters == 1){
      # ignore it
    }else{
      # is more than one cluster
      gtreg$start = gtreg$start_orig + dat.out[1,2]
      gtreg$end = gtreg$start_orig + dat.out[1,3]
      if(gtreg$start <= gtreg$prev_end){
        gtreg$note = "YES"}
      else{
        gtreg$note = "NO"
      }
      ### write fasta
      t$seq = msa$seq[, gtreg$start:gtreg$end] # alignment chunk
      t$outfile_fasta = paste("out/3-clean/", gtreg$number,"-", gtreg$start, "-", gtreg$end , ".fasta", sep = "")
      write.FASTA(t$seq, file = t$outfile_fasta)
      
      # get reference position
      gtreg$ref_start = get_relative_ref_pos_safe(t$seq, run$ref)[1]
      gtreg$ref_start = get_relative_ref_pos(t$seq, run$ref, run$ref.pattern)[1]
      gtreg$ref_end = get_relative_ref_pos(t$seq, run$ref, run$ref.pattern)[2]
      
      # write data
      cat(paste(gtreg$number, gtreg$start_orig, gtreg$start, gtreg$end, gtreg$clusters, gtreg$note, gtreg$ref_start, gtreg$ref_end, sep = ","),sep = "\n")
      cat(paste(gtreg$number, gtreg$start_orig, gtreg$start, gtreg$end, gtreg$clusters, gtreg$note, gtreg$ref_start, gtreg$ref_end, sep = ","),sep = "\n",file = "out/3-clean/gt-regions.csv",append = T)
      gtreg$number = gtreg$number + 1
      gtreg$prev_end = gtreg$end
    }
  }
}