#' For a given sequence alignment, fasta format, file name. Return the start and stop positions relative to a reference sequence
#'
#' @param alignment Multiple Sequence Alignment fasta file location
#' @param ref.file Fasta file containing Reference genome
#' @param ref.pattern Unique string to identify reference sequence in alignment labels
#' @param run_width how many base pairs should each genome scan be
#' @param run_java command to run java executable, required openjdk 13 or equivalent
#' @param run_model model selection either 'AIC' or 'BIC. Default is AIC
#' @param run_cores how many cores for parallel processes to use
#' @return vector representing start and stop positions relative to reference
#' @export
#'
get_regions_parallel = function(alignment = system.file("testdata/msa_cmv.fasta",package = "hmmcluster"), ref.file = system.file("testdata/ref_cmv.fasta",package = "hmmcluster"), ref.pattern = "Merlin", run_width = 200, run_java = "java", run_model = "AIC", run_cores = 1){
  #-----------------
  # setup
  #-----------------
  run = list();  msa = list();  t = list();  gtreg = list()
  options(scipen=999) # disable scientific number format
  msa$file = alignment
  run$ref = ref.file
  run$ref.pattern = ref.pattern
  run$width = run_width
  run$flank = 0
  run$twidth = run$width + 2*run$flank
  run$java = run_java
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
  
  cat(paste(paste(c("clusters", "start_pos", "end_pos", "LL", "paramaters", "AIC", "AIC_relative"),collapse = "\t"), "\n", sep = ""))
  
  # parallelisation 1
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
  command = paste0('find out/1-raw/ -name "*.fasta" | parallel --halt now,fail=1 --j ',run_cores, ' "', run$java,' -jar ' , run$jar , ' {} > {}.out" ')
  system(command)
  
  
  # parallelisation 2
  # identify regions to run trimmed
  print(paste("step 3 -" , date() ,"- identifying coarse regions of genomic structure & generating fasta files"))
  print(paste("region", "msa-start", "msa-end"))
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
      ape::write.FASTA(t$seq, file = t$outfile_fasta)
      print(paste(gtreg$number, gtreg$start, gtreg$end))
      gtreg$number = gtreg$number + 1
    }
    t$previous = t$current
  }
  
  print(paste("step 4 -" , date() ,"- running hmmcluster on each coarse structural region"))
  command = paste0('find out/2-gt/ -name "*.fasta" | parallel --halt now,fail=1 --j ',run_cores, ' "', run$java,' -jar ' , run$jar , ' {} -trim > {}.out" ')
  system(command)
  
  
  # get final regions
  print(paste("step 5 -" , date() ,"- generating table and alignemnts of refined structural regions"))
  gtreg$prev_end = 1
  gtreg$number = 1
  files = gtools::mixedsort(list.files("out/2-gt", pattern = "*.out"))
  cat(paste("region", "start", "end", "clusters", "ref_start", "ref_end", "LLikelihood", "Parameters", "AIC", "AIC_relative", sep = ","),sep = "\n")
  cat(paste("region", "start", "end", "clusters", "ref_start", "ref_end", "LLikelihood", "Parameters", "AIC", "AIC_relative", sep = ","),sep = "\n",file = "out/3-clean/gt-regions.csv",append = F)
  cat(paste("region", "genotype", "index", "sample", sep = ","),sep = "\n",file = "out/3-clean/gt-assignment.csv",append = F)
  for(i in 1:length(files)){
    infile = files[i]
    infile = paste0("out/2-gt/", infile)
    t$text = readLines(infile)
    t$empty_lines = grep(pattern = "^$", t$text) # which lines have empty lines
    dat.out = utils::read.table(text = t$text[t$empty_lines[1] + 1], sep = "\t") # key output line
    gtreg$clusters = dat.out[1, 1] # value = number of clusters
    gtreg$start_orig = as.numeric(stringr::str_extract(infile, "[0-9]{3,9}"))
    gtreg$start = as.numeric(gtreg$start_orig + dat.out[2])
    gtreg$end = as.numeric(gtreg$start_orig + dat.out[3])
    gtreg$LogLikelihood = as.numeric(dat.out[4])
    gtreg$params = as.numeric(dat.out[5])
    gtreg$AIC = as.numeric(dat.out[6])
    gtreg$AIC_relative = as.numeric(dat.out[7])
    if(gtreg$clusters == 1){
      # ignore it
    }else{
      # is more than one cluster
      
      ### write assignment table
      dat.genotype = data.frame(genotype = 1, V1 = "1", V2 = "1")[-1,]
      # for each genotype
      t$geno = 1:(length(t$empty_lines) - 1)
      for(gt in t$geno){
        if(gt == max(t$geno)){
          dat.temp = utils::read.table(text = t$text[(t$empty_lines[gt+1]+2):length(t$text)], sep = "\t")
        }else{
          dat.temp = utils::read.table(text = t$text[(t$empty_lines[gt+1]+2):t$empty_lines[gt+2]], sep = "\t")
        }
        dat.temp = data.frame(genotype = gt, dat.temp)
        
        dat.genotype = rbind(dat.genotype, dat.temp)
      }
      dat.genotype = data.frame(region = gtreg$number, dat.genotype)
      utils::write.table(dat.genotype,file = "out/3-clean/gt-assignment.csv", sep = "," , append = T,col.names = F,row.names = F)
      
      
      
      ### write fasta
      t$seq = msa$seq[, gtreg$start:gtreg$end] # alignment chunk
      t$outfile_fasta = paste("out/3-clean/", gtreg$number,"-", gtreg$start, "-", gtreg$end , ".fasta", sep = "")
      ape::write.FASTA(t$seq, file = t$outfile_fasta)
      
      # get reference position
      gtreg$ref_ss = get_relative_ref_pos_safe(t$seq, run$ref, run$ref.pattern)
      gtreg$ref_start = gtreg$ref_ss[1]
      gtreg$ref_end = gtreg$ref_ss[2]
      
      # write data
      cat(paste(gtreg$number, gtreg$start, gtreg$end, gtreg$clusters, gtreg$ref_start, gtreg$ref_end, gtreg$LogLikelihood, gtreg$params, gtreg$AIC, gtreg$AIC_relative, sep = ","),sep = "\n")
      cat(paste(gtreg$number, gtreg$start, gtreg$end, gtreg$clusters, gtreg$ref_start, gtreg$ref_end, gtreg$LogLikelihood, gtreg$params, gtreg$AIC, gtreg$AIC_relative, sep = ","),sep = "\n",file = "out/3-clean/gt-regions.csv",append = T)
      gtreg$number = gtreg$number + 1
      gtreg$prev_end = gtreg$end
    }
  }
  #return a dataframe
  out = utils::read.csv("out/3-clean/gt-regions.csv")
  
  print(paste("step X -" , date() ,"- hmmcluster has finished"))
  
  #
  # optional step to refine between all regions?
  #
  return(out)
}
