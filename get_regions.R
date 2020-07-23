# script to run hmmcluster, given a reference genome and MSA.
# automatically spit out the start and stop of genotypic regions.
# as a csv file.


# param1.jar is the original file
# param2.jar is richards update


#-----------------
# setup
#-----------------
library(ape)              # distance analysis & sikmple DNA object
library(stringr)          # string handling, regex
library(Biostrings)       # Provides DNAString, DNAStringSet, etc
library(GenomicRanges)    # Provides GRanges, etc
library(GenomicFeatures)  # provides Txdb
run = list()
msa = list()
t = list() # temp for each loop
gtreg = list() # temp gt data

# handle command line
if(length(commandArgs(trailingOnly = TRUE)) > 0){
  args <- commandArgs(trailingOnly = TRUE)
  msa$file = args[1]
  run$ref = args[2]
}else{ # run default for testing
  msa$file = "all_raw_best_msa_man.fasta"
  run$ref = "ref/NC_006273.2.fasta"
}

# hardcoded maybe make arguments
run$width = 200 # 300 # //todo make smaller
run$flank = 0 # 100
run$twidth = run$width + 2*run$flank
run$jar = "param1.jar"
run$hmmclust = T


# dna characters regex
pattern = "[^ACTG-]"

## functions
# for MSA -> reference genome
get_relative_ref_pos = function(alignment, ref.file = run.ref){
  # takes alignemnt and reference merlin file
  # returns start, end location relative to merlin
  
  ref.seq = read.dna(ref.file,format = "fasta", as.matrix = T)
  ref.seq.string = paste(as.character(ref.seq),collapse = "")
  
  alignemnt.ref.num = grep(x = labels(alignment), pattern = "Merlin")
  alignment.ref.seq = alignment[alignemnt.ref.num,]
  
  t = 1; t.start20 = c()
  while(length(t.start20) < 20){
    t.string = as.character(alignment.ref.seq[t])
    if(grepl(pattern, t.string)){
      t.start20 = c(t.start20,as.character(alignment.ref.seq[t]))
    }
    t = t+1
    #print(t.start20)
  }
  t.start20 = paste(t.start20, collapse = "")
  
  
  
  t = 500; t.end20 = c()
  while(length(t.end20) < 20){
    t.string = as.character(alignment.ref.seq[t])
    if(grepl(pattern, t.string)){
      t.end20 = c(as.character(alignment.ref.seq[t]), t.end20)
    }
    t = t-1
    #print(t.end20)
  }
  t.end20 = paste(t.end20, collapse = "")
  
  
  ## outputs
  t.start.merlinpos = str_locate(pattern = t.start20, string = ref.seq.string)[2]
  t.end.merlinpos = str_locate(pattern = t.end20, string = ref.seq.string)[1]
  
  return(c(t.start.merlinpos, t.end.merlinpos))
}


# create dir
if(dir.exists("out")){
}else{dir.create("out")}



#-----------------
# run
#-----------------
msa$seq <- read.dna(msa$file, format = "fasta", as.matrix = T)
msa$len = length(msa$seq[1,])
msa$genomes = length(labels(msa$seq))
run$iter = seq(run$flank, msa$len - (run$width + run$flank+1), run$width)
run$iter = run$iter[257:267] # debug

t$tbases = run$twidth * msa$genomes # total bases
t$previous = FALSE
t$loops = 0 ; gtreg$loops = 0

for(i in 1:length(run$iter - 1)){
  t$start = (run$iter[i] - run$flank)
  t$end = (run$iter[i+1] +run$flank - 1)
  t$seq = msa$seq[,t$start:t$end] # alignment chunk
  

  
  #---------- ignore if it is mostly garbage!
  t$num_indel = str_count(paste(as.character(t$seq),collapse = ""), "-")
  t$frac_indel = t$num_indel/ t$tbases
  if(t$frac_indel > 0.6){next}
  t$loops = t$loops + 1
  
  
  if(run$hmmclust == T){
    #---------- write fasta
    t$outfile_fasta = paste("out/",run$iter[i],".fasta", sep = "")
    write.FASTA(t$seq, file = t$outfile_fasta)
    
    
    #--------------- reference relative position
    t$ref_start = get_relative_ref_pos(t$seq, run$ref)[1]
    t$ref_end = get_relative_ref_pos(t$seq, run$ref)[2]
    
    
    
    
    #------------------- hmmclustering
    # extract all data, then see how any unique homologous blocks
    t$command = paste("java -jar", run$jar, t$outfile_fasta)
    t$hmmcluster = system(t$command,intern = T)
    t$outfile = paste("out/",run$iter[i],".out", sep = "")
    writeLines(t$hmmcluster, t$outfile)
    
    # best way to deal with this will be to split the data by chunks.
    t$empty_lines = grep(pattern = "^$", t$hmmcluster) # which lines have empty lines
    
    # key output line
    dat.out = read.table(text = t$hmmcluster[t$empty_lines[1] + 1], sep = "\t")
    colnames(dat.out) = c("cluters", "start_pos", "end_pos", "LL", "paramaters", "AIC", "AIC_relative")
    dat.out$start_pos[1] = t$start; dat.out$end_pos[1] = t$end 
    dat.out = cbind(dat.out, data.frame(ref_start = t$ref_start, ref_end = t$ref_end))
    # LL of data given that number of clusters 
    dat.genotype = data.frame(genotype = 1, V1 = "1", V2 = "1")[-1,]
    # for each genotype
    t$geno = 1:(length(t$empty_lines) - 1)
    for(gt in t$geno){
      if(gt == max(t$geno)){
        dat.temp = read.table(text = t$hmmcluste[(t$empty_lines[gt+1]+2):length(t$hmmcluste)], sep = "\t")
      }else{
        dat.temp = read.table(text = t$hmmcluste[(t$empty_lines[gt+1]+2):t$empty_lines[gt+2]], sep = "\t")
      }
      dat.temp = data.frame(genotype = gt, dat.temp)
      
      dat.genotype = rbind(dat.genotype, dat.temp)
    }
    
    
    
     
    
    #------------------- hmmclustering - do we trim?
    # if prev 1 now more, trim.
    # if was more than 1, now 1 trim.
    # else do nothing.
    t$num_clusters = dat.out$cluters[1]
    ifelse(t$num_clusters == 1, t$current <- FALSE, t$current <- TRUE) # if num clust > 1 this is a cluster block
    
    if(t$current == T & t$previous != t$current){ # hom -> gt: rerun this region with trim. store start pos
      print(paste(t$start, "running start hmm"))
      # need to look back at last chuck as well - cant be sure it started this chunk
      gtreg$seq = msa$seq[,(t$start - run$width):t$end] # alignment chunk
      gtreg$outfile_fasta = paste("out/temp.fasta", sep = "")
      write.FASTA(gtreg$seq, file = gtreg$outfile_fasta)
      t$command = paste("java -jar", run$jar, gtreg$outfile_fasta, "trim") # now trim
      t$hmmcluster = system(t$command,intern = T)
      cat(t$hmmcluster,file = paste("out-gt/", gtreg$loops, "start.out", sep = ""),sep = "\n",append = T)
      # best way to deal with this will be to split the data by chunks.
      t$empty_lines = grep(pattern = "^$", t$hmmcluster) # which lines have empty lines
      # key output line
      gtreg$dat.out = read.table(text = t$hmmcluster[t$empty_lines[1] + 1], sep = "\t")
      colnames(gtreg$dat.out) = c("cluters", "start_pos", "end_pos", "LL", "paramaters", "AIC", "AIC_relative")
      gtreg$start = t$start - run$width + gtreg$dat.out$start_pos[1] # store start pos
      
    }else if(t$current == F & t$previous != t$current){ # if gt -> hom: rerun this region with trim store end pos
      print(paste(t$start, "running end hmm"))
      # need to look back at last chuck as well - cant be sure it started this chunk
      gtreg$seq = msa$seq[,(t$start - run$width):t$end] # alignment chunk
      gtreg$outfile_fasta = paste("out/temp.fasta", sep = "")
      write.FASTA(gtreg$seq, file = gtreg$outfile_fasta)
      t$command = paste("java -jar", run$jar, gtreg$outfile_fasta, "trim") # now trim
      t$hmmcluster = system(t$command,intern = T)
      cat(t$hmmcluster,file = paste("out-gt/", gtreg$loops, "end.out", sep = ""),sep = "\n",append = T)
      # best way to deal with this will be to split the data by chunks.
      t$empty_lines = grep(pattern = "^$", t$hmmcluster) # which lines have empty lines
      # key output line
      gtreg$dat.out = read.table(text = t$hmmcluster[t$empty_lines[1] + 1], sep = "\t")
      colnames(gtreg$dat.out) = c("cluters", "start_pos", "end_pos", "LL", "paramaters", "AIC", "AIC_relative")
      gtreg$end = t$start + gtreg$dat.out$end_pos[1] # store end pos
      
      
      ##------ have now defined a genotypic region, and need to write it to file.
      # get ref relative positions
      gtreg$seq = msa$seq[,(gtreg$start):(gtreg$end - 1)] # alignment chunk
      gtreg$ref_start = get_relative_ref_pos(gtreg$seq, run$ref)[1]
      gtreg$ref_end = get_relative_ref_pos(gtreg$seq, run$ref)[2]
      
      
      gtreg$dat = data.frame(which = gtreg$loops,
                             msa_start = gtreg$start,
                             msa_end = gtreg$end,
                             ref_start = gtreg$ref_start,
                             ref_end = gtreg$ref_end)

      
      
      ##------- cat output
      if(gtreg$loops == 0){
        write.table(gtreg$dat, "out/gt_regions.csv",sep = ",", append=TRUE,row.names = F)
      }else{
        write.table(gtreg$dat, "out/gt_regions.csv",sep = ",", append=TRUE,col.names = F, row.names = F)
      }
      
      gtreg$loops = gtreg$loops + 1
    }else{}    # else do nothing
    
    
    
  
    # store all genotype data. for gt regions
    if(t$loops == 0){
      write.table(cbind(data.frame(iter = run$iter[i]), dat.genotype), "out/genotype_assignment.csv",sep = ",", append=TRUE)
    }else{
      write.table(cbind(data.frame(iter = run$iter[i]), dat.genotype), "out/genotype_assignment.csv",sep = ",", append=TRUE,col.names = F)
    }
  

    t$previous = t$current # update
    } # hmmcluster
    
  

  if(t$loops == 0){
    write.table(dat.out, "out/scan_res.csv",sep = ",", append=TRUE)
  }else{
    write.table(dat.out, "out/scan_res.csv",sep = ",", append=TRUE,col.names = F)
  }
  
  cat(paste(paste(dat.out,collapse = "\t"), "\n", sep = ""))
    
}

