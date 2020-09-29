# script 2
# script to produce fasta of each GT region identified by hmmclust and implemented in get_Regions.R
# also returns cleaned data for determined genotypes


#-----------------
#inputs
#-----------------
dir = ""


#-----------------
# setup
#-----------------
library(ape)              # distance analysis & sikmple DNA object
library(stringr)          # string handling, regex
library(Biostrings)       # Provides DNAString, DNAStringSet, etc
library(GenomicRanges)    # Provides GRanges, etc
library(GenomicFeatures)  # provides Txdb

msa = list(); gt = list(); t = list(); run = list() # essentially classes
run$outdir = paste(dir, "out-clean", sep = "")
suppressWarnings( dir.create("out-clean"))
suppressWarnings( dir.create("out-clean/fasta"))

#-------------- run
gt$file = paste(dir, "out/gt_regions.csv", sep = "")
gt$regions = read.csv(gt$file,header = T)
colnames(gt$regions) = c("region", "msa_start", "msa_end", "ref_start", "ref_end")

msa$file = paste(dir,"all_raw_best_msa_man_noSG.fasta", sep = "")
msa$seq <- read.dna(msa$file, format = "fasta", as.matrix = T)

# hardcoded maybe make arguments
run$jar = "BIC.jar"
run$hmmclust = T

## functions
# for MSA -> reference genome
source("R/get_relative_ref_pos.R")

t$iter = 1
for(i in 1:length(gt$regions$region)){
#for(i in 1:10){

  t$start = gt$regions$msa_start[i]
  t$end = gt$regions$msa_end[i]
  t$seq = msa$seq[,t$start:t$end] # alignment chunk

  t$outfile_fasta = paste(dir,"out-clean/fasta/", t$iter, ".fasta", sep = "")
  write.FASTA(t$seq,t$outfile_fasta )
  

  #------------------- hmmclustering
  # extract all data, then see how any unique homologous blocks
  t$command = paste("java -jar", run$jar, t$outfile_fasta) # we have the fine detailed sequence start/ end just run now.
  t$hmmcluster = system(t$command,intern = T)
  t$empty_lines = grep(pattern = "^$", t$hmmcluster) # which lines have empty lines
  dat.out = read.table(text = t$hmmcluster[t$empty_lines[1] + 1], sep = "\t")
  colnames(dat.out) = c("clusters", "start_pos", "end_pos", "LL", "paramaters", "AIC", "AIC_relative")
  t$out1 = cbind(gt_region = t$iter, gt$regions[i,c(2,3,4,5)],dat.out[,c(1,4,5,6,7)]) # //todo output
  
  # we visually inspect that occasionally GT regions, when previously considered as start-middle-end were GTpic ()usually v weak AIC here)
  # but when considered as a whole region, are best explained as no genotypes. (usually visually appear as weak signal / divergence)
  # so here they may flag up as 0 GT's. and will be ignored in cleaned output.
  if(dat.out$clusters == 1){
    print(t$out1)
    # ignore. so del fasta and next iteration
    file.remove(t$outfile_fasta)
    next
  }
  
  # LL of data given that number of clusters 
  dat.genotype = data.frame(genotype = 1, V1 = "1", V2 = "1")[-1,]
  # for each genotype
  t$geno = 1:(length(t$empty_lines) - 1)
  for(j in t$geno){
    if(j == max(t$geno)){
      dat.temp = read.table(text = t$hmmcluste[(t$empty_lines[j+1]+2):length(t$hmmcluste)], sep = "\t")
    }else{
      dat.temp = read.table(text = t$hmmcluste[(t$empty_lines[j+1]+2):t$empty_lines[j+2]], sep = "\t")
    }
    dat.temp = data.frame(genotype = j, dat.temp)
    
    dat.genotype = rbind(dat.genotype, dat.temp)
  }
  colnames(dat.genotype) = c("genotype", "seq_index", "seq_name")
  t$out2 = cbind(gt_region = t$iter, dat.genotype)
  
  
  
  
  t$iter = t$iter + 1
  #------------------ output data
  if(i == 1){
    write.table(t$out1, paste(run$outdir, "/gt_regions.csv",sep = ""),sep = ",", append=TRUE, col.names = T,row.names = F)
    write.table(t$out2, paste(run$outdir, "/gt_assignment.csv",sep = ""),sep = ",", append=TRUE, col.names = T, row.names = F)
  }else{
    write.table(t$out1, paste(run$outdir, "/gt_regions.csv",sep = ""),sep = ",", append=TRUE, col.names = F, row.names = F)
    write.table(t$out2, paste(run$outdir, "/gt_assignment.csv",sep = ""),sep = ",", append=TRUE, col.names = F, row.names = F)
  }
  
}




#---------------------- Gene regions
dat.gene = read.csv(paste(run$outdir, "/gt_regions.csv",sep = ""))




