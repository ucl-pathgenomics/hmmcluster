# hmmcluster



## Overview
This repo contains the code for a nucleotide sequence clustering algorithm based on Hidden Markov Models. This approach effectively addresses the challenge of non-homologous segments within the alignment. 


## Getting started
The easiest way to use this program is via the packaged jar file available for download from the release page to the right here. 
Then assuming Java is installed you run the program by providing a single argument which is a multi-sequences alignment (MSA) file in fasta format e.g. `java -jar hmmcluster.jar test/msa.fa`


By default, the program considers the alignment verbatim, but with the `-trim` option the underlying statistical approach can be used to find the optimal start and end base positions of the cluster-able region. e.g. `java -jar hmmcluster.jar test/msa.fa -trim `


## Method
For full details see our paper ![Genomic and geographical structure of human cytomegalovirus](https://www.pnas.org/doi/10.1073/pnas.2221797120)

From a given sequence alignment, hmmcluster determines the optimal number of sequence clusters that best explain the diversity across the genomes. This approach groups together sequences based on the statistical likelihood that they come from the same underlying source. 

## Output
xx

# Development and contacts
Richard Goldstein UCL & ![Oscar Charles UCL](mailto:oscar.charles.18@ucl.ac.uk) & ![Cristina Venturini UCL](mailto:c.venturini@ucl.ac.uk)
