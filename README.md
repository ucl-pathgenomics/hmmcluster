# hmmcluster



## Overview
A nucleotide sequence clustering algorithm based on Hidden Markov Models that overcomes the issue where parts of alignments may not be homologous. 



## Getting started
The easiest way to use this program is via the packaged jar file available for download from the release page to the right here. 
Then assuming java is installed you run the program by providing a single argument which is a fasta file e.g. `java -jar hmmcluster.jar test/msa.fa`


By default the program considers the alignment verbatim, but with the `-trim` option the underlying statistical approach can be used to find the optimal start and end base positions of the cluster-able region. e.g. `java -jar hmmcluster.jar test/msa.fa -trim `



## Method
For full details see our paper on the ![Genomic and geographical structure of human cytomegalovirus](https://www.pnas.org/doi/10.1073/pnas.2221797120)

Coarsely put:
100 HMM's will perfectly describe 100 sequences, but 50 HMM's may do just an as good a job but with half the parameters. How many HMM's is optimal?

From a given sequence alignment where not every positions needs to be homologous, hmmcluster derives n HMM's from the n sequences present. it then attempts to find the optimal way to combine HMM's in a way that optimises the Log Likelihood of sequences given the current set of models and minimises the number of parameters. We use AIC to find the best number of HMM clusters for the given alignment.



# Development and contacts
Richard Goldstein UCL & ![Oscar Charles UCL](mailto:oscar.charles.18@ucl.ac.uk)