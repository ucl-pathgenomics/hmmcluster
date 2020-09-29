
# hmmcluster
A method to identify homologous sequence clusters from a multiple sequence alignment of nucleotide sequences in FASTA format.
Standalone java programs and implementation to find such structural blocks in a MSA.

## method
Treats each sequence as it's own HMM, then finds the optimal way to cluster sequences using the Log Likelihood and Parameters of such a model of HMM's to determine the optimal such case.


## Repository structure

 - AIC/BIC.jar - A set of Hidden Markov Model sequencing clustering models. 
	 - BIC is the most penalising model
	 - AIC param2 is the middle model.
	 - AIC param1 is the least penalising model.
 - get_regions.R - Runs the model coarsely scanning along an MSA to find regions of interest.
 - get_regions_refined - Refines the determined regions, to provide the best start & stop positions.

## Examples
#### Running the models from a terminal
java -jar BIC.jar my_alignment.fasta
#### Running the models from a terminal with trimming. (to determine optimal start / stop positions via Maximum Likelihood)
java -jar BIC.jar my_alignment.fasta -trim


#### Automate finding regions using R
Rscript get_regions.R my_alingment.fasta my_reference_genome.fasta

Rscript get_regions_refined.R my_alingment.fasta  my_reference_genome.fasta

< outputs in ./out-clean >