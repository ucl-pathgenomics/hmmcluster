# this is a set of unit tests for the alignment -> regions of genomic structure processes

alignment = system.file("extdata/cmv_msa.fasta", package = "hmmcluster")
alignment = ape::read.dna(alignment,format = "fasta", as.matrix = T)
ref.file = system.file("ref/NC_006273.2.fasta",package = "hmmcluster")
ref.pattern = "Merlin"



# #------------------------ test get regions
# running on a very small subset of the msa to test method
# test_msa = alignment[,6000:7200]
# ape::write.FASTA(test_msa, "del.fasta")
# get_regions_parallel(alignment = "del.fasta",ref.file = ref.file,ref.pattern = ref.pattern,run_width = 200,run_java = "/opt/jdk-13/bin/java",run_model = "AIC")
# 
# file.remove("del.fasta")










#-------------------------------------- helper function

test_that("helper functions", {
# relative_reference_position
rrp = get_relative_ref_pos(alignment, ref.file, ref.pattern)
expect_equal(rrp[1],42122)
expect_equal(rrp[2],49391)


# relative_reference_position safe
# test with no viable kmer to reference mapping
alignment_corrupted = alignment[,c(1:10, 50:60, 500:510)]
rrps = get_relative_ref_pos(alignment_corrupted, ref.file, ref.pattern)
expect_equal(rrps[1],NA_integer_)
expect_equal(rrps[2],NA_integer_)

# test with sequence too short to build kmer, should hang and be caught and return NA
alignment_corrupted = alignment[,1:10]
rrps = get_relative_ref_pos_safe(alignment_corrupted, ref.file, ref.pattern)
expect_equal(rrps,NA)


})
