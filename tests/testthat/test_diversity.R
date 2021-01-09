# this is a set of unit tests for the alignment -> diversity

alignment_file = system.file("extdata/cmv_msa.fasta", package = "hmmcluster")



test_that("diversity is functional", {
  # relative_reference_position
  dat = hmmcluster::alignment_heterozygosity(alignment_file)
  
  expect_equal(mean(dat$hzyg),0.01736416)
  
  
})
