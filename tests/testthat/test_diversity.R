# this is a set of unit tests for the alignment -> diversity

alignment = system.file("extdata/cmv_msa.fasta", package = "hmmcluster")



test_that("diversity is functional", {
  # relative_reference_position
  dat = hmmcluster::alignment_heterozygosity(alignment)
  
  expect_equal(mean(dat$hzyg),0.02454773)
  
  
})
