if(requireNamespace("testthat", quietly = TRUE)) {
  library("robustHD", quietly = TRUE)
  testthat::test_check("robustHD")
}
