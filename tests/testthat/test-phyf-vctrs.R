# Test for pf_as_pfc.phylo
test_that("pf_as_pfc.phylo returns a pfc object", {
  tree <- ape::rtree(10)
  pfc_obj <- pf_as_pfc(tree)
  expect_s3_class(pfc_obj, "pfc")
})

# Test for pf_as_pfc.Matrix
test_that("pf_as_pfc.Matrix returns a pfc object", {
  mat <- Matrix::Matrix(0, nrow = 5, ncol = 5)
  mat[lower.tri(mat)] <- runif(10)
  pfc_obj <- pf_as_pfc(mat, is_tip = 1:5)
  expect_s3_class(pfc_obj, "pfc")
})

# Test for pf_as_pfc.matrix
test_that("pf_as_pfc.matrix returns a pfc object", {
  mat <- matrix(0, nrow = 5, ncol = 5)
  mat[lower.tri(mat)] <- runif(10)
  pfc_obj <- pf_as_pfc(mat, is_tip = 1:5)
  expect_s3_class(pfc_obj, "pfc")
})

# Test for new_pfc
test_that("new_pfc returns a pfc object", {
  pfc_obj <- new_pfc()
  expect_s3_class(pfc_obj, "pfc")
})

# Tests for pfp
test_that("pfp creates a new pfp object", {
  expect_s3_class(pfp(), "pfp")
})

test_that("is_pfp detects if an object is a pfp object", {
  x <- pfp(ape::nodepath(ape::rtree(100)))
  expect_true(pf_is_pfp(x))
})

test_that("format.pfp returns the expected output with default arguments", {
  path <- pfp(list(c(1, 2, 3)))
  expect_equal(format(path), "1->2->3")
})

test_that("format.pfp returns the expected output with the 'ansi' argument set to TRUE", {
  path <- pfp(list(c(1, 2, 3)))
  expect_equal(format(path, ansi = TRUE), "1->2->3")
})

test_that("format.pfp returns the expected output with the 'colour' argument set to TRUE", {
  path <- pfp(list(c(1, 2, 3)))
  expect_equal(format(path, colour = TRUE), "1->2->3")
})

test_that("format.pfp returns the expected output with the 'ansi' and 'colour' arguments set to TRUE", {
  path <- pfp(list(c(1, 2, 3)))
  expect_equal(format(path, ansi = TRUE, colour = TRUE), "1->2->3")
})

test_that("pfp returns a list", {
  expect_type(pfp(), "list")
})

# Tests for pfc
test_that("pfc creates a new pfc object", {
  expect_s3_class(pfc(), "pfc")
})

test_that("pf_is_pfc detects if an object is a pfc object", {
  x <- pfc()
  expect_true(pf_is_pfc(x))
})

test_that("format() formats the pf, pfc, and pfp objects correctly", {
  # Create a test pfc object
  tree <- ape::read.tree(text = '((A:0.1,B:0.2):0.3,(C:0.4,(D:0.5,E:0.6):0.7):0.8);')
  pfc_test <- pf_as_pfc(tree)
  pf_test <- pf_as_pf(tree)
  pfp_test <- pfp(ape::nodepath(tree))
  
  # Check that the output of format.pfc() matches the snapshot
  expect_snapshot(format(pfc_test))
  expect_snapshot(format(pfp_test))
  expect_s3_class(pfp_test, "pfp")
  expect_s3_class(pfc_test, "pfc")
  expect_s3_class(pf_test, "pf")
  expect_s3_class(pf_test, "tbl_df")
  #expect_snapshot(format(pf_test))
})


