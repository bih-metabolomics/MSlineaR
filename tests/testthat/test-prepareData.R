testthat::test_that("function checkData works", {

  dat = data.table::data.table(data.frame(ID = rep(-1,10), X = 1:10, Y = 1:10))
  dat2 = data.table::data.table(data.frame(ID = rep(1,10), X = -1:-10, Y = 1:10))
  dat3 = data.table::data.table(data.frame(ID = rep(1,10), X = 1, 1:10, Y = -1:8))

  COLNAMES = c(ID = "ID", REPLICATE = NULL, X = "X", Y = "Y")
  nCORE = 1
  LOG_TRANSFORM = TRUE
  MIN_FEATURE = 3
  testthat::expect_vector(checkData(dat, COLNAMES = COLNAMES)$ID, ptype = character(), size= 10)
  testthat::expect_vector(checkData(dat, COLNAMES = COLNAMES)$REPLICATE, ptype = character(), size= 10)
  testthat::expect_vector(checkData(dat, COLNAMES = COLNAMES)$X, ptype = double(), size= 10)
  testthat::expect_vector(checkData(dat, COLNAMES = COLNAMES)$Y, ptype = double(), size= 10)
  testthat::expect_error(checkData(dat2, COLNAMES = COLNAMES))
  testthat::expect_error(checkData(dat3, COLNAMES = COLNAMES))
  testthat::expect_error(checkData(dat, MIN_FEATURE = 2, COLNAMES = COLNAMES))
  testthat::expect_error(checkData(dat, nCORE = -1, COLNAMES = COLNAMES))
  testthat::expect_error(checkData(dat, LOG_TRANSFORM = "TRUE", COLNAMES = COLNAMES))

  })


testthat::test_that("function prepareData works", {

  dat1 = data.table::data.table(data.frame(IDintern = paste0("s",1 : 10), ID = rep(-1,10), rep = 1, X = 1:10, Y = 1:10))
  dat2 = data.table::data.table(data.frame(IDintern = paste0("s",1 : 10), ID = rep(-1,10), rep = 1, X = -1:8, Y = 1:10))

  testthat::expect_equal(unique(prepareData(dat1)$pch), 19)
  testthat::expect_equal(unique(prepareData(dat1)$color), "black")
  testthat::expect_equal(prepareData(dat1)$XLog[2], log(2))
  testthat::expect_equal(prepareData(dat1)$YLog[10], log(10))
  testthat::expect_warning(prepareData(dat2))
})


