testthat::test_that("function trimEnds works", {

  y. = "Y"
  x. = "X"

  dat1 <- data.table::data.table(data.frame(groupIndices = rep(1,10),
                                            IDintern = paste0("s", 1:10),
                                            ID = rep(1,10),
                                            REPLICATE = rep(1,10),
                                            X = 1:10,
                                            Y = c(10, 3,2:9),
                                            DilutionPoint = 1:10,
                                            color = rep("black",10)))

  datOut1 <- chooseModel(dat1, y., x.,  abbr = "FOD")$`1`$dat

  testthat::expect_true(all(trimEnds(dats = datOut1,y = "Y_FOD", x., thresh = 0)$trim[2:4]))
  testthat::expect_true(all(trimEnds(datOut1,y = "Y_FOD", x., thresh = 0.5)$trim[2:3]))
  testthat::expect_equal(trimEnds(datOut1,y = "Y_FOD", x.)$trim[1], NA)
  testthat::expect_match(trimEnds(datOut1,y = "Y_FOD", x.)$Comment[2], "trim: firstPoint")
  testthat::expect_match(trimEnds(datOut1,y = "Y_FOD", x.)$Comment[3], "trim: <firstPoint")


  dat2 <- data.table::data.table(data.frame(groupIndices = rep(1,10),
                                            IDintern = paste0("s", 1:10),
                                            ID = rep(1,10),
                                            REPLICATE = rep(1,10),
                                            X = 1:10,
                                            Y = c(10, 3,2:5,9,6,8,7),
                                            DilutionPoint = 1:10,
                                            color = rep("black",10)))

  datOut2 <- chooseModel(dat2, y., x., abbr = "FOD")$`1`$dat

  testthat::expect_true(all(trimEnds(datOut2,y = "Y_FOD", x)$trim[2:4]))
  testthat::expect_equal(trimEnds(datOut2,y = "Y_FOD", x.)$trim[1], NA)
  testthat::expect_match(trimEnds(datOut2,y = "Y_FOD", x.)$Comment[2], "trim: firstPoint")
  testthat::expect_match(trimEnds(datOut2,y = "Y_FOD", x.)$Comment[3], "trim: <firstPoint")
  testthat::expect_match(trimEnds(datOut2,y = "Y_FOD", x.)$Comment[9], "trim: >lastPoint")
  testthat::expect_match(trimEnds(datOut2,y = "Y_FOD", x.)$Comment[10], "trim: lastPoint")

})
