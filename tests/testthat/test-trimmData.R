testthat::test_that("function trimEnds works", {

  y = "Y"
  x = "X"

  dat1 <- data.table::data.table(data.frame(groupIndices = rep(1,10),
                                            IDintern = paste0("s", 1:10),
                                            ID = rep(1,10),
                                            REPLICATE = rep(1,10),
                                            X = 1:10,
                                            Y = c(10, 3,2:9),
                                            DilutionPoint = 1:10))

  datOut1 <- chooseModel(dat5, tidyselect::all_of(y), tidyselect::all_of(x), model)$`1`$dat
  datOut1$OutlierFOD <- datOut1$outlier
  datOut1 <- datOut1[,-"outlier"]


  testthat::expect_true(all(trimEnds(datOut1,y, x)$trim[2:4]))
  testthat::expect_equal(trimEnds(datOut1,y, x)$trim[1], NA)
  testthat::expect_match(trimEnds(datOut1,y, x)$Comment[2], "trim: firstPoint")
  testthat::expect_match(trimEnds(datOut1,y, x)$Comment[3], "trim: <firstPoint")


  dat2 <- data.table::data.table(data.frame(groupIndices = rep(1,10),
                                            IDintern = paste0("s", 1:10),
                                            ID = rep(1,10),
                                            REPLICATE = rep(1,10),
                                            X = 1:10,
                                            Y = c(10, 3,2:5,9,6,8,7),
                                            DilutionPoint = 1:10))

  datOut2 <- chooseModel(dat2, tidyselect::all_of(y), tidyselect::all_of(x), model)$`1`$dat
  datOut2$OutlierFOD <- datOut2$outlier
  datOut2 <- datOut2[,-"outlier"]

  testthat::expect_true(all(trimEnds(datOut2,y, x)$trim[2:4]))
  testthat::expect_equal(trimEnds(datOut2,y, x)$trim[1], NA)
  testthat::expect_match(trimEnds(datOut2,y, x)$Comment[2], "trim: firstPoint")
  testthat::expect_match(trimEnds(datOut2,y, x)$Comment[3], "trim: <firstPoint")
  testthat::expect_match(trimEnds(datOut2,y, x)$Comment[9], "trim: >lastPoint")
  testthat::expect_match(trimEnds(datOut2,y, x)$Comment[10], "trim: lastPoint")

})
