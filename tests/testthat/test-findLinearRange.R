testthat::test_that("function findLinearRange works", {

  y. = "logY"
  x. = "logX"

  dat1 <- data.table::data.table(data.frame(groupIndices = rep(1,10),
                                            IDintern = paste0("s", 1:10),
                                            ID = rep(1,10),
                                            REPLICATE = rep(1,10),
                                            X = 1:10,
                                            logX = log(c(1:10)*3),
                                            Y = c(2.1, 2.7,3:10) + sample(x = seq(from = 0.1, to = 1, by = 0.1) * 0.4,size = 10,replace = F),
                                            logY = log(c(NA, NA,3:10) + sample(x = seq(from = 0.1, to = 1, by = 0.1) ,size = 10,replace = F)),
                                            DilutionPoint = 1:10,
                                            color = c("grey", "grey" , rep("black",8))))

  datOut1 <- findLinearRange(dat1,  x., y., max_res = 3, min_feature = 5, real_x = "X", slope_tol = 0.15, delta_tol = 0.182 )

  testthat::expect_true(round(datOut1[[2]]$R2) == 1)
  testthat::expect_true(datOut1[[2]]$RangeLength == 8)
  testthat::expect_true(datOut1[[2]]$slope == 1)
  testthat::expect_true(datOut1[[2]]$spearman_rho == 1)
  testthat::expect_true(datOut1[[2]]$Slope_within_Tolerance)
  testthat::expect_true(datOut1[[2]]$Linearity_Criterion_Deviation)



  dat2 <- data.table::data.table(data.frame(groupIndices = rep(1,10),
                                            IDintern = paste0("s", 1:10),
                                            ID = rep(1,10),
                                            REPLICATE = rep(1,10),
                                            X = 1:10,
                                            logX = log(c(1:10)*3),
                                            Y = c(3.5, 3.2, 5:10, 10.3, 9.7),
                                            logY = log(c(3.5, 3.2, 5:10, 10.3, 9.7)),
                                            DilutionPoint = 1:10,
                                            color = c(rep("black",10))))

  datOut2 <- findLinearRange(dat2,  x., y., max_res = 3, min_feature = 5, real_x = "X", slope_tol = 0.15, delta_tol = 0.182 )

  testthat::expect_false(all(datOut2[[1]]$positiveSlope[c(1,2,10)]))
  testthat::expect_true(datOut2[[2]]$RangeLength == 8)
  testthat::expect_true(datOut2[[2]]$enoughPointsWithinRange)
  testthat::expect_false(datOut2[[2]]$Slope_within_Tolerance)
  testthat::expect_true(datOut2[[2]]$Linearity_Criterion_Deviation)

})
