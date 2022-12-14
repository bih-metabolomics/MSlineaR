testthat::test_that("checks function chooseModel",{
  dat1 <- data.table::data.table(data.frame(groupIndices = rep(1,10),
                                            ID = rep(1,10),
                                            REPLICATE = rep(1,10),
                                            X = 1:10,
                                            Y = 1:10,
                                            DilutionPoint = 1:10))

  dat2 <- data.table::data.table(data.frame(groupIndices = rep(1,10),
                                            ID = rep(1,10),
                                            REPLICATE = rep(1,10),
                                            X = 1:10,
                                            Y = c(1:5,4:0),
                                            DilutionPoint = 1:10))

  dat3 <- data.table::data.table(data.frame(groupIndices = rep(1,10),
                                            ID = rep(1,10),
                                            REPLICATE = rep(1,10),
                                            X = 1:10,
                                            Y = c(1,1,1,2,3,4,5,6,6,6),
                                            DilutionPoint = 1:10))

  dat4 <- data.table::data.table(data.frame(groupIndices = rep(1,10),
                                            ID = rep(1,10),
                                            REPLICATE = rep(1,10),
                                            X = 1:10,
                                            Y = c(1:3,10,5:10),
                                            DilutionPoint = 1:10))

  dat5 <- data.table::data.table(data.frame(groupIndices = rep(1,10),
                                            ID = rep(1,10),
                                            REPLICATE = rep(1,10),
                                            X = 1:10,
                                            Y = c(10, 2:10),
                                            DilutionPoint = 1:10))

  dat6 <- data.table::data.table(data.frame(groupIndices = rep(1,10),
                                            ID = rep(1,10),
                                            REPLICATE = rep(1,10),
                                            X = 1:10,
                                            Y = c(1:5,4:2,8,0),
                                            DilutionPoint = 1:10))


  dat7 <- data.table::data.table(data.frame(groupIndices = rep(1,10),
                                            ID = rep(1,10),
                                            REPLICATE = rep(1,10),
                                            X = 1:10,
                                            Y = c(1,1,6,2,3,4,5,6,6,6),
                                            DilutionPoint = 1:10))



  y = "Y"
  x = "X"

  testthat::expect_equal(chooseModel(dat1, tidyselect::all_of(y), tidyselect::all_of(x), model)[[1]]$model.name, "linear")
  #testthat::expect_equal(chooseModel(dat1, tidyselect::all_of(y), tidyselect::all_of(x), model)[[1]]$cor, 1)

  testthat::expect_equal(chooseModel(dat2, tidyselect::all_of(y), tidyselect::all_of(x), model)[[1]]$model.name, "quadratic")
  #testthat::expect_equal(round(chooseModel(dat2, tidyselect::all_of(y), tidyselect::all_of(x), model)[[1]]$cor, 3), 0.962)

  testthat::expect_equal(chooseModel(dat3, tidyselect::all_of(y), tidyselect::all_of(x), model)[[1]]$model.name, "logistic")
  #testthat::expect_equal(round(chooseModel(dat4, tidyselect::all_of(y), tidyselect::all_of(x), model)[[1]]$cor, 3), 0.992)

  testthat::expect_true(chooseModel(dat4, tidyselect::all_of(y), tidyselect::all_of(x), model)[[1]]$dat$outlier[4])

  testthat::expect_equal(chooseModel(dat5, tidyselect::all_of(y), tidyselect::all_of(x), model)[[1]]$model.name, "linear")
  testthat::expect_true(chooseModel(dat5, tidyselect::all_of(y), tidyselect::all_of(x), model)[[1]]$dat$outlier[1])

  testthat::expect_equal(chooseModel(dat6, tidyselect::all_of(y), tidyselect::all_of(x), model)[[1]]$model.name, "quadratic")
  testthat::expect_true(chooseModel(dat6, tidyselect::all_of(y), tidyselect::all_of(x), model)[[1]]$dat$outlier[9])

  testthat::expect_equal(chooseModel(dat7, tidyselect::all_of(y), tidyselect::all_of(x), model)[[1]]$model.name, "logistic")
  testthat::expect_true(chooseModel(dat7, tidyselect::all_of(y), tidyselect::all_of(x), model)[[1]]$dat$outlier[3])



})


# testthat::test_that("checks function outlier", {
#
#   dat1 <- data.table::data.table(data.frame(groupIndices = rep(1,10),
#                                             ID = rep(1,10),
#                                             REPLICATE = rep(1,10),
#                                             X = 1:10,
#                                             Y = 1:10,
#                                             DilutionPoint = 1:10))
#
#   dat2 <- data.table::data.table(data.frame(groupIndices = rep(1,10),
#                                             ID = rep(1,10),
#                                             REPLICATE = rep(1,10),
#                                             X = 1:10,
#                                             Y = c(1:3,10,5:10),
#                                             DilutionPoint = 1:10))
#
#   dat3 <- data.table::data.table(data.frame(groupIndices = rep(1,10),
#                                             ID = rep(1,10),
#                                             REPLICATE = rep(1,10),
#                                             X = 1:10,
#                                             Y = c(10, 2:10),
#                                             DilutionPoint = 1:10))
#
#
#   y = "Y"
#   x = "X"
#   model=c("logistic", "linear", "quadratic")
#
#
#   MODEL <- chooseModel(dat1, tidyselect::all_of(y), tidyselect::all_of(x), model)
#   MODEL2 <- chooseModel(dat2, tidyselect::all_of(y), tidyselect::all_of(x), model)
#   MODEL3 <- chooseModel(dat3, tidyselect::all_of(y), tidyselect::all_of(x), model)
#
#
#   testthat::expect_vector(outlier(dat = dat1, modelObject =  MODEL[[1]][[2]], res = 2, count = 1)$outlier, ptype = logical(), size = 10)
#   testthat::expect_false(all(outlier(dat = dat1, modelObject =  MODEL[[1]][[2]], res = 2, count = 1)$outlier))
#
#   testthat::expect_vector(outlier(dat = dat2, modelObject =  MODEL2[[1]][[2]], res = 2, count = 1)$outlier, ptype = logical(), size = 10)
#   testthat::expect_true(outlier(dat = dat2, modelObject =  MODEL2[[1]][[2]], res = 2, count = 1)$outlier[4])
#
#   testthat::expect_true(outlier(dat = dat3, modelObject =  MODEL3[[1]][[2]], res = 2, count = 1)$outlier[1])
#
# })
