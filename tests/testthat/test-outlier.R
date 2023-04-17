testthat::test_that("function chooseModel works",{
  dat1 <- data.table::data.table(groupIndices = rep(1,10),
                                            ID = rep(1,10),
                                            IDintern = paste0("s", 1:10),
                                            REPLICATE = rep(1,10),
                                            color = "black",
                                            X = 1:10,
                                            Y = 1:10,
                                            DilutionPoint = 1:10)

  dat2 <- data.table::data.table(groupIndices = rep(1,10),
                                            ID = rep(1,10),
                                            IDintern = paste0("s", 1:10),
                                            REPLICATE = rep(1,10),
                                            color = "black",
                                            X = 1:10,
                                            Y = c(1:5,4:0),
                                            DilutionPoint = 1:10)

  dat3 <- data.table::data.table(groupIndices = rep(1,10),
                                            ID = rep(1,10),
                                            IDintern = paste0("s", 1:10),
                                            REPLICATE = rep(1,10),
                                            color = "black",
                                            X = 1:10,
                                            Y = c(1,1,1,2,3,4,5,6,6,6),
                                            DilutionPoint = 1:10)

  dat4 <- data.table::data.table(groupIndices = rep(1,10),
                                            ID = rep(1,10),
                                            IDintern = paste0("s", 1:10),
                                            REPLICATE = rep(1,10),
                                            color = "black",
                                            X = 1:10,
                                            Y = c(1:3,10,5:10),
                                            DilutionPoint = 1:10)

  dat5 <- data.table::data.table(groupIndices = rep(1,10),
                                            ID = rep(1,10),
                                            IDintern = paste0("s", 1:10),
                                            REPLICATE = rep(1,10),
                                            color = "black",
                                            X = 1:10,
                                            Y = c(10, (2:10)+0.2),
                                            DilutionPoint = 1:10)

  dat6 <- data.table::data.table(groupIndices = rep(1,10),
                                            ID = rep(1,10),
                                            IDintern = paste0("s", 1:10),
                                            REPLICATE = rep(1,10),
                                            color = "black",
                                            X = 1:10,
                                            Y = c(1:5,4:2,8,0),
                                            DilutionPoint = 1:10)


  dat7 <- data.table::data.table(groupIndices = rep(1,10),
                                            ID = rep(1,10),
                                            IDintern = paste0("s", 1:10),
                                            REPLICATE = rep(1,10),
                                            color = "black",
                                            X = 1:10,
                                            Y = c(1,1,6,2,3,4,5,6,6,6),
                                            DilutionPoint = 1:10)


  model. = c("logistic", "linear", "quadratic")
  y. = "Y"
  x. = "X"
  SDRES_MIN. = 1
  STDRES. = 2
  abbr. = "FOD"

  testthat::expect_equal(chooseModel(
    dats = dat1,
    y = y.,
    x = x.,
    model =  model.,
    SDRES_MIN = SDRES_MIN.,
    STDRES = STDRES.,
    abbr = abbr.)[[1]]$model.name, "linear")

  testthat::expect_equal(which(!chooseModel(
    dats = dat1,
    y = y.,
    x = x.,
    model =  "linear",
    SDRES_MIN = SDRES_MIN.,
    STDRES = STDRES.,
    abbr = abbr.)[[1]]$model$RMSE %in% NA), 2)

  testthat::expect_equal(chooseModel(
    dats = dat2,
    y = y.,
    x = x.,
    model =  model.,
    SDRES_MIN = SDRES_MIN.,
    STDRES = STDRES.,
    abbr = abbr.)[[1]]$model.name, "quadratic")

  testthat::expect_equal(chooseModel(
    dats = dat3,
    y = y.,
    x = x.,
    model =  model.,
    SDRES_MIN = SDRES_MIN.,
    STDRES = STDRES.,
    abbr = abbr.)[[1]]$model.name, "logistic")

  testthat::expect_true(chooseModel(
    dats = dat4,
    y = y.,
    x = x.,
    model =  model.,
    SDRES_MIN = SDRES_MIN.,
    STDRES = STDRES.,
    abbr = abbr.)[[1]]$dat$OutlierFOD[4])

  testthat::expect_equal(chooseModel(
    dats = dat5,
    y = y.,
    x = x.,
    model =  model.,
    SDRES_MIN = SDRES_MIN.,
    STDRES = STDRES.,
    abbr = abbr.)[[1]]$model.name, "linear")

  testthat::expect_true(chooseModel(
    dats = dat5,
    y = y.,
    x = x.,
    model =  model.,
    SDRES_MIN = SDRES_MIN.,
    STDRES = STDRES.,
    abbr = abbr.)[[1]]$dat$OutlierFOD[1])

  testthat::expect_equal(sum(chooseModel(
    dats = dat5,
    y = y.,
    x = x.,
    model =  model.,
    SDRES_MIN = SDRES_MIN.,
    STDRES = STDRES.,
    abbr = abbr.)[[1]]$dat$OutlierFOD), 1)

  testthat::expect_equal(chooseModel(
    dats = dat6,
    y = y.,
    x = x.,
    model =  model.,
    SDRES_MIN = SDRES_MIN.,
    STDRES = STDRES.,
    abbr = abbr.)[[1]]$model.name, "quadratic")

  testthat::expect_true(chooseModel(
    dats = dat6,
    y = y.,
    x = x.,
    model =  model.,
    SDRES_MIN = SDRES_MIN.,
    STDRES = STDRES.,
    abbr = abbr.)[[1]]$dat$OutlierFOD[9])

  testthat::expect_equal(chooseModel(
    dats = dat7,
    y = y.,
    x = x.,
    model =  model.,
    SDRES_MIN = SDRES_MIN.,
    STDRES = STDRES.,
    abbr = abbr.)[[1]]$model.name, "logistic")

  testthat::expect_true(chooseModel(
    dats = dat7,
    y = y.,
    x = x.,
    model =  model.,
    SDRES_MIN = SDRES_MIN.,
    STDRES = STDRES.,
    abbr = abbr.)[[1]]$dat$OutlierFOD[3])

})

