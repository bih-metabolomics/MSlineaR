#' serviceTools.
#'
#' @name serviceTools
#' @docType package
NULL

p_trans <- function(p){
  ifelse(p < 0.001, "<0.001", p)
}


##' Builds the inference table 
##'
##' This function builds the inference table for beta, upper, and lower 
##' @title Build inference table
##' @param infStat Data.frame of beta, lower, upper, and p. The column names must include BETA, LOWER, UPPER and P)
##' @return data.frame
##' @author Jochen Kruppa
##' @export
buildInference <- function(infStat, pCol = "p", id_row = FALSE){
    ## build up the inference
    betaCol <- grep("beta", names(infStat))
    lowerCol <- grep("lower", names(infStat))
    upperCol <- grep("upper", names(infStat))
    infStat[betaCol] <- ifelse(infStat[,betaCol] > 1000, Inf, infStat[,betaCol])
    infStat[upperCol] <- ifelse(infStat[,upperCol] > 1000, Inf, infStat[,upperCol])
  if(!id_row){
    inference <- data.frame(OR = paste0(
                              sprintf("%.5f", infStat[,betaCol]), " [",
                              sprintf("%.5f", infStat[,lowerCol]), ", ",
                              sprintf("%.5f", infStat[,upperCol]), "]"),
                            p = infStat[,pCol],
                            row.names = row.names(infStat))
  } else {
    inference <- data.frame(id = infStat$id,
                            OR = paste0(
                              sprintf("%.5f", infStat[,betaCol]), " [",
                              sprintf("%.5f", infStat[,lowerCol]), ", ",
                              sprintf("%.5f", infStat[,upperCol]), "]"),
                            p = infStat[,pCol])
  }
  return(inference)
}


##' Generate output table for pandoc
##'
##' This function generates a descriptive table and adds the inference
##' statistics (must be given). Due to many changes this function is really big.
##' @title Descriptive and inference statistics for a given data set 
##' @param data Data set of descriptive variables. No response!
##' @param infStat Inference data set with the column names, including: id, beta, lower, upper, p 
##' @param continous_flag Which of the variables should be continous?
##' If not given, automatically decision by cate_threshold
##' @param cate_threshold Cutoff of classes when a variable should be considered has continous?
##' @param grouped.by Should the descriptive statistics be grouped by this factor?
##' @param ratio.by.col Should be the ratio of the grouped variable calculated by column?
##' @param caption Caption of the generated pandoc table
##' @param colnames The column names of the generated pandoc table
##' @return pandoc.table
##' @author Jochen Kruppa
##' @export
##' @examples
##' data(iris)
##' irisInfStat <- data.frame(id = names(iris)[-5],
##'                           beta = c(1, 2, 3, 4),
##'                           lower.conf = c(0, 1, 2, 3),
##'                           upper.conf = c(2, 3, 4, 5),
##'                           p = c(0.11111, 0.22222, 0.44443, 0.4434343))
##' prepareDeskStatTable(data = iris, infStat = irisInfStat)
##' https://cran.r-project.org/web/packages/kableExtra/vignettes/awesome_table_in_html.html

## pacman::p_load(tidyverse, magrittr)
## infStat <- data.frame(id = names(iris),
##                       beta = c(1, 2, 3, 4, 4, 4, 5, 5),
##                       lower.conf = c(0, 1, 2, 3, 1, 1, 1, 5),
##                       upper.conf = c(2, 3, 4, 5, 10, 10, 10, 5),
##                       p = c(0.11111, 0.22222, 0.44443, 0.4434343, 0, 0, 0, 0))
## iris %<>% mutate(POD = sample(c(0,1,NA), size = nrow(iris), replace = TRUE),
##                  sex = sample(c("man", "woman", NA),
##                               size = nrow(iris), replace = TRUE),
##                  sep_na = sample(c(Sepal.Length, rep(NA, 100)),
##                                  size = nrow(iris), replace = TRUE))

## sum(is.na(iris$sep_na))
## tally(sep_na ~ POD, iris)

## require(mosaic)
## tally(iris[,6] ~ iris[,7], iris)
## grouped.by <- as.factor(iris$POD) 

## prepareDeskStatTable(iris, infStat)

prepareDeskStatTable <- function(data,
                                 infStat = NULL,
                                 continous_flag = NULL,
                                 cate_threshold = 12,
                                 grouped.by = NULL,
                                 ratio.by.col = TRUE,
                                 caption = "This is the caption",
                                 colnames = NULL,
                                 color_row = FALSE){
    require(plyr)
    require(mosaic)
    require(pander)
    require(knitr)
    require(kableExtra)
    ## check the infStat data for consistency
    if(!is.null(infStat)){
        id_flag <- any(grepl("id", names(infStat)))
        if(id_flag){
            row.names(infStat) <- gsub("\\.", "-", infStat[,grep("id", names(infStat))])
            names(data) <- gsub("\\.", "-", names(data))
        } else {
            stop("The infStat data must have a 'id'-named column, which has the same entries like names(data)")
        }
    }
    ## find the continuous variables by threshold if not given
    if(is.null(continous_flag)){
        continous_flag <- laply(data, function(x) length(unique(x)) >= cate_threshold)
    } else {
        continous_flag <- continous_flag
    }
    ## must we group by a treatment or other?
    if(is.null(grouped.by)){
        if(any(continous_flag)){
            ## get mean and sd for continous data
            continousStats <- ldply(data[continous_flag],
                                    function(x) paste0(sprintf("%.2f",round(mean(x, na.rm = TRUE), 1)), " +/- ",
                                                       sprintf("%.2f",round(sd(x, na.rm = TRUE),1)), " (",
                                                       sprintf("%.0f",round(min(x, na.rm = TRUE),1)), ", ",
                                                       sprintf("%.0f",round(max(x, na.rm = TRUE),1)), ")"
                                                       ))
            row.names(continousStats) <- continousStats$.id
            continousStats$.id <- NULL
        } else {
            continousStats <- data.frame()
        }
        if(any(!continous_flag)){
            ## get the category tables
            makeCategoryTable <- function(var){
                count <- ldply(table(var))[,2]
                relVal <- round(ldply(table(var)/sum(table(var)))[,2]*100,0)
                return(c(paste0(count, " (", relVal, "%)"), sum(is.na(var))))
            }   
            categoryTable <- llply(data[!continous_flag],
                                   function(x) data.frame(V1 = makeCategoryTable(x),
                                                          row.names = c(ldply(table(x))[,1], "ZZZMissing")))
        } else {
            categoryTable <- list()
        }
    } else {
        if(!is.factor(grouped.by) | is.null(levels(grouped.by)))
            stop("grouped.by must be a factor with levels! Use: factor(x, labels = c('foo', 'bar'))")
        if(any(continous_flag)){
          ## get the *grouped* continous tables
          continousStats <- ldply(data[continous_flag],
                                  function(x) {
                                    small_df <- data.frame(x, grouped.by)
                                    groupedContStatDf <- ddply(small_df, ~grouped.by, summarise,
                                                               mean = sprintf("%.2f",
                                                                              round(mean(x, na.rm = TRUE), 1)),
                                                               sd = sprintf("%.2f", round(sd(x, na.rm = TRUE), 1)))
                                    ## if(any(is.na(x))){
                                    ##   cont_miss <- tally(x ~ grouped.by, small_df)
                                    ##   miss_cont_vec <- as.numeric(cont_miss[nrow(cont_miss),])
                                    ## } else {
                                    ##   miss_cont_vec <- rep("0", length(levels(grouped.by)) + 1)
                                    ## }
                                    if(any(is.na(x)) & !any(is.na(grouped.by))){
                                      cont_miss <- tally(x ~ grouped.by, small_df)
                                      miss_cont_vec <- as.numeric(cont_miss[nrow(cont_miss),])
                                    }
                                    if(any(is.na(x)) & any(is.na(grouped.by))) {
                                      miss_cont_vec <- rep("0", length(levels(grouped.by)) + 1)
                                    }
                                    if(!any(is.na(x)) & !any(is.na(grouped.by))) {
                                      miss_cont_vec <- rep("0", length(levels(grouped.by)))
                                    }
                                    if(!any(is.na(x)) & any(is.na(grouped.by))) {
                                      miss_cont_vec <- rep("0", length(levels(grouped.by)) + 1)
                                    }
                                    ## paste0(groupedContStatDf[,2], " +/- ", groupedContStatDf[,3])
                                    tmp_df <- rbind(paste0(groupedContStatDf[,2], " +/- ",
                                                           groupedContStatDf[,3]),
                                                    miss_cont_vec
                                                    )
                                    row.names(tmp_df)[2] <- c("Missing")
                                    return(tmp_df)
                                  })
          continousStats$.id[duplicated(continousStats$.id)] <- str_c(unique(continousStats$.id), ".ZZZmissing")
          row.names(continousStats) <- continousStats$.id
          continousStats$.id <- NULL
          names(continousStats) <- levels(grouped.by)
          if(any(is.na(names(continousStats)))){
            miss_vec <- apply(continousStats[str_detect(row.names(continousStats), "ZZZ"),],
                              1, function(x) sum(as.numeric(x)))
            tmp_na_df <- continousStats[which(is.na(names(continousStats)))]
            tmp_na_df[str_detect(tmp_na_df[,1], "\\+\\/\\-"), ] <- miss_vec 
            continousStats[which(is.na(names(continousStats)))] <- tmp_na_df
          }
        } else {
            continousStats <- data.frame()
        }
        if(any(!continous_flag)){
            ## get the *grouped* category tables
            categoryTable <- llply(data[!continous_flag], function(x) {
                absVal <- cbind(table(x, grouped.by, useNA = "always"))
                denom <- rep(apply(cbind(table(x, grouped.by, useNA = "always")), 2, sum),
                             each = nrow(absVal))
                relVal <- cbind(table(x, grouped.by, useNA = "always"))/denom
                naVal <- ddply(data.frame(x, grouped.by), ~grouped.by, summarise, sum(is.na(x)))
                if(ratio.by.col == TRUE){
                    catTable <- ldply(1:ncol(absVal), function(i){
                        paste0(absVal[,i], " (", round(relVal[,i]*100, 0), "%)")})
                } else {
                    catTable <- ldply(1:ncol(absVal), function(i){
                        paste0(absVal[,i], " (", round(relVal[i,]*100, 0), "%)")})
                }
                catTableOut <- data.frame(rbind(t(catTable)## , naVal[,2]
                                                ))
                row.names(catTableOut) <- c(na.omit(row.names(absVal)), "ZZZMissing")
                names(catTableOut) <- levels(grouped.by)
                ## complex na check last column
                na_check_last_column <- catTableOut[, length(levels(grouped.by)) + 1] %>%
                  str_replace("(^\\d+)\\s.*", "\\1") %>%
                  as.numeric  %>%
                  sum
                if(!na_check_last_column){
                  catTableOut <- catTableOut[,-(length(levels(grouped.by)) + 1)]
                }
                return(catTableOut)
            })                
        }else {
            categoryTable <- list()
        }
    }
    ## put the data together
    preData <- rbind(continousStats, do.call(rbind, categoryTable))
    preData$idx_sort <- row.names(preData)
    preData$idx <- laply(preData$idx_sort, function(x) strsplit(x, "\\.")[[1]][1])
    ## add the inference
    ## build up the inference
    if(is.null(infStat)){
        inference <- data.frame(OR = rep(NA, ncol(data)), p = NA, row.names = names(data))
    } else {
        inference <- buildInference(infStat)
    }
    ## inference$idx <- gsub("2", "", row.names(inference))
    inference$idx <- row.names(inference)
    ## merge everything
    outDf <- merge(preData, inference, by.x = "idx", by.y = "idx", all.x = TRUE)
    ## prepare for ordering
    orderDf <- data.frame(orgNames = outDf$idx_sort,
                          sortNames = tolower(outDf$idx_sort),
                          idx = 1:length(outDf$idx))
    orderDf <- orderDf[order(orderDf$sortNames),]
    ##finalize
    outDf$Level <- laply(outDf$idx_sort, function(x) strsplit(x, "\\.")[[1]][2])
    if(is.null(grouped.by)){
      outDf <- outDf[orderDf$idx, c("idx", "Level", "V1", "OR", "p")]
    } else {
      if(any(is.na(grouped.by))) {
        levels(grouped.by) <- c(levels(grouped.by), "missing")
        names(outDf)[which(is.na(names(outDf)))] <- "missing"
      }      
      outDf <- outDf[orderDf$idx, c("idx", "Level", levels(grouped.by), "OR", "p")]
    }
    outDf$Level <- gsub("ZZZ", "", outDf$Level)
    outDf$OR[duplicated(outDf$idx)] <- NA
    outDf$p[duplicated(outDf$idx)] <- NA
    outDf$p <- sprintf("%.3f", round(outDf$p,3))
    outDf$p <- ifelse(outDf$p == "0.000", "<0.001", outDf$p)
    outDf$idx[duplicated(outDf$idx)] <- NA
    row.names(outDf) <- NULL
    outDf$idx_sort <- NULL
    ## remove all NA and make the data.frame clean!
    outDf$OR <- as.character(outDf$OR)
    outDf$p[outDf$p == "NA"] <- ""
    outDf[is.na(outDf)] <- ""
    ## write outDf out
    ##
    if(is.null(grouped.by)){
        if(is.null(colnames)){
            names(outDf) <- c("Risc factor", "Level", "Desc. Stat.", "OR", "p")
        } else {
            names(outDf) <- colnames
        }
        if(is.null(infStat)){
            outDf <- outDf[,1:3]
            ## pandoc.table(outDf,
            ##              justify = c('left', 'right', 'right'),
            ##              caption = caption,
            ##              row.names = FALSE,
            ##              keep.trailing.zeros = TRUE,
            ##              style = 'rmarkdown')
            kable(outDf,
                  align = c('l', 'r', 'r'),
                  digits = 3,
                  caption = caption) %>%
              kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                        full_width = FALSE)
        } else {
            ## kable(outDf, align = c('l', 'r', 'r', 'r', 'r'), digits = 3)
            ## pandoc.table(outDf,
            ##              justify = c('left', 'right', 'right', 'right', 'right'),
            ##              caption = caption,
            ##              row.names = FALSE,
            ##              keep.trailing.zeros = TRUE,
            ##              style = 'rmarkdown')
            kable(outDf,
                  align = c('l', 'r', 'r', 'r', 'r'),
                  digits = 3,
                  caption = caption) %>%
              kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                        full_width = FALSE)
        }
    } else {
        if(is.null(colnames)){
            names(outDf) <- c("Risc factor", "Level", paste("Stat", levels(grouped.by)), "OR", "p")
        } else {
            names(outDf) <- colnames
        }
        ## kable(outDf, align = c('l', 'r', rep('r', times = length(levels(grouped.by))), 'r', 'r'),
        ##       digits = 3)
        ## pandoc.table(outDf,
        ##              justify = c('left', 'right', rep('right', times = length(levels(grouped.by))),
        ##                  'right', 'right'),
        ##              caption = caption,
        ##              row.names = FALSE,
        ##              keep.trailing.zeros = TRUE,
        ##              style = 'rmarkdown')
        k_out <- kable(outDf,
                       align = c('l', 'r', rep('r', times = length(levels(grouped.by))), 'r', 'r'),
                       digits = 3,
                       caption = caption) %>%
          kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                        full_width = FALSE)
        if(!color_row){
          return(k_out)
        } else {
          k_out %>%
            row_spec(which(outDf$p < 0.05), bold = TRUE, color = "white",
                     background = "#D7261E")
        }
    }
}
