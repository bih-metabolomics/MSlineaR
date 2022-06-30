
kable_style <- function(tbl) {
  kable_html <- tbl %>%
    kable %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                  full_width = FALSE)
  return(kable_html)
}



##' A wrapper how saves a list to a xls file
##'
##' A wrapper how saves a list to a xls file
##' @title A wrapper how saves a list to a xls file
##' @param list The data frames in a list which should be written
##' @param file Save location
##' @param oneFile One file with sheets or many files?
##' @return NULL
##' @author Jochen Kruppa
##' @export
write_xlsx <- function(list, file, oneFile = TRUE){
    require(openxlsx)
    wb <- createWorkbook()
    if (is.null(names(list))) 
        stop("The list must be named by: names(list) <- c('nameList1',...)")
    if (oneFile) {
        for (i in seq_along(names(list))) {
            addWorksheet(wb, names(list)[i])
            writeData(wb, sheet = i, x = list[[i]], colNames = TRUE)
        }
        saveWorkbook(wb, file, overwrite = TRUE)
    } else {
        if(length(file) != length(list)) stop ("If you do not want 'oneFile' provide number of files equal list length")
        for (i in seq_along(names(list))) {
            wb <- createWorkbook()
            addWorksheet(wb, "Results")
            writeData(wb, sheet = 1, x = list[[i]], colNames = TRUE)
            saveWorkbook(wb, file[i], overwrite = TRUE)
        }
    }
}


##' Message function with time
##'
##' Message function with time
##' @title Message function with time 
##' @param ... 
##' @return Message with time
##' @author Jochen Kruppa
##' @export
talk <- function(...) {
  require(stringr)
  message(str_c(TIME(), " ..... "), ...)
}


##' Lazy wrapper for system
##'
##' Lazy wrapper for system
##' @title Lazy wrapper for system
##' @param cmd pasted command to run
##' @return NULL
##' @author Jochen Kruppa
##' @export
runCMD <- function(cmd, intern = FALSE) {
    try(system(cmd, wait = TRUE, intern = intern))
}

##' Small wrapper for dir.create
##'
##' Small wrapper for dir.create
##' @title Small wrapper for dir.create 
##' @param dir_path 
##' @return NULL
##' @author Jochen Kruppa
##' @export
mk_dir <- function(dir_path){
  if(!file.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
}


##' Time function
##'
##' Time function
##' @title Time function 
##' @param ... 
##' @return Time in specified format
##' @author Jochen Kruppa
TIME <- function(...) format(Sys.time(), "%b %d %X")


##' Alias for glimpse
##'
##' Alias for glimpse
##' @title Alias for glimpse
##' @param ... 
##' @return NULL
##' @author Jochen Kruppa
##' @export
gli <- function(...) {
  require(dplyr)
  glimpse(...) 
}

##' Small wrapper function for dir
##'
##' Small wrapper function for dir
##' @title Small wrapper function for dir
##' @param ... 
##' @return file.path
##' @author Jochen Kruppa
##' @export
dir0 <- function (...){ 
  dir(..., full.names = TRUE, recursive = TRUE)
}

##' Wrapper function for log reg
##'
##' Wrapper function for log reg
##' @title Wrapper function for log reg
##' @param fit 
##' @param kable_lgl 
##' @return 
##' @author Jochen Kruppa
print_log_reg <- function(fit, kable_lgl = FALSE) {
  p_trans <- function(p){
    p <- sprintf("%.3f", round(p,3))
    p <- ifelse(p == "0.000", "<0.001", p)
    ifelse(p < 0.001, "<0.001", p)
  }
  ## get the coefficients
  sum_fit <- summary(fit)
  ## get the confidence intervals
  lower.conf <- try(exp(confint(fit))[,1][-1])
  if(class(lower.conf) == "try-error")
    lower.conf <- 0
  upper.conf <- try(exp(confint(fit))[,2][-1])
  if(class(upper.conf) == "try-error")
    upper.conf <- 0
  ## get the summary data
  coef_tbl <- coef(sum_fit)[,c(1,4)] %>%
    as.data.frame() %>%
  extract(-1,) %>%
    mutate(risk_factor = rownames(.),
           Estimate = exp(Estimate),
           lower = round(lower.conf, 3),
           upper = round(upper.conf, 3),
           CI = str_c("[", lower, "; ", upper, "]")) %>%
    set_colnames(c("OR", "p_value", "risk_factor", "lower", "upper", "CI")) %>%
    tbl_df() %>%
    arrange(p_value) %>%
    mutate(p_value = p_trans(p_value)) %>%
    select(risk_factor, OR, CI, p_value) 
    if(kable_lgl){
    kable_print <- coef_tbl %>%
      kable(.,
            align = c('l', 'r', 'r', 'r'),
            digits = 5,
            caption = "") %>%
      kable_styling(kable_input = ., bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                    full_width = FALSE)
    return(kable_print)
  } else {
    return(coef_tbl)
  }
}

