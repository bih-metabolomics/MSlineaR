library(pacman)
pacman::p_load(MASS,
               readxl,
               plyr,
               magrittr,
               knitr,
               haven,
               mosaic,
               knitr,
               kableExtra,
               grid,
               gridExtra,
               tidyverse,
               tableone,
               broom)
source("source/dataTools.R")
source("source/descStatTable.R")
rename <- function(...) dplyr::rename(...)
select <- function(...) dplyr::select(...)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

main_dir <- "C:/Users/wiebachj/Desktop/"
project_dir <- file.path(main_dir, "sigmodial/scripts")

img_dir <- file.path(project_dir, "img")

data_dir <- file.path(main_dir, "sigmodial/data")
in_dir <- file.path(data_dir, "input")
proc_dir <- file.path(data_dir, "proc")
out_dir <- file.path(data_dir, "output")
