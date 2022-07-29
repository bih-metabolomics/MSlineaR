#' Title
#'
#' @param dat
#' @param x
#' @param y
#'
#' @return
#' @export
#'
#' @examples
plotSignals <- function(dat, x="DP", y = "IntensityNorm"){
  #browser()
  g <- dat %>%
    ggplot(aes(get(x), get(y))) +
    geom_point(aes(color = color)) +#, shape = pch)) +
    scale_shape_identity() +
    scale_colour_identity() +
    geom_hline(aes(yintercept = get(y)), data = . %>% filter(str_detect(Comment, "last|first")), linetype = "dashed", color = "blue") +
    theme_bw() +
    labs(x = x, y = y) +
    scale_x_continuous(labels = c(min(dat[x]): max(dat[x])), breaks = seq(min(dat[x]), max(dat[x]))) +
    facet_wrap(~ID, scales = "free_y", ncol = 5)

  if (any(is.na(dat[y]))){
    g <- g + geom_vline(data = dat %>% filter(is.na(get(y))), aes(xintercept = get(x), color = "darkgrey"))
  }

  if("enoughPoints"  %in% colnames(dat)){
    if (any(!dat$enoughPoints)){
      g <- g + geom_segment(aes(y = 0, x = min(dat[x]), yend = max(dat[y]), xend = max(dat[x]), color = "red"), data = . %>% filter(enoughPoints %in% FALSE))
    }
  }

  if("fitted"  %in% colnames(dat)){
    if (any(!is.na(dat$fitted))){
      datsub <- dat |>
        dplyr::group_by(ID) |>
        dplyr::mutate( start = DilutionPoint %in% linearRangeStart) |>
        dplyr::mutate( end = DilutionPoint %in% unique(linearRangeEnd))

      g <- g +
        geom_vline(data = datsub |> filter(start %in% TRUE), aes(xintercept = .data[[x]]), linetype = "dashed") +
        geom_vline(data = datsub |> filter(end %in% TRUE), aes(xintercept = .data[[x]]), linetype = "dashed") +
        geom_line(aes(x = get(x), y = fittedLM), col = "darkgreen") +
        geom_line(aes(x = get(x), y = fitted), col = "blue") +
        ylim(0,max(dat[[y]]) + 20)

    }
  }

  return(g)
}
