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
plotSignals <- function(dat, x="Concentration", y = "Intensity", Replicate = "Batch", nrcol = 2){
  #browser()
  facetNames = paste(dat$groupIndices, dat[[Replicate]], dat$mz, sep = "_")
  names(facetNames) = dat$groupIndices
  g <- dat %>%
    ggplot(aes(get(x), get(y))) +
    geom_point(aes(color = color)) +#, shape = pch)) +
    scale_shape_identity() +
    scale_colour_identity() +
    geom_hline(aes(yintercept = get(y)), data = . %>% filter(str_detect(Comment, ": last|: first")), linetype = "dashed", color = "blue") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = x, y = y) +
    scale_x_continuous(labels = c(min(dat[[x]]): max(dat[[x]])), breaks = seq(min(dat[[x]]), max(dat[[x]]))) +
    facet_wrap(~groupIndices,  ncol = nrcol,scales = "free_y",
               labeller = labeller(groupIndices = facetNames))

  if (any(is.na(dat[[y]]))){
    g <- g + geom_vline(data = dat %>% filter(is.na(get(y))), aes(xintercept = get(x), color = "darkgrey"))
  }

  if("enoughPeaks"  %in% colnames(dat)){
    if (any(!dat$enoughPeaks %in% TRUE)){
      g <- g + geom_segment(aes(y = 0, x = min(dat[[x]]), yend = max(dat[[y]], na.rm = T), xend = max(dat[[x]]), color = "red"),
                            data = . %>% filter(!enoughPeaks %in% TRUE))
    }
  }
  if("aboveCorFit"  %in% colnames(dat)){
    if (any(dat$aboveCorFit %in% FALSE)){
      g <- g + geom_segment(aes(y = 0, x = min(dat[[x]]), yend = max(dat[[y]], na.rm = T), xend = max(dat[[x]]), color = "red"),
                            data = . %>% filter(aboveCorFit %in% FALSE))
    }
  }



  if("IslinearRange"  %in% colnames(dat)){
    if (any(!is.na(dat$IslinearRange))){
      datsub <- dat |>
        dplyr::group_by(groupIndices) |>
        dplyr::mutate( start = DilutionPoint %in% linearRangeStart,
                       end = DilutionPoint %in% unique(linearRangeEnd),
                       startpos = DilutionPoint %in% min(DilutionPoint[IsPositivAssociated %in% TRUE]),
                       endpos = DilutionPoint %in% max(DilutionPoint[IsPositivAssociated %in% TRUE]))

      g <- g +
        geom_vline(data = datsub |> group_by(groupIndices) |> filter(start %in% TRUE), aes(xintercept = .data[[x]]), linetype = "dashed") +
        geom_vline(data = datsub |> group_by(groupIndices) |> filter(end %in% TRUE), aes(xintercept = .data[[x]]), linetype = "dashed") +
        geom_line(aes(x = get(x), y = ablineFit), col = "orange") +
        geom_vline(data = datsub |> group_by(groupIndices)|> filter(startpos %in% TRUE), aes(xintercept = .data[[x]]), linetype = "dotted", col = "darkgrey") +
        geom_vline(data = datsub |> group_by(groupIndices)|> filter(endpos %in% TRUE), aes(xintercept = .data[[x]]), linetype = "dotted", col = "darkgrey") #+
        #xlim(0, max(dat[[x]]))

        #geom_line(aes(x = get(x), y = modelFit), col = "blue") +
        #ylim(0,max(dat[[y]]) + 20)

    }
  }

  return(g)
}
