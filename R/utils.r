#' @title Olink Dist Plot
#' @description Create a Distance plot
#' @param npx NPX dataset (output from OlinkAnalyze::read_NPX)
#' @param panel Panel Name
#' @param filename filename, the suffix will determine the filetype
#' @param width Figure width,  Default: 7
#' @param height Figure height,  Default: 7
#' @return distance plot
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  olink_dist(plasma_npx, "Olink CARDIOVASCULAR III", filename="test.png")
#'  }
#' }
#' @seealso
#'  \code{\link[OlinkAnalyze]{olink_dist_plot}}
#' @rdname olink_dist
#' @export
#' @import OlinkAnalyze
#' @import ggplot2
#' @import dplyr
olink_dist <- function(npx, panel, filename, width = 7, height = 7) {
  olink_dist_plot(npx %>% dplyr::filter(Panel == panel)) +
  scale_fill_manual(values=c("green", "red")) +
  theme(text = element_text(size=6),
        axis.text.x = element_text(angle = 90, vjust = 0.5))
  ggsave(filename, width = width, height = height)
}


