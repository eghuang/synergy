#' @title Loads CSV dataset as a data.frame, dropping entries with zero.
#'
#' @description Meant to handle single-ion CSV datasets.
#'
#' @param filename String of the CSV filename.
#'
#' @details The filename includes the path from the current working directory.
#'          CSV should have a numeric "dose" column. Numeric olumns with "Z"
#'          for atomic number and "NWeights" are required for compatibility
#'          with built in modeling functions.
#'
#' @return Data.frame of the selected data file.
#'
#' @examples
#' # Data used is that in 16Chang plus later NSRL data scored before 3/31/2019.
#' ion_data <- load_ion_data("one_ion.csv")
#' mix_data <- load_ion_data("mix_ion.csv")
#' # The two .csv files contain all input HG data except Y_0
#'
#' @author Edward Greg Huang <eghuang@@berkeley.edu>
#' @export

load_ion_data <- function(filename) {
  df <- data.frame(read.csv(filename))
  return(df[!(is.na(df$dose)), ]) # Requires a numeric "dose" column
}


