#' @name get_reference_db
#'
#' @title Function to download reference taxonomic trees, traits dataset,
#' and biomonitoring index scores.
#'
#' @description
#' This function allows the user to download several data sets useful for `biomonitoR` analyses.
#' See details to discover which data sets you can download.
#'
#' @param db_name dataset name to download. See details to know the available data sets.
#'
#'
#' @details Available data set:
#' 1) mi_cuba: Taxonomic data set of macroinvertebrates from Cuba (Torres-Cambas et al., 2023).
#'
#' @keywords reference data sets download
#'
#' @export
#'
#' @examples
#' db_reference <- get_reference_db(db_name = "mi_cuba_2023")
#'
#' @importFrom utils read.csv

get_reference_db <- function(db_name = "mi_cuba_2023"){

  error_message <-  db_name %in% dfMain$dataset

  if(!error_message){
    stop("Please provide a valid db_name")
  }

  db <- read.csv(dfMain$link[dfMain$dataset == db_name])

  # message(paste("Please to cite ‘biomonitoR’ in your publications use:\n
  #               Laini A., Guareschi S., Bolpagni R., Burgazzi G., Bruno D., Gutiérrez-Cánovas C.,
  #               ... & Cancellario T. (2022). biomonitoR: an R package for managing ecological
  #               data and calculating biomonitoring indices. PeerJ, 10, e14183. \n
  #               Moreover to cite the reference dataset use: \n \n",
  #               dfMain$citation_APA[dfMain$dataset == db_name]))

  db[is.na(db)] <- ""
  return(db)
}
