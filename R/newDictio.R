#' @importFrom utils write.table
newDictio <- function(x) {
  taxa.temp <- as.character(x$Taxa)
  taxa.temp <- trimws(taxa.temp)
  taxa.temp <- gsub(" ", "_", taxa.temp, fixed = T)
  taxa.temp <- as.data.frame(c(length(taxa.temp), taxa.temp))
  write.table(taxa.temp, "custom_dictio.dic", col.names = F, row.names = F, quote = F)
  file.create(paste0(getwd(), "/custom_dictio.aff"), showWarnings = FALSE)
}
