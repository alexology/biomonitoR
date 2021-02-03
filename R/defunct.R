#' Defunct functions
#'
#' \Sexpr[results=rd, stage=render]{lifecycle::badge("defunct")}
#' Executing these functions will tell you which function replaces them.
#'
#' @keywords internal
#' @name defunct
NULL

#' @export
#' @rdname defunct
asBiomonitor <- function(.variables, drop = FALSE) {
  lifecycle::deprecate_stop("0.9.1", "asBiomonitor()", "as_biomonitor()")
}


#' @export
#' @rdname defunct
aggregatoR <- function(.variables, drop = FALSE) {
  lifecycle::deprecate_stop("0.9.1", "aggregatoR()", "aggregate_taxa()")
}


#' @export
#' @rdname defunct
abuTax <- function(.variables, drop = FALSE) {
  lifecycle::deprecate_stop("0.9.1", "abuTax()", "get_taxa_abundance()")
}


#' @export
#' @rdname defunct
addBiasToTraits <- function(.variables, drop = FALSE) {
  lifecycle::deprecate_stop("0.9.1", "addBiasToTraits()", "add_bias_to_traits()")
}


#' @export
#' @rdname defunct
ambiguousSolver <- function(.variables, drop = FALSE) {
  lifecycle::deprecate_stop("0.9.1", "ambiguousSolver()", "solve_ambiguous()")
}


#' @export
#' @rdname defunct
combTaxa <- function(.variables, drop = FALSE) {
  lifecycle::deprecate_stop("0.9.1", "combTaxa()", "combine_taxa()")
}


#' @export
#' @rdname defunct
convertovegan <- function(.variables, drop = FALSE) {
  lifecycle::deprecate_stop("0.9.1", "convertovegan()", "convert_to_vegan()")
}


#' @export
#' @rdname defunct
convertobiotic <- function(.variables, drop = FALSE) {
  lifecycle::deprecate_stop("0.9.1", "convertobiotic()", "convert_to_biotic()")
}


#' @export
#' @rdname defunct
totInfo <- function(.variables, drop = FALSE) {
  lifecycle::deprecate_stop("0.9.1", "totInfo()", "general_info()")
}


#' @export
#' @rdname defunct
manageTraits <- function(.variables, drop = FALSE) {
  lifecycle::deprecate_stop("0.9.1", "manageTraits()", "manage_traits()")
}


#' @export
#' @rdname defunct
quickRename <- function(.variables, drop = FALSE) {
  lifecycle::deprecate_stop("0.9.1", "quickRename()", "quick_rename()")
}


#' @export
#' @rdname defunct
refFromTree <- function(.variables, drop = FALSE) {
  lifecycle::deprecate_stop("0.9.1", "refFromTree()", "ref_from_tree()")
}


#' @export
#' @rdname defunct
sampleTraits <- function(.variables, drop = FALSE) {
  lifecycle::deprecate_stop("0.9.1", "sampleTraits()", "sample_traits()")
}


#' @export
#' @rdname defunct
selectPcoaAxes <- function(.variables, drop = FALSE) {
  lifecycle::deprecate_stop("0.9.1", "selectPcoaAxes()", "select_pcoa_axes()")
}


#' @export
#' @rdname defunct
showscores <- function(.variables, drop = FALSE) {
  lifecycle::deprecate_stop("0.9.1", "showscores()", "show_scores()")
}


#' @export
#' @rdname defunct
ricTax <- function(.variables, drop = FALSE) {
  lifecycle::deprecate_stop("0.9.1", "ricTax()", "get_taxa_richness()")
}


#' @export
#' @rdname defunct
traitsMean <- function(.variables, drop = FALSE) {
  lifecycle::deprecate_stop("0.9.1", "traitsMean()", "average_traits()")
}


#' @export
#' @rdname defunct
traitScaling <- function(.variables, drop = FALSE) {
  lifecycle::deprecate_stop("0.9.1", "traitScaling()", "assign_traits()")
}







