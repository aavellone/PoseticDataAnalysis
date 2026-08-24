library(methods)
#' An S4 class to represent a POSet.
#'
#' @slot ptr an external pointer to C++ data
setClass("POSet", slots = list(ptr = "externalptr"))

#' An S4 class to represent a Linear POSet.
#'
#' @slot ptr an external pointer to C++ data
setClass("LinearPOSet", contains = "POSet")

#' An S4 class to represent a virtual class for POSet extention.
#'
#' @slot ptr an external pointer to C++ data
setClass("FromPOSet", contains = c("VIRTUAL", "POSet"))

#' An S4 class to represent a Product POSet.
#'
#' @slot ptr an external pointer to C++ data
setClass("ProductPOSet", contains = "FromPOSet")

#' An S4 class to represent a Lexicographic Product POSet.
#'
#' @slot ptr an external pointer to C++ data
setClass("LexicographicProductPOSet", contains = "FromPOSet")

#' An S4 class to represent a Binary Variable POSet.
#'
#' @slot ptr an external pointer to C++ data
setClass("BinaryVariablePOSet", contains = "POSet")

#' An S4 class to represent the exact linear extension generator.
#'
#' @slot ptr an external pointer to C++ data
setClass("LEGenerator", slots = list(ptr = "externalptr"))

#' An S4 class to represent the linear extension generator based on the Bubley-Dyer procedure.
#'
#' @slot ptr an external pointer to C++ data
#' @slot seed the 64-bit seed actually used, as a string of decimal digits
setClass("BubleyDyerGenerator", contains = "LEGenerator",
         slots = list(seed = "character"))

#' An S4 class to represent the exact MRP generator.
#'
#' @slot ptr an external pointer to C++ data
setClass("ExactMRPGenerator", slots = list(ptr = "externalptr"))

#' An S4 class to represent the MRP generator based on the Bubley-Dyer procedure.
#'
#' @slot ptr an external pointer to C++ data
#' @slot seed the 64-bit seed actually used, as a string of decimal digits
setClass("BubleyDyerMRPGenerator", slots = list(ptr = "externalptr", seed = "character"))

#' An S4 class to represent function evaluation based on the Bubley-Dyer procedure.
#'
#' @slot ptr an external pointer to a C++ data
#' @slot seed the 64-bit seed actually used, as a string of decimal digits
setClass("BubleyDyerEvaluationGenerator", slots = list(ptr = "externalptr", seed = "character"))

#' An S4 class to represent function separation based on the Bubley-Dyer procedure.
#'
#' @slot ptr an external pointer to a C++ data
#' @slot types list of separation to be computed
#' @slot seed the 64-bit seed actually used, as a string of decimal digits
setClass("BubleyDyerSeparationGenerator", slots = list(ptr = "externalptr", types="list", seed = "character"))
