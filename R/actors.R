#' Actors from a small hypothetical network
#'
#' The related data files 'events', 'same_location', 'same_culture' contain
#' information on the  event sequence and the two event statistics respectively.
#'
#'
#' @name actors
#' @docType data
#' @usage data(actors)
#' @keywords datasets
#' @format dataframe (25 rows, 4 columns)
#'
#' \tabular{lll}{
#'    \strong{actors$id} \tab \code{integer} \tab ID of the employee, corresponding to
#'    the sender and receiver IDs in the events dataframe \cr
#'    \strong{actors$location} \tab \code{numeric} \tab Location of the actor,
#'    ranging from 1-4 \cr
#'    \strong{actors$culture} \tab \code{character} \tab Categorical variable, indicating the
#'    culture of the employee \cr
#' }
#'
#' @keywords datasets
NULL
