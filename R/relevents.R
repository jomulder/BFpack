#' A sequence of innovation-related e-mail messages
#'
#' A time-ordered sequence of 247 communication messages between 25 actors.
#'
#' The related data files 'actors', 'same_location', 'same_culture' contain information
#' on the actors and three event statistics respectively.
#'
#'
#' @name relevents
#' @docType data
#' @usage data(relevents)
#' @format dataframe (247 rows, 3 columns)
#'
#' \tabular{lll}{
#'    \strong{relevents$time} \tab \code{numeric} \tab Time of the e-mail message,
#'    in seconds since onset of the observation \cr
#'    \strong{relevents$sender} \tab \code{integer} \tab ID of the sender, corresponding to
#'    the employee IDs in the actors dataframe \cr
#'    \strong{relevents$receiver} \tab \code{integer} \tab ID of the receiver \cr
#' }
#'
#' @keywords datasets
NULL
