#' Same location event statistic
#'
#' A matrix coding whether senders of events (in the rows) and receivers of events
#' (in the column) have the same location. Related to the 'events' data object,
#' that contains a relational event sequence, and the 'actors' object, that contains
#' information on the 25 actors involved in the relational event sequence.
#'
#'
#' @name same_location
#' @docType data
#' @usage data(same_location)
#' @format dataframe (25 rows, 4 columns)
#'
#' \tabular{lll}{
#'    \strong{same_location} \tab \code{integer} \tab Event statistic. Matrix with senders in the
#'    rows and receivers in the columns. The event statistic is 1 if sender and receiver have
#'    the same location and 0 otherwise. \cr
#' }
#'
#' @keywords datasets
NULL
