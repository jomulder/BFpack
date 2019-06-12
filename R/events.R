#' A sequence of innovation-related e-mail messages
#'
#' A time-ordered sequence of e-mail messages between employees of a 
#' consultancy firm and information on the actors in the relational event 
#' sequence. The data is orignally analyzed by Mulder & Leenders (2019), to 
#' find drivers of innovation-related e-mail messages exchanged between  
#' employees of a large consultancy firm. Originally, the data consist of 2081
#' e-mail messages exchanged between 70 employees over the course o fa year. 
#' The current data is a sample of a simulated data set, based on estimates of 
#' the model parameters in Mulder & Leenders (2019). 
#' 
#' The related data files actors', 'same_building', 'same_division' and 
#' 'same_hierarchy' contain information on the actors and three event statistics 
#' respectively. 
#' 
#' 
#' @name events
#' @docType data
#' @usage data(events)
#' @format dataframe (227 rows, 3 columns)
#'
#' \tabular{lll}{
#'    \strong{events$time} \tab \code{numeric} \tab Time of the e-mail message, 
#'    in seconds since onset of the observation \cr
#'    \strong(events$sender) \tab \code{integer} \tab ID of the sender, corresponding to
#'    the employee IDs in the actors dataframe \cr
#'    \strong{events$receiver} \tab \code{integer} \tab ID of the receiver \cr
#' }
#' 
#' @references Mulder, J., & Leenders, R. T. (2019). Modeling the evolution of 
#' interaction behavior in social networks: A dynamic relational event approach 
#' for real-time analysis. Chaos, Solitons and Fractal Nonlinear, 119, 73-85,
#' https://doi.org/10.1016/j.chaos.2018.11.027
#' \href{https://doi.org/10.1016/j.chaos.2018.11.027}{
#' doi:10.1016/j.chaos.2018.11.027}
#' @source \href{https://doi.org/10.1016/j.chaos.2018.11.027}{
#' doi:10.1016/j.chaos.2018.11.027}
#' @keywords datasets
NULL
