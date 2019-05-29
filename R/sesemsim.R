#' Sesame Street Data
#'
#' This is a simulated counterpart of part of the Sesame street data presented by Stevens (1996, Appendix A) concerning
#' the effect of the first year of the Sesame street series on the knowledge of 240 children in the age range 
#' 34 to 69 months. We will use the following variables: Age in months, Bnum and Anum, the knowledge of
#' numbers before and after watching Sesame street for a year; analogously we will use Bbody/Abody, Blet/Alet, 
#' Bform/Aform, Brel/Arel, Bclas/Aclas, representing the knowledge of body parts, letters,  forms, relations, 
#' and classifications, respectively (score ranges vary between 1 to 20 and 1 to 70); Pea, the Peabody test 
#' measuring the mental age of children (score range 15 to 89); and Sex, where the score 1 denotes the boys 
#' and 2 the girls. Based on the original Sesame street data a data set with the same characteristics 
#' was simulated using a multivariate normal distribution for all the variables in both sex groups. 
#'
#' \tabular{lll}{
#'    \strong{age} \tab \code{numeric} \tab Participant age in months\cr
#'    \strong{pea} \tab \code{numeric} \tab Participant metal age score\cr
#'    \strong{Bbody} \tab \code{numeric} \tab Knowledge of body parts before\cr
#'    \strong{Blet} \tab \code{numeric} \tab Knowledge of letters before\cr
#'    \strong{Bform} \tab \code{numeric} \tab Knowledge of forms before \cr
#'    \strong{Bnum} \tab \code{numeric} \tab Knowledge of numbers before\cr
#'    \strong{Brel} \tab \code{numeric} \tab Knowledge of relations before\cr
#'    \strong{Bclass} \tab \code{numeric} \tab Knowledge of classifications before\cr
#'    \strong{Abody} \tab \code{numeric} \tab Knowledge of body parts after\cr
#'    \strong{Alet} \tab \code{numeric} \tab Knowledge of letters after\cr
#'    \strong{Aform} \tab \code{numeric} \tab Knowledge of forms after\cr
#'    \strong{ABnum} \tab \code{numeric} \tab Knowledge of numbers after\cr
#'    \strong{Arel} \tab \code{numeric} \tab Knowledge of relations after\cr
#'    \strong{Aclass} \tab \code{numeric} \tab Knowledge of classifications after\cr
#'    \strong{sex} \tab \code{factor} \tab sex\cr
#' }
#' @docType data
#' @keywords datasets
#' @name sesemsim
#' @usage data(sesemsim)
#' @references Stevens, J. (1996). Applied Multivariate Statistics for the Social Sciences. Mahwah NJ: Lawrence Erlbaum.
#' @format A data.frame with 240 rows and 15 columns.
NULL
