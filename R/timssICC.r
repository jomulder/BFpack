#' Trends in International Mathematics and Science Study (TIMSS) 2011-2015
#'
#' A stratified sample was drawn by country and school to obtain a balanced
#' sample of p = 15 grade-4 students
#' per school for each of four countries (The Netherlands (NL), Croatia (HR),
#' Germany
#' (DE), and Denmark (DK)) and two measurement occasions (2011, 2015).
#' Achievement scores
#' (first plausible value) of overall mathematics were considered. Performances
#' of fourth
#' and eight graders from more than 50 participating countries around the world
#' can be found at (https://www.iea.nl/timss)
#' The TIMSS achievement scale is centered at 500 and the standard deviation is
#' equal to 100 scale score points.
#' The TIMSS data set has a three-level structure, where students are nested
#' within classrooms/schools, and
#' the classrooms/schools are nested within countries. Only one classroom was
#' sampled per school.
#' Changes in the mathematics achievement can be investigated by examining the
#' grouping of
#' students in schools across countries. Changes in country-specic intraclass
#' correlation coefficient
#' from 2011 to 2015, representing heterogeneity in mathematic achievements
#' within and between schools across years,
#' can be tested. When detecting a decrease in average performance together
#' with an increase
#' of the intraclass correlation, a subset of schools performed worse. For a
#' constant
#' intraclass correlation across years the drop in performance applied to the
#' entire population
#' of schools. For different countries, changes in the intraclass correlation
#' across years
#' can be tested concurrently to examine also differences across countries.
#'
#' \tabular{lll}{
#'    \strong{math} \tab \code{numeric} \tab math score child\cr
#'    \strong{groupNL11} \tab \code{numeric} \tab
#'    Indicator for child from NL in 2011\cr
#'    \strong{groupNL15} \tab \code{numeric} \tab
#'    Indicator for child from NL in 2015\cr
#'    \strong{groupHR11} \tab \code{numeric} \tab
#'    Indicator for child from HR in 2011\cr
#'    \strong{groupHR15} \tab \code{numeric} \tab
#'    Indicator for child from HR in 2015\cr
#'    \strong{groupDE11} \tab \code{numeric} \tab
#'    Indicator for child from DE in 2011\cr
#'    \strong{groupDE15} \tab \code{numeric} \tab
#'    Indicator for child from DE in 2015\cr
#'    \strong{groupDR11} \tab \code{numeric} \tab
#'    Indicator for child from DK in 2011\cr
#'    \strong{groupDR15} \tab \code{numeric} \tab
#'    Indicator for child from DK in 2015\cr
#'    \strong{gender} \tab \code{numeric} \tab Female=0,Male=1 \cr
#'    \strong{weight} \tab \code{numeric} \tab Child sampling weight \cr
#'    \strong{yeargender} \tab \code{numeric} \tab
#'    Interaction for occassion and gender \cr
#'    \strong{lln} \tab \code{numeric} \tab
#'    total number of children in school-class \cr
#'    \strong{groupschool} \tab \code{factor} \tab
#'    Nested indicator for school in country\cr
#'    \strong{schoolID} \tab \code{factor} \tab
#'    Unique indicator for school
#' }
#' @docType data
#' @keywords datasets
#' @name timssICC
#' @usage data(timssICC)
#' @references Mulder, J. & Fox, J.-P. (2019). Bayes factor testing of multiple
#' intraclass correlations. Bayesian Analysis. 14, 2, p. 521-552.
#' @format A data.frame with 16770 rows and 15 columns.
NULL

