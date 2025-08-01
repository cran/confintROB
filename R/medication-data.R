#' The medication data set
#'
#' Often used for didactic
#' purposes (Singer & Willett, 2003), and originally discussed in Tomarken, Shelton, Elkins, and Anderson (1997).
#' During seven days, three times a day (from time = 0 to time = 6.667), a sample of n = 64
#' adults (identified by the variable id) were randomly assigned to either a treatment group (treat=1)
#' or a control group (treat=0) and were required to report their mood (pos).
#' @format `medication`
#' a data.frame with 5 columns and 1242 rows:
#' \describe{
#' \item{obs}{row number}
#' \item{id}{participant number}
#' \item{treat}{treatment assignment, 1= treatment; 0= control}
#' \item{time}{time from 0 to 6.667, with increments of 0.333}
#' \item{pos}{the positive mood score}}
#' @docType data
#'
#' @usage data(medication)
#'
#' @keywords datasets
#'
#' @references Tomarken, A. J., Shelton, R. C., Elkins, L., & Anderson, T. (1997). Sleep deprivation and anti-depressant medication: unique effects on positive and negative affect. In American Psychological Society Meeting, Washington, DC.
"medication"
