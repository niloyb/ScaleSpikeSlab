#' Riboflavin GWAS dataset
#'
#' Dataset of riboflavin production by Bacillus subtilis containing n = 71 observations 
#' of a one-dimensional response (riboflavin production) and p = 4088 predictors (gene expressions).
#' The one-dimensional response corresponds to riboflavin production. 
#'
#' @docType data
#'
#' @usage data(riboflavin)
#'
#' @format A data frame containing a vector y of length 71 (responses) and 
#' a matrix X of dimension 71 by 4088 (gene expressions)
#'
#' @keywords datasets
#'
#' @references 
#' Buhlmann, P., Kalisch, M. and Meier, L. (2014) 
#' \emph{High-dimensional statistics with a view towards
#' applications in biology}. Annual Review of Statistics and
#' its Applications \bold{1}, 255--278
#'
#' @details 
#' The processed dataset is the same as in the R packages qut and hdi.
#' 
#' @examples
#' data(riboflavin)
#' y <- as.vector(riboflavin$y)
#' X <- as.matrix(riboflavin$x)
"riboflavin"