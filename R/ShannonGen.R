#' @title ShannonGen function \cr (modified code from \code{SpadeR})
#'
#' @keywords ShannonGen
#' 
#' @param gInd a \emph{genind} object. For details on \emph{genind} objects see 
#'       \href{https://github.com/thibautjombart/adegenet}{\code{adegenet}} documentation.
#' @param estimator a character string, a vector or a list of \emph{H} estimators that
#'              will be estimated:
#'              \itemize{
#'                \item\emph{ch}{ - Chao et al 2013,}
#'                \item\emph{cs}{ - Chao and Shen 2003,}
#'                \item\emph{sh}{ - the original Shannon diversity index (Shannon 1948),}
#'                \item\emph{z}{ - Zahl 1977 (\strong{default}).}
#'                }
#' @details For the details on performance of Shannon \emph{H} estimators
#'  see Konopinski MK 2020 (\emph{in press})
#'  The functions \emph{MLE}, \emph{Z}, \emph{CS} and \emph{Ch} were copied from \code{Diversity} function
#'  from \href{https://github.com/AnneChao/SpadeR}{\code{SpadeR}} package
#'  (WARNING: the function was tested only on microsatellite data)
#'
#' @return The function returns a list of data frames (one list element per each estimator)
#'   containing selected Shannon \emph{H} estimators for all loci. Rows represent loci, while columns
#'   represent populations.
#'
#' @references Chao A., Ma K.H., Hsieh T.C. and Chiu C.H. 2016. SpadeR (Species-richness
#'      Prediction And Diversity Estimation in R): an R package in CRAN. Program and User's
#'      Guide also published at \href{https://github.com/AnneChao/SpadeR}{https://github.com/AnneChao/SpadeR}
#' @references Chao A. and Shen T.J. 2003. Nonparametric estimation of Shannon's index of
#'      diversity when there are unseen species. \emph{Environmental and Ecological Statistics}
#'      10:429-443. DOI: 10.1023/A:1026096204727
#' @references Chao A., Wang Y.T., Jost L. 2013. Entropy and thespecies accumulation curve:
#'      a nearly unbiased estimator ofentropy via discovery rates of new species.
#'      \emph{Methods in Ecology and Evolution} 4:1091-1100. DOI: 10.1111/2041-210X.12108
#' @references Konopinski MK. 2020.  Shannon diversity index: a call to replace the original Shannon’s 
#' 		formula with unbiased estimator in the population genetics studies. \emph{PeerJ, 8:e9391}
#' 		\href{https://doi.org/10.7717/peerj.9391}
#' @references Shannon C.E. 1948. A Mathematical Theory of Communication.
#'      \emph{The Bell System Technical Journal}, 27:379-423,623-656.
#' @references Zahl S. 1977. Jackknifing An Index of Diversity. \emph{Ecology} 58 (4): 907-913.
#'      DOI: 10.2307/1936227
#' @examples
#' data(Dzik) ## sample genotypes from wildboars
#' 
#' ## Default usage:
#' ShannonGen(Dzik) ## returns Zahl 1977 unbiased estimator
#'
#' ## Calculating all the estimators:
#' ShannonGen(Dzik, estimator = c("z", "sh", "cs", "ch"))
#' 
#' @export

ShannonGen <- function(gInd, estimator = NULL) {
  MLE = function(X) {
    X = X[X > 0]
    n = sum(X)
    - sum(X / n * log(X / n))
  }
  Z = function(X){
    X=X[X>0]
    Y=X[X>1]
    n=sum(X)
    -n*sum(X/n*log(X/n))-(n-1)/n*sum((n-X)*(-X/(n-1)*log(X/(n-1))) )-(n-1)/n*sum(-Y*(Y-1)/(n-1)*log((Y-1)/(n-1)))
  }
  CS = function(X) {
    x = X
    x = x[x > 0]
    n = sum(x)
    f1 = sum(x == 1)
    C_head = 1 - f1 / n
    a = -sum(C_head * (x / n) * log(C_head * (x / n)) / (1 - (1 - C_head *
          (x / n)) ^ n))
    a
  }
  Ch = function(X) {
    x = X
    x = x[x > 0]
    n = sum(x)
    UE <- sum(x / n * (digamma(n) - digamma(x)))
    f1 <- sum(x == 1)
    f2 <- sum(x == 2)
    if (f1 > 0)
    {
      A <-
        1 - ifelse(f2 > 0, (n - 1) * f1 / ((n - 1) * f1 + 2 * f2), (n - 1) * f1 /
                     ((n - 1) * f1 + 2))
      B = sum(x == 1) / n * (1 - A) ^ (-n + 1) * (-log(A) - sum(sapply(1:(n - 1), function(k) {
        1 / k * (1 - A) ^ k})))
    }
    if (f1 == 0) {
      B = 0
    }
    if (f1 == 1 & f2 == 0) {
      B = 0
    }
    UE + B
  }
  inputs <- lapply(adegenet::seppop(gInd),
                   function (i)
                     sapply(adegenet::seploc(i),
    function(j) colSums(j@tab,na.rm = TRUE)))
  out <- list()
  if (is.null(estimator))
    estimator <- "z"
  
  for (i in estimator){
    if (i %in% c("z","sh","cs","ch")) 
           print (paste(i,"Letter code OK")) else {
             print(paste("Wrong letter code:",i))
    stop("<estimator> option is not correct. 
          For details see help: ?ShannonGen")}}

  if ("sh" %in% estimator) {
    out[["Shannon_1949"]] <-
      as.data.frame(lapply(inputs, function (i)
        sapply(i, MLE)))
  }
  if ("z" %in% estimator) {
    out[["Zahl_1977"]] <-
      as.data.frame(lapply(inputs, function (i)
        sapply(i, Z)))
  }
  if ("cs" %in% estimator) {
    out[["Chao_Shen_2003"]] <-
      as.data.frame(lapply(inputs, function (i)
        sapply(i, CS)))
  }
  if ("ch" %in% estimator) {
    out[["Chao_et_al_2013"]] <-
      as.data.frame(lapply(inputs, function (i)
        sapply(i, Ch)))
  }
  return(out)
}
