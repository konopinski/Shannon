#' ShannonGen (package)
#' 
#' @docType package
#' @name ShannonGen.package
#' @description 
#' WARNING! The package was tested only on microsatellite data!\cr
#' It is my first package so, please, let me know if something went wrong - mkonop(at)wp.pl\cr
#' 
#' @details 
#' The package has only one function: \code{\link{ShannonGen}}. The function calculates 
#' four different Shannon Diversity Index estimators. For the details see funcion's help.\cr
#' The package uses functions \emph{MLE}, \emph{Z}, \emph{CS} and \emph{Ch} 
#'  copied from \code{Diversity} function from 
#'  \href{https://github.com/AnneChao/SpadeR}{\code{SpadeR}}
#'  and \code{seploc} and \code{seppop} functions from
#'  \href{https://github.com/thibautjombart/adegenet}{\code{ adegenet}} packages.
#' Data should be in \emph{genind} format.\cr
#' 
#' @references Chao A., Ma K.H., Hsieh T.C. and Chiu C.H. 2016. SpadeR (Species-richness
#' Prediction And Diversity Estimation in R): an R package in CRAN. Program and User's
#' Guide also published at 
#' \href{https://github.com/AnneChao/SpadeR}{https://github.com/AnneChao/SpadeR}
#' 
#' @references Jombart T. (2008) adegenet: a R package for the multivariate
#' analysis of genetic markers \emph{Bioinformatics} 24: 1403-1405. doi:
#' 10.1093/bioinformatics/btn129\cr
#' 
#' @import adegenet
#' @importFrom adegenet seppop
#' @importFrom adegenet seploc
#' 
NULL