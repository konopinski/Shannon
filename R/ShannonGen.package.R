#' ShannonGen (package)
#' 
#' @docType package
#' @name ShannonGen.package
#' @description 
#' WARNING! The package was tested only on microsatellite data!\cr
#' It is my first package so, please, let me know if something went wrong - mkonop(at)wp.pl\cr
#' 
#' @details 
#' The package has two functions: \code{\link{ShannonGen}} and \code{\link{weightedPi}}. 
#' \emph{ShannonGen} function calculates four different Shannon Diversity Index estimators. 
#' The function uses functions \emph{MLE}, \emph{Z}, \emph{CS} and \emph{Ch} 
#'  copied from \code{Diversity} function from 
#'  \href{https://github.com/AnneChao/SpadeR}{\code{SpadeR}}
#'  and \code{seploc} and \code{seppop} functions from
#'  \href{https://github.com/thibautjombart/adegenet}{\code{ adegenet}} packages.
#' Data should be in \emph{genind} format.\cr
#' \emph{weightedPi} function estimates average weighted nucleotide diversity.
#' The function imports data from fasta file using function
#' For the details see funcions' help.\cr
#' 
#' @author Maciej K. Konopinski 
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
#' @references Konopinski MK. (2020) Shannon diversity index: a call to replace the original 
#' Shannon’s formula with unbiased estimator in the population genetics studies 
#' \emph{PeerJ} 8:e9391. doi: 10.7717/peerj.9391\cr
#' 
#' @references Konopinski MK. (unpublished) Weighted nucleotide diversity is more precise 
#' than pixy in estimating true value of  from sequences containing missing data. 
#' \emph{submitted} . doi: \cr
#' 
#' @import adegenet
#' @importFrom adegenet seppop
#' @importFrom adegenet seploc
#' @import ape
#' @importFrom ape read.FASTA
#' 
NULL