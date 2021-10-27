#' @title weightedPi function \cr (function to calculate nucleodie diversity from fasta file)
#'
#' @keywords weightedPi
#' 
#' @param fasta path to \emph{fasta} file.
#'
#' @return The function returns a value of average weighted nucleotide diversity for the sequences 
#'  from the provided fasta file
#' 
#' @references Konopinski MK. (unpublished) Weighted nucleotide diversity is more precise 
#'  than pixy in estimating true value of  from sequences containing missing data. 
#'  \emph{submitted} . doi: \cr
#' 
#' @export

weightedPi <- function(fasta){
  sequences <- ape::read.FASTA(fasta,type = "DNA")
  dnas <- as.data.frame(as.character(sequences))
  freqs <- apply(dnas,1,table,exclude = "n")
  nDiff <- sapply(freqs,function(z){ 
    tabs <- outer(z,z)
    sum(tabs[lower.tri(tabs)])})
  nComp <- sapply(freqs,sum)
  nComp <- (nComp*(nComp-1))/2
  return(mean(nDiff/nComp,na.rm = TRUE))
}

