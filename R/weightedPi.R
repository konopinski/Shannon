#' @title weightedPi function \cr
#'
#' @keywords weightedPi
#' 
#' @description This function estimates nucleotide diversity \emph{sensu} Nei and Li 1979 
#'   from \emph{DNAbin} objects. \emph{Fasta} files can be easily converted to \emph{DNAbin} using
#'   \code{ape::read.FASTA()} function with \code{type = "DNA"}. Check \code{ape} documentation for details.
#'   Should also work with fastq files but I didn't test that yet.
#' 
#' @param sequences \emph{DNAbin} object.
#'
#' @return The function returns a value of average weighted nucleotide diversity for the sequences 
#'   for the samples included in \emph{DNAbin}  object.
#' 
#' @references Konopinski MK. (unpublished) Weighted nucleotide diversity is more precise 
#'   than pixy in estimating true value of \eqn{\pi} from sequences containing missing data. 
#'   \emph{???} . DOI: \cr
#' @references Nei M and Li WH. 1979. Mathematical model for studying genetic variation 
#'   in terms of restriction endonucleases. \emph{PNAS}, 76(10), 5269–5273. 
#'   DOI: 10.1073/pnas.76.10.5269
#' 
#' @export

weightedPi <- function(sequences){
  dnas <- as.data.frame(as.character(sequences))
  freqs <- apply(dnas,1,table,exclude = "n")
  nDiff <- sapply(freqs,function(z){ 
    tabs <- outer(z,z)
    sum(tabs[lower.tri(tabs)])})
  nComp <- sapply(freqs,sum)
  nComp <- (nComp*(nComp-1))/2
  return(mean(nDiff/nComp,na.rm = TRUE))
}

