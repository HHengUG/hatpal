#' Extracting the chromosomes information from the indexed BAM file, including names and length.
#'
#' @param in_bamfile A indexed BAM file, and the index files (.bai) must be in the same directory of the indexed BAM
#'
#' @return A data.table with chromosomes names and length
#' @export
#'
#' @import data.table
#' @rawNamespace import(GenomicAlignments, except = c(first, second, last))
#'
#' @examples
#' chr_info <- ExtractChrinfo("input/test.sorted.bam")
#'
ExtractChrinfo <- function(in_bamfile){
  if (!file.exists(in_bamfile)){
    warning(paste0(in_bamfile, " does not exist"))
    return(message("ERROR: check the file"))
  }
  print("Scan bam to extract the info of Chrs")
  para <- Rsamtools::ScanBamParam(what = c("rname", "pos"))
  chr_get <- Rsamtools::scanBam(in_bamfile, param = para)
  chr_info <- data.table(
    chr = chr_get[[1]]$rname,
    pos = chr_get[[1]]$pos
  )
  chr_info <- chr_info[, max(pos), by = chr]
  chr_info <- chr_info[complete.cases(chr,V1),]  # remove na
  rm(chr_get)
  gc()
  print(paste0("There are ", nrow(chr_info), " Chrs"))
  return(chr_info)
}
