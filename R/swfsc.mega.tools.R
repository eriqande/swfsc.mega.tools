#' Streamlining SWFSC Molecular Genetics and Analysis Team workflows with R
#'
#' \code{swfsc.mega.tools} is a nascent effort to bring together bits of
#' R code that we have written for various applications into a single
#' R package that everyone in the lab should find easy to use.
#' 
#' Utlimately we hope that it will include functions that will assist
#' in a variety of analyses we do in the lab.  The several areas 
#' that the package addresses are detailed below.
#' 
#' \strong{Standardizing Files of Genotypes}
#' 
#' Our lab has a frightening mound of genetic data sets that were collected
#' over various time periods and by various people in the lab.  Some of
#' these data may remain in that state, but in the effort to extract
#' data from these files for loci that will continue to use going forward 
#' we've got some tools here.
#' \itemize{
#'  \item{Locus Sets:}{The package includes names for loci that have been
#'    stanardized and should be used moving forward.  These are stored as 
#'    simple character vectors.  Currently you will find:
#'    \code{\link{Omykiss.Standard.Loci}} }
#'  \item{Functions:}{The function \code{\link{CompareLocusNames}} lets you
#'    compare the locus headers in a file (read into a data frame) with another
#'    file or with a string of locus names, like \code{\link{Omykiss.Standard.Loci}}.
#'    Some small functions related to this that might be found to be helpful are: 
#'    \code{\link{ExcelColumns}}, and  \code{\link{FindLocusColumns}}. }
#' }
#' 
#' \strong{Datasets}
#' 
#' Currently we don't have any standardized data sets in this package, but we 
#' may include some (like our baseline files, etc).  Until that time however,
#' the package includes some small data sets, like \code{\link{sthd.geno.A}}
#' and \code{\link{sthd.geno.B}} for use in function examples.
#' 
#'
#' @import stringdist
#' @docType package
#' @name swfsc.mega.tools
NULL