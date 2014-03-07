#' do stringdist but handle NAs in a way that is appropriate to my task
#' 
#' It seems that stringdist returns NA when any strings being compared to
#' are NA.  I needed somethign else.
#' @param x character vector
#' @param y character vector of same length as x. 
#' @param ... other arguments to pass to the function \code{stringdist} 
#' @return returns the string distance between the i-th element in x and the i-th element in y
#'   and returns a vector with length=length(x) containing NAs wherever x or y was NA.
#' @export
#' @examples
#' stringdist.handle.NAs(c("Big", "Bad", "Bob", NA), c(NA, "Bed", "Bob", NA))
stringdist.handle.NAs <- function(x, y, ...) {
  if(length(x) != length(y)) stop("x and y must be the same length")
  ret <- rep(NA, length(x))
  keep <- !is.na(x) & !is.na(y)
  res <- stringdist(x[keep], y[keep], ...)
  ret[keep] <- res
  ret
}




#' Check for locus names that are not consistent between the two alleles of the same locus
#' 
#' @param D A data frame in which the loci start in column 1 and the column for the second allele
#'   at a locus named "LocusX" is "LocusX.1".  Thus, it is the sort of 
#'   data frame that would result from reading in a two-column (ToolKit) format file while setting
#'   the first column to be the row names, and using R's default method for making column names unique
#' @return A matrix of the two instances of each locus name that are not identical in the
#' file.   
#' @details Note that if someone already added ".1" to the second occurrence of each locus name, then
#' that will work fine as R will already find it to be unique (however adding ".a" or ".b" would break things).
#' @export
#' @examples
#'   silly.data <- data.frame(LocA=c(2,2), LocA=c(3,2), LocB=c(45,42), LocBB=c(48,42))
#'   rownames(silly.data) <- c("Fish1", "Fish2")
#'   silly.data  # just look at it. Notice that LocB and LocBB aren't right!
#'   CheckLocusNames(silly.data)
#'   
#'   # note that the example steelhead data are OK
#'   CheckLocusNames(sthd.geno.A)
CheckLocusNames <- function(D) {
  matrix(names(D)[which(!(names(D)[c(T,F)] == sub("\\.1$", "", names(D)[c(F,T)]))) * 2 + c(-1,0)], ncol=2)
}



#' give the column header (A or J or AC, etc) name in an excel file for 
#' 
#' @param C a vector of integers giving the column numbers
#' @param Shift the amount to add to C before returning the column name. This
#'   is useful when columns were stripped from the Excel file to get it into 
#'   an R data frame.  For example, if there were 15 columns in the Excel file
#'   that were dropped before reading the remaining columns, Shift would be 15
#' @details As currently implemented, this can give values out to 18278 columns
#' @export
#' @examples
#'   ExcelColumns(11:18)
#'   ExcelColumns(11:18, Shift=10)
ExcelColumns <- function(C, Shift=0) {
  c(LETTERS, paste(rep(LETTERS, each=26), LETTERS, sep=""),  paste(rep(LETTERS, each=26^2), paste(rep(LETTERS, each=26), LETTERS, sep=""), sep="") )[C+Shift]
}





#' Find the excel columns holding specific loci.
#' 
#'  Given a vector X of locus names for which you want positions, and a vector of locus names N that are in 
#'  an excel file somewhere this returns the column names a la letters
#'  @param X names of loci you want to try to find in the file
#'  @param N names of the loci in the order they appear in in the excel file
#'  @param Shift the number of the column in the excel file that the first locus in N appears
FindLocusColumns <- function(X, N, Shift=1) {
  tmp <- sapply(X, function(z) which(z==N))
  tmp2 <- sapply(tmp, function(x) {if(length(x)>0) { x } else NA})  # this has NAs where they should be
  tmp3 <- rep(tmp2, each=2) + c(0,1)
  tmp4 <- ExcelColumns(tmp3, Shift=Shift)
  tmp4[is.na(tmp4)] <- "-"
  apply(matrix(tmp4, ncol=2, byrow=T), 1, paste, collapse=",")
}


#' find closely matching locus names between two files
#' 
#' @param First Either a vector of locus names, or a data frame in which the loci start in column 1 and the column for the second allele
#'   at a locus named "LocusX" is "LocusX.1".  Thus, it is the sort of 
#'   data frame that would result from reading in a two-column (ToolKit) format file while setting
#'   the first column to be the row names, and using R's default method for making column names unique
#' @param Second a vector of locus names or a data frame formatted like \code{First}.
#' @param maxD  The maximum string distance between two locus names that are still considered an approximate match
#' @return A list with three components:
#'  \itemize{
#'  \item{\code{InFirstButNotSecond}: }{Loci appearing in First but with no \emph{exact} matches in Second}
#'  \item{\code{RemainingInSecondButNotInFirst}: }{Any loci in second that had neither approximate nor exact matches
#'     to any loci in First.}
#'	\item{\code{InBoth}: }{Loci that have an exact matching name in both First and Second.}
#' } 
#' Each component of the returned list has a matrix with 5 columns with headers as follows:
#' \itemize{
#' \item{\code{NameInFirst}: }{ Locus name as it appears in First  }
#' \item{\code{NameInSecond}: }{ Locus name as it appears in Second }
#' \item{\code{StringSimilarity}: }{ The distance between the two approximately matching names.  Lower means more similar. }
#' \item{\code{ColumnsInFirst}: }{ If First is a data.frame, these correspond to the columns in an excel file that hold the original data table. }
#' \item{\code{ColumnsInSecond}: }{ See above, but for second. }
#' }
#'  @export
#'  @examples
#'  # look for closely matching locus names between the two steelhead example data sets
#'  CompareLocusNames(sthd.geno.A, sthd.geno.B)  
#'
#'  # look for the standard Omykiss locus names in sthd.geno.B
#'  CompareLocusNames(sthd.geno.B, Omykiss.Standard.Loci)
CompareLocusNames <- function(First, Second, maxD=7) {
  
  
  ## Deal with whether a and b are data frames or character vectors
  if(is.character(First)) {
      na<-First;
  } else if(is.data.frame(First)) {
      na<-names(First)[c(T,F)];
  } else {
    stop("Argument First must be either a character vector or a data.frame");
  }
  
  
  if(is.character(Second)) {
    nb<-Second; 
  } else if(is.data.frame(Second)) {
    nb<-names(Second)[c(T,F)]; 
  } else {
    stop("Argument Second must be either a character vector or a data.frame");
  }
  
  if(any(duplicated(na))) stop("Duplicated locus names in First")
  if(any(duplicated(nb))) stop("Duplicated locus names in Second")
   
   
 	both <- intersect(na, nb)  # names common to both
  X <- setdiff(na, nb)  # names in a that are not in b
  Y <- setdiff(nb, na)  # locus names in b that are not in a
  
  YsLikeX <- Y[amatch(X,Y, maxDist=maxD)]  # closest approximate matches to X that are in Y
  XsLikeY <- X[amatch(Y,X, maxDist=maxD)]  # closest approximate matches to Y that are in X
  
  Y.extra <-  setdiff(Y, YsLikeX)  # remaining Y's
  XsLikeY.extra  <- X[amatch(Y.extra,X, maxDist=maxD)]  # closest approximate matches in X to Y.extra
  
  # now, prepare some things to ouput the results
  matnames <- c("NameInFirst", "NameInSecond", "StringSimilarity", "ColumnsInFirst", "ColumnsInSecond")

  
  # a function to create the output matrices. G can be X or Y or both, HsLikeG can be XsLikeY or YsLikeXs or both, 
  makeMats <- function(G, HsLikeG) {
  	if(length(G)==0) {
  		z<-matrix(NA, nrow=0, ncol=length(matnames)); 
  		colnames(z)<-matnames;
  	}
  	else {
  		Pos1<-Pos2<-NA
  		if(is.data.frame(First)) Pos1 <-  FindLocusColumns(G, names(First))
  		if(is.data.frame(Second)) Pos2 <-  FindLocusColumns(HsLikeG, names(Second))
  		z<-cbind(
  			G,
  			HsLikeG,
  			stringdist.handle.NAs(G, HsLikeG),
  			Pos1,
  			Pos2
  			)
  		colnames(z) <- matnames	
  	}
  	z
  } 

	# get the list we want to return:
	res <- list(
		InFirstButNotSecond = makeMats(X, YsLikeX),
		RemainingInSecondButNotInFirst = makeMats(XsLikeY.extra, Y.extra),
		InBoth = makeMats(both, both)
		)
		
		res
	
}

