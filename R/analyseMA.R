## analyseMA.R -- program by Jobst Landgrebe and Frank Bretz to analyse
## microarray experiment designs
# $Id: analyseMA.R,v 1.2 2003/10/02 16:11:37 jgentry Exp $

# INPUT:
# 	data    : G x N data matrix (rows: genes; columns: arrays)
# 	design  : N x (K+2) design matrix X (rows: observations; columns: parameters)
#	id      : ID vector for genes
#	cmat    : P x (K+2) contrasts matrix (rows: experimental questions; columns: parameters)
#	cinfo   : flag vector, indicating which contrasts define a single experimental question
#	padj    : method for multiplicity ajustment (none, bonferroni or fdr)
#	tol     : tolerance parameter for estimability check (numerical)
#
# OUTPUT: matrix with elements
# 	ID, estimates of the linear functions (where cinfo == 1), first degrees of freedom
# 	for the F-test (where cinfo > 1), statistics (t or F), associated raw and adjusted
# 	P-values, means square error and residual degrees of freedom	

require(MASS)



########################

analyseMA <- function( data, design, id, cmat, cinfo, padj=c("none","bonferroni","fdr"), tol=1e-06 ) {

	# general checks
	if (sum(cinfo)   != nrow(cmat))   {stop("number of contrasts in cinfo and cmat unequal")}
	if (ncol(design) != ncol(cmat))   {stop("number of parameters in design and cmat unequal")}
	if (ncol(data)   != nrow(design)) {stop("number of arrays in data and design unequal")}
	if (length(id)   != nrow(data))   {stop("lengths of id and data unequal")}
	if (nrow(data) < 2 && padj != "none") {padj <- "none"}
	if (!is.vector(id))		  {stop("id not a vector")}
	if (!is.numeric(id))		  {stop("id must be numeric")}

    	len  <- length(cinfo)
	
	keep <- apply(data,1,function(x) sum(!is.na(x)) ) 
	data <- data[keep>1,]
	id   <- id[keep>1]
	cat("Deleting ",length(keep[keep<=1])," genes with all but one or all measurements missing \n" )

	res <- t(apply(data, 1, core, design, cmat, cinfo, tol))

	result     <- matrix(NA, nrow(data), 3*len + 3)
	result[,1] <- id
	result[,2:(3*len+3)] <- res

	names.a <- vector(mode="character", length=len)
	names.b <- vector(mode="character", length=len)

	for (i in 1:len){
		if(cinfo[i] == 1) { 
			names.a[i] <- paste("cont.est",i, sep = "");
			names.b[i] <- paste("T.stat",i, sep = "");
		}
		else {
			names.a[i] <- paste("df.num",i,sep="");
			names.b[i] <- paste("F.stat",i,sep="");
		}
	}
	
	if (padj != "none"){
                result <- cbind(result, apply(as.data.frame(res[,(2*len+1):(3*len)]), 2, p.adjust, padj))
                colnames(result) <- c("ID",names.a,names.b,paste("raw.p",1:len,sep=""),
                                    "MSE","df.residual",paste(paste("p.adj",".",sep=""),padj,1:len,sep=""))
        } else {
                colnames(result) <- c("ID",names.a,names.b,paste("raw.p",1:len,sep=""),
                                    "MSE","df.residual")
        }
        return(result)
}

### core function
core <- function( vector, design, cmat, cinfo, tol){

	len <- length(cinfo)
	out <- as.vector(rep(NA, 3*len+2))

	# drop NA arrays and extract data
	x <- design[!is.na(vector),]
	z <- as.matrix(vector[!is.na(vector)])

       	# compute variance estimate
       	n       <- nrow(x)
       	xt      <- t(x)
	xtx     <- xt %*% x
	gxx     <- ginv(xtx)
	gxxxtx  <- gxx %*% xtx
	f       <- n - round(sum(diag(gxxxtx)))
       	sse     <- t(z) %*% (diag(n) - x %*% gxx %*% xt) %*% z
       	mse     <- sse/f

	estim   <- abs(cmat %*% gxxxtx - cmat)
	count   <- 0

	# loop over contrasts
	for (i in 1:len) {  
		# check for estimability of every contrast
		if ( max( estim[ (sum(cinfo[1:i-1])+1) : sum(cinfo[1:i]), ] ) < tol ) { 
			
			# read out the i-th contrast matrix or vector, i.e., the i-th experimental question
		 	cmati  <- cmat[ (sum(cinfo[1:i-1])+1) : sum(cinfo[1:i]), , drop = FALSE ]
		 	cxxc   <- cmati %*% gxx %*% t(cmati)
                	
			# compute T-matrix
               		tmat   <- cmati %*% gxx %*% xt
                	
			# compute test statistic and raw P-value
               		if (cinfo[i]==1) { 
				estpar     <- tmat %*% z
               			stat       <- estpar / sqrt(mse * cmati %*% gxx %*% t(cmati) ) 
               			pval       <- 2*(1-pt(abs(stat), f))
				count      <- count + 1
				out[count] <- estpar 
               		}
               		else {
              			f2         <- round(sum(diag(ginv(t(tmat) %*% tmat) %*% (t(tmat) %*% tmat))))
               			stat       <- ( t(tmat %*% z) %*% ginv(cxxc) %*% (tmat %*% z) ) / (f2*mse)
               			pval       <- 1-pf(stat, f2, f)
				count      <- count + 1
				out[count] <- f2
               		}
       		
       		# contrast dependent output
       		out[  len+i] <- stat
       		out[2*len+i] <- pval
		} else {count <- count+1}
       	}

	# contrast independent output
       	out[3*len+1] <- mse
       	out[3*len+2] <- f

	return(out)
}

