## designMA.R -- program by Frank Bretz and Jobst Landgrebe to compare
## microarray experiment designs


# $Id: designMA.R,v 1.3 2003/11/07 16:19:15 jgentry Exp $


# INPUT:
# 	design: design matrices (list)
#	cmat  : contrasts in each row (matrix of dimension 2)
#	cinfo : flag, indicating which contrasts define
#			a single experimental question (vector)
#	type  : optimality criterion (character)
#	tol   : tolerance parameter for estimability check (numerical)
#
# OUTPUT: list with elements
#  	alleff    : individual efficiencies (matrix of dimension 2)
#	alleffrel : relative efficiencie (matrix of dimension 2)
#	alleffave : average efficiencies (matrix of dimension 1)
#	effdesign : most efficient design (character)
#	df        : degrees of freedom (matrix of dimension 1)


########################

#predefine some objects needed

# 1. build up list of design matrices


#Basic Designs

designs.basic <- structure( list (

matrix( c(
1, -1,  1, -1,  0,  0,  0,  0,  0,
1, -1, -1,  1,  0,  0,  0,  0,  0,
1, -1,  0,  0,  1, -1,  0,  0,  0,
1, -1,  0,  0, -1,  1,  0,  0,  0,
1, -1,  0,  0,  0,  0,  1, -1,  0,
1, -1,  0,  0,  0,  0, -1,  1,  0 ), nr = 6, byrow = TRUE),

matrix( c(
1, -1,  1,  0, -1,  0,  0,  0,  0,
1, -1,  0,  0,  1,  0, -1,  0,  0,
1, -1, -1,  0,  0,  0,  1,  0,  0,
1, -1,  0,  1,  0, -1,  0,  0,  0,
1, -1,  0,  0,  0,  1,  0, -1,  0,
1, -1,  0, -1,  0,  0,  0,  1,  0 ), nr = 6, byrow = TRUE),

matrix( c(
1, -1,  1,  0,  0, -1,  0,  0,  0,
1, -1,  0,  0,  0,  1, -1,  0,  0,
1, -1,  0, -1,  0,  0,  1,  0,  0,
1, -1,  0,  1, -1,  0,  0,  0,  0,
1, -1,  0,  0,  1,  0,  0, -1,  0,
1, -1, -1,  0,  0,  0,  0,  1,  0 ), nr = 6, byrow = TRUE),

matrix( c(
1, -1,  1,  0, -1,  0,  0,  0,  0,
1, -1,  0,  0,  1,  0, -1,  0,  0,
1, -1,  0,  0,  0,  0,  1, -1,  0,
1, -1,  0,  0,  0, -1,  0,  1,  0,
1, -1,  0, -1,  0,  1,  0,  0,  0,
1, -1, -1,  1,  0,  0,  0,  0,  0 ), nr = 6, byrow = TRUE),

matrix( c(
1, -1,  1,  0,  0,  0,  0, -1,  0,
1, -1, -1,  0,  0,  0,  0,  1,  0,
1, -1,  0,  1,  0,  0, -1,  0,  0,
1, -1,  0, -1,  0,  0,  1,  0,  0,
1, -1,  0,  0,  1, -1,  0,  0,  0,
1, -1,  0,  0, -1,  1,  0,  0,  0 ), nr = 6, byrow = TRUE),

matrix( c(
1, -1,  1,  0,  0,  0, -1,  0,  0,
1, -1,  0,  0,  0, -1,  1,  0,  0,
1, -1, -1,  0,  0,  1,  0,  0,  0,
1, -1,  0,  1, -1,  0,  0,  0,  0,
1, -1,  0,  0,  1,  0,  0, -1,  0,
1, -1,  0, -1,  0,  0,  0,  1,  0), nr = 6, byrow = TRUE),

matrix( c(
1, -1,  1,  0,  0,  0,  0,  0, -1,
1, -1,  0,  1,  0,  0,  0,  0, -1,
1, -1,  0,  0,  1,  0,  0,  0, -1,
1, -1,  0,  0,  0,  1,  0,  0, -1,
1, -1,  0,  0,  0,  0,  1,  0, -1,
1, -1,  0,  0,  0,  0,  0,  1, -1 ), nr = 6, byrow = TRUE) ) ,

names = c("BS", "AL","XL","CL","RS","TL","CR")  )

#Composite Designs
attach(designs.basic)

designs.composite <- structure( list (

matrix(data=kronecker(BS,rbind(1,1,1)), nrow=18,ncol=9,
dimnames=list(paste("array_", 1:18, sep = ""),
 c("G", "R", paste("C", 1:7, sep = "") ) ) ) ,


matrix(data=rbind(BS,BS,AL), nrow=18,ncol=9,
dimnames=list(paste("array_", 1:18, sep = ""),
 c("G", "R", paste("C", 1:7, sep = "") ) ) ) ,

matrix(data=kronecker(AL,rbind(1,1,1)),  nrow=18,ncol=9,
dimnames=list(paste("array_", 1:18, sep = ""),
 c("G", "R", paste("C", 1:7, sep = "") ) ) ) ,


matrix(data=rbind(XL,AL,AL),  nrow=18,ncol=9,
dimnames=list(paste("array_", 1:18, sep = ""),
 c("G", "R", paste("C", 1:7, sep = "") ) ) ) ,

matrix(data=rbind(XL,BS,BS), nrow=18,ncol=9,
dimnames=list(paste("array_", 1:18, sep = ""),
 c("G", "R", paste("C", 1:7, sep = "") ) ) ) ,


matrix(data=rbind(CL,CL,TL), nrow=18,ncol=9,
dimnames=list(paste("array_", 1:18, sep = ""),
 c("G", "R", paste("C", 1:7, sep = "") ) ) ) ,


matrix(data= rbind(CL,CL,RS), nrow=18,ncol=9,
dimnames=list(paste("array_", 1:18, sep = ""),
 c("G", "R", paste("C", 1:7, sep = "") ) ) ) ,


matrix(data=kronecker(XL,rbind(1,1,1)),  nrow=18,ncol=9,
dimnames=list(paste("array_", 1:18, sep = ""),
 c("G", "R", paste("C", 1:7, sep = "") ) ) ) ,


matrix(data=rbind(CL,CL,XL), nrow=18,ncol=9,
dimnames=list(paste("array_", 1:18, sep = ""),
 c("G", "R", paste("C", 1:7, sep = "") ) ) ) ,


matrix(data=kronecker(CR,rbind(1,1,1)), nrow=18,ncol=9,
dimnames=list(paste("array_", 1:18, sep = ""),
 c("G", "R", paste("C", 1:7, sep = "") ) ) ) ),

names = c("BSBSBS", "BSBSAL", "ALALAL",  "XLALAL", "XLBSBS",
"CLCLTL","CLCLRS","XLXLXL","CLCLXL","CRCRCR") )

detach(designs.basic)


# 2. build up contrast matrices and info vectors

cmat <- matrix( data = rbind( 	cbind( kronecker( diag(rep(1,4)), cbind(1,-1) ), rep(0,4) ),
				# dye effect and simple B effect
				cbind( rep(0,3), rep(0,3), kronecker((diag(rep(1,3)) - matrix(1,3,3)/3), cbind(1,1)/2), rep(0,3)),
				# A effect
				cbind( 0, 0, kronecker(cbind(1,1,1)/3, cbind(1,-1)), 0 ),
				# B effect
				cbind( rep(0,3), rep(0,3), kronecker((diag(rep(1,3)) - matrix(1,3,3)/3) , cbind(1,-1)), rep(0,3))
		 		# AB interaction effect
		       ),
		nrow=11,ncol=9 ,
		dimnames=list(	c( "R/G", "11-12", "21-22", "31-32", "A", "A", "A", "B", "AB", "AB", "AB"),
				c("G", "R", paste("C", 1:7, sep = "") )
			 )
		)


cmatB.AB <- cmat[8:11,]
cinfo <- c(1, 1, 1, 1, 3, 1, 3)
cinfoB.AB <- c(1,3)


#####################################

designMA <- function( design.list, cmat, cinfo, type=c("d", "e", "t"), tol=1e-06 ) {

	len    <- length(design.list) # number of designs
	alleff <- matrix(nrow = length(cinfo), ncol = len)
	df     <- matrix(nrow = 1, ncol = len)

	# general checks
	if (sum(cinfo) != nrow(cmat)) {stop("number of contrasts incorrect")}
	for (i in 1:len) {
		if (ncol(design.list[[i]]) != ncol(cmat)) {
		 stop("number of columns unequal in design and cmat")
		}
		if (ncol(design.list[[i]]) != ncol(design.list[[1]]) ){
		 stop("design matrices have unequal numer of rows")
		}
	}

	#loop over design matrices
	for (k in 1:len) {
		x <- design.list[[k]]

		xtx     <- t(x) %*% x
		gxx     <- ginv(xtx)
		gxxxtx  <- gxx %*% xtx
		estim   <- abs(cmat %*% gxxxtx - cmat)

		# loop over contrasts
		for (i in 1:length(cinfo)) {
			# check for estimability of every contrast
			if ( max( estim[ (sum(cinfo[1:i-1])+1) : sum(cinfo[1:i]), ] ) > tol ) { eff <- NA }
			else {
				# read out the i-th contrast matrix or vector, i.e. the i-th experimental question
			 	cmati <- cmat[ (sum(cinfo[1:i-1])+1) : sum(cinfo[1:i]), , drop = FALSE ]
				cmati <- t(cmati) %*% ginv(cmati %*% t(cmati)) %*% cmati
			 	cxxc  <- cmati %*% gxx %*% t(cmati)

			 	# choice of optimality criteria
			 	switch( type,
					t = {
						s      <- round(sum(diag(ginv(cxxc) %*% cxxc)))
						eff    <- s/sum(diag(cxxc))
					},
					d = {
						eigval <- eigen(cxxc)$values
			     			eigval <- sapply(eigval, function(x) ( if(x<tol) {x <- 1} else {x} ) )
						s      <- round(sum(diag(ginv(cxxc) %*% cxxc)))
						eff    <- 1/prod(eigval)**(1/s)
				 	},
					e = {eff <- max( eigen(ginv(cxxc))$values )}
			 	)
			}
			alleff[i,k] <- eff
		}
		df[k] <- nrow(x)-round(sum(diag(gxxxtx)))
	}

	rownames(alleff) <- paste("C", 1:nrow(alleff), sep = "")
	name <- names(design.list)
	if (!is.null(name)) {
		colnames(alleff) <- name
		colnames(df)     <- name
	}
	alleffrel <- alleff/apply(alleff,1,max,na.rm=TRUE)
        alleffave <- colMeans(alleffrel)
	effdesign <- names(alleffave[alleffave == max(alleffave,na.rm=TRUE)])
	effdesign <- effdesign[!is.na(effdesign)]

	effres    <- list(alleff=alleff,alleffrel=alleffrel,alleffave=alleffave,effdesign=effdesign,df=df)
}

