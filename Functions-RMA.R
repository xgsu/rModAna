# rm(list=ls(all=TRUE))

# install.packages("sandwich")
# help(package="sandwich")
library(sandwich); 
# install.packages("stringr")
library(stringr)
# install.packages("gee")
# library(gee); 
# install.packages("geepack")
library(geepack)
# install.packages("lme4")
# library(lme4)
library(tidyverse)


# ===============================
# FUNCTION rdat() GENERATES DATA
# ===============================

# install.packages("MultiRNG")
require(MASS) # TO USE FUNCTION mvrnorm()
library(MultiRNG)

is.even <- function(x) x %% 2 == 0
expit <- function(x) (tanh(x/2)+1)/2  	# CORRECT & VERIFIED
# expit0 <- function(x) exp(x)/(1+exp(x))

rdat <- function(n=100, b0=c(-0.5, 2, 2, -2, 0, -2, 0, 0), link.function="logistic",
	rho=0.2, # UNIFORM COVARIATES
	observational=FALSE, trt.p=0.5, a0=c(-2, 2, 2, 0), details=FALSE)   # TRT ASSIGNMENT
{
	# GENERATE X
	if (!is.even(length(b0))) stop("Length of b0 should be an even number!")
	p <- (length(b0)-2)/2  # NUMBER OF COVARIATES
	S <- matrix(1, p, p)
	for (i in 1:p){
		for (j in 1:p){
			S[i, j] <- rho^(abs(i-j))
		}
	}	
	X <- draw.d.variate.uniform(no.row=n,d=p,cov.mat=S) 	

	# TREATMENT (RANDOM)
	if (observational) {
		if (length(a0)!=(p+1)) stop("Lenght of a0 must be consistent with b0.")
		eta.trt <- as.vector(cbind(1, X)%*%a0); 
		pi.trt <- expit(eta.trt)
		trt <- rbinom(n=n, size=1, prob=pi.trt)
		if (details) print(cbind(p.treatment=mean(pi.trt), p.trt=mean(trt)))
	} else {trt <- rbinom(n=n, size=1, prob=trt.p)}
	X0 <- cbind(1, trt, X, trt*X)

	# MODEL
	eta <- as.vector(X0%*%b0); 
	if (details) print(cbind(eta.mean=mean(eta), eta.min=max(eta), prop.positive=sum(eta>0)/length(eta)))
	if (link.function=="exp") pi0 <- exp(eta*(eta<=0))  # LOGLINEAR 
	else pi0 <- exp(eta)/(1+exp(eta))  # LOGISTIC 
	y <- rbinom(n, size=1, prob=as.vector(pi0)); 
	if (details) print(cbind(p.respone=mean(pi0), p.y=mean(y)))
		
 	dat <- data.frame(cbind(y, trt, X))
	colnames(dat) <- c("y", "trt", paste("x", 1:NCOL(X), sep=""))
	dat
}

# -----------------------------------
# FUNCTION swap() SWAPS TWO COLUMNS 
# -----------------------------------
swap <- function(DF, n, m)
{
  n <- if (class(n)=="character" & is.na(suppressWarnings(as.integer(n)))) which(colnames(DF)==n) else as.integer(n)
  m <- if (class(m)=="character" & is.na(suppressWarnings(as.integer(m)))) which(colnames(DF)==m) else as.integer(m)

  if (!(1<=n & n<=length(DF))) stop( "`n` represents invalid index!" )
  if (!(1<=m & m<=length(DF))) stop( "`m` represents invalid index!" )

  return (DF[ if (n==m) 1:length(DF) else c( (if (min(n,m)==1) c() else 1:(min(n,m)-1) ), (if (min(n,m)+1 == max(n,m)) (min(n,m)+1):(max(n,m)-1) else c( max(n,m), (min(n,m)+1):(max(n,m)-1), min(n,m))), (if (max(n,m)==length(DF)) c() else (max(n,m)+1):length(DF) ) ) ])
}

arrange.data <- function(dat, y="y", trt="trt") {
	n <- NROW(dat)
	dat1 <- swap(dat, y, trt)
	names(dat1) <- names(dat)
	dat <- data.frame(rbind(dat, dat1))
	dat$group <- rep(c(0,1), c(n, n))
	dat$ID <- rep(1:n, 2)
	return(dat)
}

as.numeric.factor <- function(x){as.numeric(levels(x))[x]}




# ====================================================
# REFINED MODERATE ANALYSIS - FUNCTION rModAna()
# ====================================================


rModAna <- function(formula, data, robust.VCOV=FALSE, adjust=TRUE, detail=FALSE)
{
  # ONE EXAMPLE OF FORMULA y~trt*(x1+x2+x3)+x4
  form <- formula
  vars <- all.vars(form)
  y.var <- vars[1]
  trt.var <- vars[2]
  form.direct <- form
  form.inverse <- as.formula(paste(trt.var, "~", sub(trt.var, y.var, form)[3], sep=" ")) 
  control0.glm <- glm.control(epsilon = 1e-8, maxit=50, trace = FALSE)

  # DIRECT MODEL
  # --------------
  fit.direct <- glm(form.direct, data=data, family=binomial(link = "logit"), 
                    control=control0.glm, x=TRUE)
  if (detail) {
    cat("The direct modeling fitting: \n")
    cat("=================================================================\n")
    print(summary(fit.direct))
  }
  terms.int.trt <- grep(pattern=paste(trt.var, ":", sep=""), names(coef(fit.direct)))  # TERMS FOR TRT BY COVARIATES INTERACTIONS
  vars.int <- str_remove(string=names(coef(fit.direct))[terms.int.trt], pattern=paste(trt.var, ":", sep=""))
  Beta.direct <- coef(fit.direct); beta.direct <- Beta.direct[terms.int.trt]
  P <- length(Beta.direct);  p <- length(beta.direct); n <- NROW(data) 
  # se.direct <- summary(fit.direct)$"coefficients"[terms.int.trt,2]
  # pvalue.direct <- summary(fit.direct)$"coefficients"[terms.int.trt,4]
  x.direct <- fit.direct$x 			# DESIGN MATRIX FOR DIRECT MODEL
  x.d2 <- x.direct[, terms.int.trt] 	# DESIGN MATRIX WITH INTERACTION TERMS ONLY  - DIRECT
  x.moderators <- x.direct[,vars.int] 		# DESIGN MATRIX OF MODERATORS
  pi.d <- fit.direct$fitted
  cov.d <- summary(fit.direct)$cov.scaled
  v.d <- diag(cov.d)[terms.int.trt];  # VARIANCE FOR INDIVIDUAL ESTIMATES (DIRECT)
  out.direct <- summary(fit.direct)$"coefficients"
  
  # INVERSE MODEL
  # ----------------
  fit.inverse <- glm(form.inverse, data=data, family=binomial(link = "logit"), 
                     control =control0.glm, x=TRUE)
  if (detail) {
    cat("The inverse modeling fitting: \n")
    cat("=================================================================\n")
    print(summary(fit.inverse))
  }
  terms.int.y <- grep(pattern=paste(y.var, ":", sep=""), names(coef(fit.inverse)))  # TERMS FOR RESPONSE BY COVARIATES INTERACTIONS
  Beta.inverse <- coef(fit.inverse); beta.inverse <- Beta.inverse[terms.int.y]
  # se.inverse <- summary(fit.inverse)$"coefficients"[terms.int.y,2]
  # pvalue.inverse <- summary(fit.inverse)$"coefficients"[terms.int.y,4]
  x.inverse <- fit.inverse$x # DESIGN MATRIX FOR INVERSE MODEL
  x.i2 <- x.inverse[, terms.int.y] # DESIGN MATRIX WITH INTERACTION TERMS ONLY - INVERSE
  pi.i <- fit.inverse$fitted
  cov.i <- summary(fit.inverse)$cov.scaled
  v.i <- diag(cov.i)[terms.int.y];  # VARIANCE FOR INDIVIDUAL ESTIMATES (INVERSE)
  out.inverse <- summary(fit.inverse)$"coefficients"

  # USING ROBUST SANDWICH ESTIMATORS OF VARIANCES?
  if (robust.VCOV){                                    
    v.d <- diag(sandwich(fit.direct))[terms.int.trt];
    v.i <- diag(sandwich(fit.inverse))[terms.int.y];
  }
  
  # COVARIANCE MATRIX VIA IEE (INDEPENDENT ESTIMATING EQUATION)
  # ------------------------------------------------------------
  # PREPARE THE MEAT MATRIX
  y <- data[, y.var]; trt <- data[, trt.var]
  r.y <- y-pi.d; r.t <- trt-pi.i
  X.d <- t(x.direct)%*%diag(r.y); X.i <- t(x.inverse)%*%diag(r.t)
  meat <- rbind(cbind(X.d%*%t(X.d), X.d%*%t(X.i)), 
		cbind(X.i%*%t(X.d), X.i%*%t(X.i)))
  if (adjust) meat <- n*meat/(n-(2*P-p))

  # PREPARE BREAD MATRIX
  W.d <- diag(pi.d*(1-pi.d)); W.i <- diag(pi.i*(1-pi.i))
  X <- x.direct[, 3:(2+p)]; 
  x.d1 <- x.direct[, 1:(2+p)]; x.i1 <- x.inverse[, 1:(2+p)]  # DESIGN MATRIX FOR ADDITIVE PARTS
  D.d <- diag(as.vector(Beta.direct[2] + X%*% Beta.direct[3:(2+p)])); 
  D.i <- diag(as.vector(Beta.inverse[2] + X%*% Beta.inverse[3:(2+p)])); 
  x.d2 <- x.moderators*pi.d; x.i2 <- x.moderators*pi.i;   # DESING MATRIX FOR INTERACTION PARTS		############### USE pi      ###############  
  # x.d2 <- x.moderators*y; x.i2 <- x.moderators*trt;   		                                      ############### USE y AND trt #######################  
  B12 <- t(x.direct)%*%W.d%*%D.d%*%W.i%*%x.d1
  B21 <- t(x.inverse)%*%W.i%*%D.i%*%W.d%*%x.i1
  bread <- rbind(cbind(t(x.direct)%*%W.d%*%x.direct,
                       B12, t(x.direct)%*%W.d%*%x.i2),
                 cbind(B21, t(x.inverse)%*%W.i%*%x.d2, 
                       t(x.inverse)%*%W.i%*%x.inverse)) 
  # SANDWICH ESTIAMTOR OF VCOV(BETA, ALPHA)
  V0.iee <-  solve(bread)%*%meat%*%solve(t(bread))    ### TRANSPOSED BREA ####  
  if (anyNA(V0.iee)) stop("You got a problem with the variance-coviance matrix!!")
  cov.d0 <- V0.iee[1:P, 1:P]
  cov.i0 <- V0.iee[(P+1):(2*P), (P+1):(2*P)] 
  cov.iee <- V0.iee[terms.int.trt, P+terms.int.y]
  v1 <- diag(cov.d0)[terms.int.trt]; v2 <- diag(cov.i0)[terms.int.y];
  # COMPUTE CORRELATION
  if (NROW(cov.iee) <=1) {
    rho <- cov.iee/sqrt(v1 *v2)  			
  } else rho <- diag(cov.iee)/sqrt(v1 *v2)  	### WATCH OUT: COULD rho BE OUTSIDE [-1, 1]?
  if (detail) {print(rho)} 
  if (sum(abs(rho)>1) >=1) {
	  print(rho); 
	  rho[abs(rho)>1] <- 0.999*sign(rho[abs(rho)>1])
	  warning("Hey, you got a correlation out of range!! It has been changed to +/- 0.999.") 
  }  
  
  # UNIVARIATELY COMBINED ESTIMATOR
  # ----------------------------------
  lambda <- v.d/v.i
  alpha <- (1-rho*lambda)/(1-2*rho*lambda+lambda^2)
  beta.unicom <- alpha*beta.direct + (1-alpha)*beta.inverse
  V.beta.unicom <- v.d*(1-rho^2)/(1-2*rho*lambda+lambda^2)
  SE.beta.unicom <- sqrt(V.beta.unicom)
  z.beta.unicom <- beta.unicom/SE.beta.unicom
  pvalue.beta.unicom <- pchisq(z.beta.unicom^2, df=1, lower.tail =FALSE)
  out.combined <- as.data.frame(cbind(beta.unicom, SE.beta.unicom, 
                                          z.beta.unicom, pvalue.beta.unicom))
  names(out.combined) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")  
  return(list(out.direct=out.direct, out.inverse=out.inverse, out.combined=out.combined, 
		V=V0.iee))
}
  






















#
