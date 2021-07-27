
# ===============================
# SIMULATION STUDY
# ===============================

rm(list=ls(all=TRUE))
source("Functions-RMA.R")

# WEAK
b0 <- c(-1, 0.5, 0.5, -1, 1, -0.5, 1, 0); a0 <- c(-0.5, 1, -0.5, 1)
# STRONG
# b0 <- c(-1.5, 1.5, 1.5, -2, 2, -1.5, 2, 0); a0 <- c(-1, 1.5, -1, 1.5)
# dat <- rdat(n=500, b0=b0, link.function="logistic", rho=0, # UNIFORM COVARIATES
#            observational=TRUE, trt.p=0.5, a0=a0, details=TRUE) 

RHO <- c(0, 0.5)
N <- (1:10)*150
# N <- c(150, 900)
OBS <- c("Experimental", "Observational")
nrun <- 1000
BETA <- SE <- array(0, dim=c(nrun, 12, length(N), length(RHO), length(OBS))) 
COV <- array(0, dim=c(nrun, 16, 16, length(N), length(RHO), length(OBS))) 
set.seed(123)
for (m in 1:length(OBS)){
	observational <- ifelse(OBS[m]=="Experimental", FALSE, TRUE)
	for (j in 1:length(RHO)){
		rho <- RHO[j]; 
		for (k in 1:length(N)){
			n <- N[k]
			for (i in 1:nrun){
				print(cbind(study=OBS[m], rho=rho, n=n, run=i))
				dat <- rdat(n=n, b0=b0, link.function="logistic", rho=rho, # UNIFORM COVARIATES
					observational=observational, trt.p=0.5, a0=a0, details=FALSE)   # TRT ASSIGNMENT	
				out <- rModAna(formula=y~trt*(x1+x2+x3), data=dat, detail=FALSE, 
					  robust.VCOV=FALSE, adjust=TRUE)
				out.direct <- out$"out.direct"; out.inverse <- out$"out.inverse"; out.combined <- out$"out.combined"
				# out.direct; out.inverse; out.combined
				COV[i,,,k, j, m] <- out$V
				beta <- c(study=OBS[m], rho=rho, n=n,  as.numeric(as.character(out.direct[6:8, 1])), 
				          as.numeric(as.character(out.inverse[6:8, 1])),
				          as.numeric(as.character(out.combined$Estimate)))
				BETA[i,, k, j, m] <- beta
				se <- c(study=OBS[m], rho=rho, n=n, as.numeric(as.character(out.direct[6:8, 2])), 
				        as.numeric(as.character(out.inverse[6:8, 2])),
				        as.numeric(as.character(out.combined$'Std. Error')))
				SE[i,, k,j,m] <- se
}}}}
save(BETA, SE, COV, file="out-Weak.Rdat")



# ==============================
# NUMERICAL SUMMARY OF RESULTS
# ==============================

# BRING IN THE RESUTLS
# rm(list=ls(all=TRUE))
# load("out-n-Weak.Rdat")
# load("out-n-Weak-06JAN2021.Rdat")
# load("out-n-Strong.Rdat")

# ls()

n.study <- dim(BETA)[[5]]
n.rho <-  dim(BETA)[[4]]
ns <-  dim(BETA)[[3]]
# RHO <- c(0, 0.2, 0.5, 0.8)
# N <- (1:10)*100
OUT <- NULL
for (i in 1:n.study) {
  study <- ifelse(i==1, "Experimental", "Observational")
  for (j in 1:n.rho){
    rho <- RHO[j]
    for (k in 1:ns){
      n <- N[k]
      print(cbind(i, j, k, study=study, rho=rho, n=n))
      BETA0 <- as.matrix(BETA[,4:12,k,j,i]);
      SE0 <-  as.matrix(SE[,4:12,k,j,i])
      storage.mode(BETA0) <- storage.mode(SE0) <-"numeric"
      beta.mean <- matrix(apply(BETA0, MARGIN=2, FUN=mean), nrow=3, byrow=FALSE)
      SD <- matrix(apply(BETA0, MARGIN=2, FUN=sd), nrow=3, byrow=FALSE)
      ASE <- matrix(apply(SE0, MARGIN=2, FUN=mean), nrow=3, byrow=FALSE)
      COR <- cor(BETA0[,1:3], BETA0[, 4:6])
      Z <- pchisq((BETA0/SE0)^2, df=1, lower.tail=FALSE)
      POWER <- matrix(apply(Z<=0.05, 2, FUN=mean), nrow=3, byrow=FALSE)
      out <- cbind(Study=study, rho=rho, n=n, beta.mean, SD, ASE, POWER, COR)
      OUT <- rbind(OUT, out)
    }
  }
}
OUT <- as.data.frame(OUT)
methods <- c("Direct", "Inverse", "Combined")
colnames(OUT) <- c("Study", "rho", "n", paste("beta.", methods, sep=""),
                   paste("SD.", methods, sep=""), paste("ASE.", methods, sep=""),
                   paste("size/power.", methods, sep=""), 
                   paste("correlation", 1:3, sep=""))
head(OUT)
write.csv(OUT, file="result-Weak.csv", row.names =FALSE)




# ====================
# PLOT OF THE RESULTS
# ====================


# rm(list=ls(all=TRUE))
# load("out-Weak.Rdat")

RHO <- c(0, 0.5)
N <- (1:10)*150

n.study <- dim(BETA)[[5]]
n.rho <-  dim(BETA)[[4]]
# RHO <- c(0, 0.2, 0.5, 0.8)
for (k in 1:n.study) {
  study <- ifelse(k==1, "exp", "obs")
  for (m in 1:n.rho){
    rho <- RHO[m]
    print(cbind(k=k, m=m, study=study, rho=rho))
    filename <- paste("fig-", study, "-", rho, "-Weak-JUN2021.eps", sep="")
    BETA0 <- BETA[,, , m, k]; SE0 <- SE[,, , m, k]
    
    # SD & CORR
    ns <- dim(BETA0)[[3]] 
    N <- (1:ns)*150
    SD <- matrix(0, ns, 9)
    CORR <- matrix(0, ns, 3)
    for (i in 1:ns){
      Beta <- as.matrix(BETA0[, 4:12,i]) 
      storage.mode(Beta) <- "numeric"
      SD[i, ] <- apply(Beta, 2, sd)
      CORR[i, 1] <- cor(Beta[,1], Beta[,7], method="pearson")
      CORR[i, 2] <- cor(Beta[,2], Beta[,8], method="pearson")
      CORR[i, 3] <- cor(Beta[,3], Beta[,9], method="pearson")
    }
    # AVERAGED SE
    AvgSe <- matrix(0, ns, 9)
    for (i in 1:ns){
      Se <- as.data.frame(SE0[, 4:12,i]) 
      Se <- as.data.frame(sapply(Se, function(x) as.numeric(as.character(x)))) 
      AvgSe[i, ] <- apply(Se, 2, mean)
    }
    # PREPARE BOXPLOT DATA WITH SE
    dat0 <- NULL
    for (i in 1:ns) dat0 <- rbind(dat0, as.matrix(SE0[, ,i])) 
    dat0 <- as.data.frame(dat0)
    dat0[, 4:12] <- as.data.frame(sapply(dat0[, 4:12], function(x) as.numeric(as.character(x)))) 
    names(dat0) <- c("study", "rho", "n", paste("x", 1:9, sep=""))
    # dat0$n <- ordered(dat0$n, levels =as.character((1:10)*100))
    dat0$n <- as.numeric(as.character(dat0$n))
    dat0$rd1 <- (dat0$x1 - dat0$x7)/dat0$x1
    dat0$rd2 <- (dat0$x2 - dat0$x8)/dat0$x2
    dat0$rd3 <- (dat0$x3 - dat0$x9)/dat0$x3
    
    # -------------
    # 2 x 3 PLOT 
    # -------------
    postscript(file=filename, horizontal=TRUE)
    par(mfrow=c(2, 3), mar=c(4, 4, 3, 2))
    # SD
    for (j in 1:3){
      # BOXPLOTS OF SE
      form <- as.formula(paste(paste("x", j, sep=""), " ~ ", "n"))
      boxplot(form, data=dat0, boxwex=0.2, border="coral1", col="coral1", xlab="n", ylab="", notch=TRUE, 
              xaxt="n", outline=FALSE, at=(1:ns)-0.15, add=FALSE, cex.lab=1.2, lty=1)
      if(j==1) mtext(text="SD & SEs", side=2, line=2.1, col="black", cex=0.8)
      if (j==3) legend(6, 2.0, fill=c("coral3", "cadetblue3"), border=c("coral3", "cadetblue3"), 
                       legend=c("direct", "combined"), cex=1.2)
      beta.j <- substitute(list(hat(beta))[list(3,j0)], list(j0=j))
      text(5, par("usr")[4] + 0.1, srt=0, adj = 0, labels=beta.j, xpd = TRUE, col="blue", cex=1.5)
      # mtext(text=beta.j, side=4, line=1, col="blue", cex=2)
      axis(1, at=1:10, labels=N, tick=FALSE , cex=0.3)
      form0 <- as.formula(paste(paste("x", j+6, sep=""), " ~ ", "n"))
      boxplot(form0, data=dat0, boxwex=0.2, border="cadetblue2", col="cadetblue2", xaxt="n", notch=TRUE, 
              add=TRUE, at=(1:ns)+0.15, outline=FALSE, lty=1)
      # Averaged SE
      # lines(1:ns, AvgSe[1:ns,j], col="tomato", lty=1, lwd=0.5)
      # lines(N, AvgSe[1:ns, j+3], col="green4", lwd=0.5, lty=2)
      # lines(1:ns, AvgSe[1:ns, j+6], col="skyblue2", lwd=0.5, lty=1)
      # SD
      lines(1:ns, SD[1:ns, j], col="coral4", type="l", lwd=1, ylab="SD", xlab="n", cex=0.5)
      grid()	
      # lines(N, SD[1:ns, j+3], col="green4", lwd=0.5)
      lines(1:ns, SD[1:ns, j+6], col="cadetblue4", lwd=1, type="l", pch=20, cex=0.5)
    }
    
    for (j in 1:3){
      # BOXPLOT OF RELATIVE DIFFERENCE IN SE
      form <- as.formula(paste(paste("rd", j, sep=""), " ~ ", "n"))
      boxplot(form, data=dat0, boxwex=0.35, border="darkseagreen", col="palegreen2", 
              xlab="n", ylab="", notch=TRUE, cex.lab=1.2, lty=1, 
              xaxt="n", outline=FALSE, at=(1:ns), add=FALSE)
      abline(h=0, col="gray35", lty=1)
      axis(1, at=1:10, labels=N, tick=FALSE , cex=0.3)
      if(j==1) mtext(text="Relative Diff in SD", side=2, line=2.1, col="black", cex=0.8)
      if(j==2) mtext(text="Correlation", side=2, line=2.2, col="gray65", cex=0.8)
      
      # RELATIVE DIFFERENCE 
      rd.sd <- (SD[1:ns, j]-SD[1:ns, j+6])  /SD[1:ns, j]
      lines(1:ns, rd.sd, col="darkseagreen4", type="b", lwd=1.5); 
      grid()
      # rd.ase <- (AvgSe[1:ns, j]-AvgSe[1:ns, j+6])  /AvgSe[1:ns, j]
      # lines(1:ns, rd.ase, col="blue", type="b", lwd=1.5)
      par(new=T)
      plot(1:ns, CORR[ , j], lwd=2, col="gray75", xlab="", ylab="", type="b", 
           ylim=c(0.925, 0.995), pch=20, cex=0.8, cex.lab=1.2, lty=1, axes=F)
      axis(4,col="gray75", col.ticks="gray75", col.axis="gray75")
      
    }
    dev.off()
  }
}






#
