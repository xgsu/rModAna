

################################################
# REFINED MODERATION ANALYSIS OF WART DATA
################################################

# rm(list=ls(all=TRUE))
dat <- read.csv(file="wart.csv", header=TRUE)
# MERGE TYPE 1&2 TOGETHER; type=type3
dat$type <- ifelse(dat$type==3, 1, 0)

source("Functions-RMA.R")
form.direct <- response ~ cryo*(age + type + nwarts)
out <- rModAna(formula=form.direct, data=dat, detail=FALSE, 
               robust.VCOV=FALSE, adjust=TRUE)
out.direct <- out$"out.direct"; 
out.inverse <- out$"out.inverse"; 
out.combined <- out$"out.combined"
out.direct; out.inverse; out.combined

out.temp <- out.direct
out.temp[1:5,] <- NA 
out.temp[6:8,] <- as.matrix(out.combined) 
RESULT <- cbind(out.direct, out.inverse, out.temp)
RESULT
write.csv(RESULT, file="result-wart.csv", row.names=TRUE)

# EXPLORE
dat$agegrp <- cut(dat$age, breaks=quantile(dat$age, probs = seq(0, 1, 0.25)), 
                  labels = 1:4)
table(dat$response, dat$cryo, dat$agegrp)
table(dat$response, dat$cryo, dat$type)

