
library(plyr)
basedir <- "~/RA/PPDK/PPDK_lw/"

pvals <- read.csv(file = paste(basedir, "GLMpvals.csv", sep = ""),
                  header = TRUE,
                  row.names = 1)

qvals <- apply(pvals, 2, p.adjust, method = "BH")


cutoffs <- c(.001, .01, .05, .1, .15, .20)

out <- matrix(NA, nrow = length(cutoffs), ncol = dim(qvals)[2])
colnames(out) <- colnames(qvals)
rownames(out) <- cutoffs

getCount <- function(dat, cut){
    return(sum(dat <= cut))
}

for(i in seq_along(cutoffs)){
    for(j in seq(dim(qvals)[2])){
        out[i, j] <- getCount(qvals[,j], cutoffs[i])
    }
}

## Plots

par(mfrow = c(2,2))
for(i in 1:4){
    hist(pvals[,i], main = colnames(pvals)[i], xlab = "P-Value")
}

windows()
par(mfrow = c(2,2))
for(i in 5:8){
    hist(pvals[,i], main = colnames(pvals)[i], xlab = "P-Value")
}
