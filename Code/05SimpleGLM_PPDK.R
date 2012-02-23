library(plyr)
basedir <- "~/RA/PPDK/PPDK_lw/"
datfile <- "~/RA/PPDK/PPDK_lw/SigDataPPDK.csv"
genedat <- read.csv(datfile, header = TRUE, row.names = 1)
genedat <- round(genedat)

## Create a data frame with the conditions for each column
nm <- data.frame(nm = colnames(genedat), stringsAsFactors = FALSE)
tmp <- as.data.frame(do.call(rbind, strsplit(nm$nm, "_")))
names(tmp) <- c("Light", "Type", "Loc")
tmp$Rep <- factor(substring(tmp$Light, 3,3))
tmp$Light <- factor(substring(tmp$Light, 1, 2))

block <- factor(paste(tmp$Light, tmp$Rep, sep = ""))

jfull <- with(tmp, paste(Light, Type, Loc))

testvals <- expand.grid(c("HL", "LL"), c("1cm", "4cm", "Base", "Tip"))
testlist <- list()
apply(testvals, 1, function(x){
    idx <- paste(x, collapse = "")
    tmp <- jfull
    tmp[grep( paste(x, collapse = ".*"), jfull)] <- idx
    testlist[[idx]] <<- factor(tmp)
    invisible()
})



## Read in all the data we can get the libsizes
## Note we're assuming that fulldat and genedat have
## the columns representing the same samples.
fulldat <- read.csv("~/RA/PPDK/PPDK_lw/FullData.csv", row.names=1, header = T)

offsetfun <- function(x){
    ans <- log(quantile(x[x != 0], .75))
    return(ans)
}

offset <- as.numeric(apply(fulldat, 2, offsetfun))

jfull <- factor(jfull)

ngenes <- dim(genedat)[1]
cnt <- 1
analyzeGLM <- function(y){
    cat(round(cnt/ngenes, digits = 4), "%\n")

    y <- as.numeric(y)
    o.full <- glm(y ~ jfull + block,
                  family = poisson(link = "log"),
                  offset = offset)

    tmpphihat <- deviance(o.full)/df.residual(o.full)
    phihat <- max(1, tmpphihat)

    pvals <- lapply(testlist,
        function(x){
            o <- glm(y ~ x + block,
                 family = poisson(link = "log"),
                 offset = offset)
            a <- anova(o, o.full)
            Fstat <- (a[2,4] / a[2,3]) / phihat
            pval <- 1 - pf(Fstat, a[2,3], a[2,1])
            return(pval)
        }
    )

    cnt <<- cnt + 1
    return(unlist(pvals))
}

start <- Sys.time()
output <- apply(genedat, 1, analyzeGLM)
end <- Sys.time()
print(end - start)

out <- t(output)
write.csv(out, file = paste(basedir, "GLMpvals.csv", sep = ""))

