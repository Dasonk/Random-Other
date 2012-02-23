basedir <- "~/RA/PPDK"
setwd(basedir)

## Read in ncr background noise for HL data
clHL <- read.table("covlist_HL", sep = "\t", header = FALSE)
clHL[, 1] <- gsub(".out", "", substring(clHL[,1], 6))
clHL <- clHL[,1:3]

## Read in ncr background noise for LL data
clLL <- read.table("covlist_LL", sep = "\t", header = FALSE)
clLL[, 1] <- gsub(".out", "", substring(clLL[,1], 6))

## Row 8 is messed up
clLL <- clLL[-8, 1:3]
## Reset row names
rownames(clLL) <- NULL

ncr <- rbind(clHL, clLL)
ncr$rate <- ncr[,3]/ncr[,2]
colnames(ncr) <- c("TREATMENT", "NCRLENGTH", "NCRCOUNT", "RATE")

write.csv(ncr, file = "NCRrates.csv", row.names = FALSE)
