basedir <- "~/RA/PPDK/PPDK_lw"
rdir <- paste(basedir, "/ReducedData/", sep = "")

files <- dir(rdir)
i <- 1
j <- read.csv(paste(rdir, files[i], sep = ""), header = TRUE)
n <- dim(j)[1]
m <- length(files)
trans <- j[,1]

data <- matrix(nrow = n, ncol = m)
for(i in seq(files)){
##for(i in 1:3){
    cat(i, "/", m, "\n")
    flush.console()
    j <- read.csv(paste(rdir, files[i], sep = ""), header = TRUE)
    ## Check to make sure transcript names match
    if(sum(j[,1] != trans) > 0){
        warning(paste("Transcript names don't match on file:", file[i]))
    }
    data[,i] <- j$COVERAGE
}


colnames(data) <- gsub(".csv", "", files)
rownames(data) <- trans

write.csv(data, file = paste(basedir, "/FullData.csv", sep = ""))
