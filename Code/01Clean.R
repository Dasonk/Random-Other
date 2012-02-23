basedir <- "~/RA/PPDK/PPDK_lw"
dir1 <- "~/RA/PPDK/PPDK_lw/lwoutput"
dir2 <- "~/RA/PPDK/PPDK_lw/lwoutput2"

files1 <- grep("gene_RPKM", dir(dir1, full.names = TRUE), value = TRUE)
files2 <- grep("gene_RPKM", dir(dir2, full.names = TRUE), value = TRUE)
files <- c(files1, files2)

fileidx <- seq(1, length(files), 2)
n <- tail(fileidx, 1)
for(i in fileidx){

    ## Grab name of stuff we're working with
    name <- gsub("_trim.*", "", gsub(".*PPDK_", "", files[i+1]))
    cat(i, "/", n, "- Operating on:", name, format(Sys.time()),"\n")
    flush.console()

    antisense <- read.table(files[i], sep="\t", header = TRUE, stringsAsFactors = FALSE)

    sense <- read.table(files[i+1], sep="\t", header = TRUE, stringsAsFactors = FALSE)

    ## Grab transcript and coverage
    ## Only grab length from antisense - assume they're the same
    ## for both
    antisense <- antisense[,c("TRANSCRIPT", "COVERAGE", "LENGTH")]
    sense <- sense[,c("TRANSCRIPT", "COVERAGE")]

    ## Merge to make sure transcript name matches
    j <- merge(sense, antisense, by.x = "TRANSCRIPT", by.y = "TRANSCRIPT")

    ## If merging ends up in different dimension we have a problem
    if(dim(j)[1] != dim(sense)[1]){
        warning(paste("Merged data has different dimension:", name))
    }

    ## Plot sense v. antisense
    png(filename = paste(basedir, "/LogPlots/", name, ".png", sep = ""))
    plot(log(j[,2]),
         log(j[,3]),
         main = name,
         xlab = "Log(Sense)",
         ylab= "Log(Antisense)",
         col=rgb(20, 10, 80, 50, maxColorValue = 255),
         pch = 16)
    abline(0, 1)
    dev.off()

    colnames(j)[2:3] <- c("COVERAGE.SENSE", "COVERAGE.ANTISENSE")

    ## We want the sum of
    j$COVERAGE <- j[, "COVERAGE.SENSE"] + j[, "COVERAGE.ANTISENSE"]
    ## j <- j[,c("TRANSCRIPT", "COVERAGE", "LENGTH")]

    write.csv(j, file = paste(basedir, "/ReducedData/", name, ".csv", sep = ""), row.names = FALSE)

}

## Check to make sure transcript names are the same
## j1 <- read.csv(paste(basedir, "/ReducedData/HL3_WT_Tip.csv", sep = ""), header = T)
## j2 <- read.csv(paste(basedir, "/ReducedData/LL1_heter_1cm.csv", sep = ""), header = T)
## sum(j1[,1] == j2[,1]) == dim(j1)[1]
