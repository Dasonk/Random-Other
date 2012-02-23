
basedir <- "~/RA/PPDK/cuffdiff2/"
setwd(basedir)

# Contains total read counts for each sample
totalcounts <- read.csv("aligned_read_count.csv", header = FALSE, col.names = c("sample", "counts"))
# Should give the gene length for each gene
lengths <- read.csv(dir("~/RA/PPDK/PPDK_lw/ReducedData", full=T)[1])[,c("TRANSCRIPT", "LENGTH")]

# get rates
rates <- read.csv("~/RA/PPDK/NCRrates.csv", header = TRUE)

samps <- grep("*_b", dir(), value = TRUE)

dat <- list()
for(i in seq_along(samps)){
    j <- read.table(paste(samps[i], "/genes.fpkm_tracking", sep = ""), header = T, sep = "\t")
    cols <- c("gene_id", grep("FPKM", names(j), value = T))
    dat[[i]] <- j[,cols]
}

getTotalCount <- function(nm){
    nm <- names(dat[[1]])[2]
    nm <- gsub("_FPKM", "", nm)
    light <- substr(nm, 1, 2)
    rep <- substr(nm, 3, 3)
    tmp <- substring(nm, 4)
    n <- nchar(tmp)
    type <- substr(tmp, 1, n-1)
    loc <- toupper(substr(tmp, n, n))
    loc <- switch(loc, "B" = "BASE", "T" = "Tip", "1" = "1cm", "4" = "4cm")
           
    j <- toupper(paste("PPDK_", light, rep, "_", type, "_", loc, sep = ""))
    k <- toupper(totalcounts$sample)
    id <- grep(j, k)
    return(totalcounts[id, "counts"])
}

k <- dat[[1]][,c(1,2)]
k$gene_id <- as.character(k$gene_id)

lengths$TRANSCRIPT <- as.character(lengths$TRANSCRIPT)

m <- merge(lengths, k, by.x = "TRANSCRIPT", by.y="gene_id")