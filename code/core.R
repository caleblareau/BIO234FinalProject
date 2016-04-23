#' Packages
library(corrplot)
library(readr)

#' Read data
#+ cache = TRUE
path.to.input <- "/Users/lareauc/Desktop/Harvard University/Coursework/Spring 2016 Coursework/BIO 234-Data Structures/BIO234FinalProject/input"
hap <- "ALL_1000G_phase1interim_jun2011_chr20_impute.hap"
samp <- "ALL_1000G_phase1interim_jun2011.sample"
legend <- "ALL_1000G_phase1interim_jun2011_chr20_impute.legend"

hapdat <- read_delim(paste(path.to.input, hap,sep="/"), col_names = FALSE, delim = " ")
#legdat <- read.table(paste(path.to.input, legend, sep="/"), header=TRUE) # don't actually need this
sampdat <- read.table(paste(path.to.input, samp, sep="/"), header=TRUE)

#' Take a random subset of 'r' variants to speed up Var/Cov
#+ cache = TRUE
r <- 1000
gRand <- hapdat[sample(ncol(hapdat), r), ]
G <- as.matrix(gRand, ncol=2188)

#' The data is formatted where each pair of columns represents one person
#' but we want to combine them for the Var/Cov matrix
genetic.matrix <- G[ , c(TRUE,FALSE)] + G[ , c(FALSE,TRUE)]
var.covar <- var(genetic.matrix)

#' Get groups per haplotype; two haplotypes (columns) per sample
levels <- c("ASW", "CEU", "CHB", "CHS", "CLM", "FIN", "GBR", "IBS", "JPT", "LWK", "MXL", "PUR", "TSI", "YRI")
haplogroups <- unlist(rep(sampdat$population, each=2),use.names=FALSE)

#' Compute Transition probabilities
top <- G[1,]
hapdat2 <-  rbind(G[-1,],top) #shift first row to bottom
sameTrans <- G == G # determine number of times i,j = i+1, j; last row won't be informative. 

trans <-sapply(levels, function(t){ #separate matrices based on groups
    idx <- t==haplogroups
    init <- 1-sum(top[idx])/sum(idx) #proportion of 0s in first position
    tP <- rowSums(sameTrans[,t==haplogroups])/sum(idx) #transition probabilities
    tP <- tP[1:(length(tP)-1)] #remove last row
    c(init, tP) # add initial probability to top with transitions following
})

write.table(trans, "transition-probabilities_G.txt", row.names = FALSE, quote = FALSE)

#' Corrplot
#+ cache=TRUE
png("subsamp.png")
corrplot(var.covar[1:10,1:10], method="color", order="hclust", addrect=1)
dev.off()
