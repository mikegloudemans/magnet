#!/usr/bin/Rscript
# Author: Mike Gloudemans
# Date: 8/23/2016

# Note: I'm trying to follow Bruna's pipeline for this

library(peer)
library(preprocessCore)

# Bruna's normalization function
InverseGaussian=function(x, mean=0, sd=1){
  L <- length(x)
  normal.x <- qnorm( rank(x)/(L+1), mean = mean, sd = 1)
}

data = read.table("/users/zappala/projects/sardinia/data/expression/Fpkm.GeneLevel.20140927.tsv", row.names=1, header=TRUE)

# Filter out lowly expressed genes. Would like to know more about best
# practices for this, but for now I'll just go with something that seems reasonable.
# This is actually pretty conservative right now.
data <- data[rowSums(data > 5) > 10,]

# Set random seed and select a random subset of ~30% of these genes. Would like
# to know more about why this is done, but for now I'll just accept it.
set.seed(1)
data <- data[sample(1:dim(data)[1], dim(data)[1]/3),]

# Save row names and column names
rnames <- rownames(data)
cnames <- colnames(data)

# Add pseudocounts for later log transformation
data <- data + 1

# Correct for library size
data <- t(t(as.matrix(data))/rowSums(t(as.matrix(data))))

# Log-transform and quantile normalize
data <- log(data)

qn <- normalize.quantiles(as.matrix(data),copy=TRUE)
qn <- apply(qn,1,InverseGaussian)

qn <- as.data.frame(qn)
rownames(qn) <- cnames
colnames(qn) <- rnames

# Run PEER to estimate covariates
model = PEER()
PEER_setPhenoMean(model,as.matrix(qn))
PEER_setNk(model,10)
#PEER_setAdd_mean(model, TRUE)

PEER_update(model)

# Get covariate factors and output to a table
factors = PEER_getX(model)
rownames(factors) <- cnames
write.table(factors, quote=FALSE, sep = "\t", col.names=FALSE, file="/users/mgloud/projects/motrpac/sardinia/output/PEER/PEER_covariates.txt")

# Get covariate weights and output to a table
weights = PEER_getW(model)
rownames(weights) <- rnames
write.table(weights, quote=FALSE, sep = "\t", col.names=FALSE, file="/users/mgloud/projects/motrpac/sardinia/output/PEER/PEER_weights.txt")
