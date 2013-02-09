rm(list=ls())

#############
## functions
#############
# barcode and likelihood

`compute.barcode` <- function(dataset, adjmatrix)
{	
	edges <- which(adjmatrix == 1) ## a vector of indices
	barcode <- do.call(rbind, mclapply(edges, function(edge, nn, theweights) {
						j <- edge%%nn ## index 1:100
						if(j == 0)
							j <- nn   ## if j == 0, replace 0 by 100 
						
						i <- ceiling(edge/nn)  ## column index of genes
						
						c(i, j,dataset[, i] > dataset[, j])   ## select columns of dataset, j is an element of 1:100 if nn=100
					}, nn=ncol(adjmatrix),mc.cores = 8))
	barcode <- barcode[barcode[ ,1] != barcode[ ,2], , drop=FALSE]  ## removes duplicated rows so gene1_gene2 and gene2_gene1 only have one
	rownames(barcode) <- paste(rownames(adjmatrix)[barcode[ ,1]], rownames(adjmatrix)[barcode[ ,2]], sep="_")  ## rownames gene1_gene2
	barcode <- barcode[ ,-c(1,2), drop=FALSE]  
	return(barcode)
}

`compute.barcode.proba` <- function(dataset, theweight, adjmatrix)
{
	edges <- which(adjmatrix == 1)
	barcode <- do.call(rbind, mclapply(edges, function(edge, nn) {
						j <- edge%%nn
						if(j == 0)
							j <- nn        
						i <- ceiling(edge/nn)
						c(i, j, weighted.mean(dataset[, i] > dataset[, j], w = theweight, na.rm=TRUE)) #weighted.mean function:  be really careful (will extract function and include it in code) because if you don't assign weights it will not crash, so you sometimes don't know what exactly it is doing (at least I did not)
					}, nn=ncol(adjmatrix),mc.cores = 8))
	barcode <- barcode[barcode[ ,1] != barcode[ ,2], , drop=FALSE]
	rownames(barcode) <- paste(rownames(adjmatrix)[barcode[ ,1]], rownames(adjmatrix)[barcode[ ,2]], sep="_")
	# barcode <- barcode[ ,-c(1,2), drop=FALSE]
	#remove rows with all zeros
	return(barcode)
}

# Returns the product of likelihoods for a given barcode and sample
`compute.likelihood` <- function(barcode, sample)
{
	likelihood <- apply(barcode, 1, function(feature) {
				if (sample[feature[1]] > sample[feature[2]])
					return(feature[3])
				else
					return(1 - feature[3])
			})
	return(sum(log10(likelihood), na.rm=TRUE) / sum(!is.na(likelihood)))
}

# Returns the difference between one subtype and the mean, median, maximum, or minimum of all other subtypes for each gene pair
`top.features` <- function(proba,n,top, m){
	cols1 <- (1:ncol(proba))[1:ncol(proba)%%n==0]
	features <- do.call(cbind,lapply(cols1, function(x) { 
						abs(proba[,x] - m(proba[,cols1[-x/n]]))
	}))
	top_features <- do.call(cbind,lapply(lapply(1:ncol(features), function(x){sort(features[,x], TRUE)[1:top]}), names))	
	return(top_features)
}

# Maximum Likelihood
`compute.maxlikelihood` <- function(barcode, subtype, features, n, weights){
	likelihood <- lapply(1:ncol(features),function(x){subtype[features[,x],x*n]*barcode[features[,x],]})
	sum_likelihood <- do.call(cbind,(lapply(likelihood,function(x){apply(log10(x),2,sum,na.rm=TRUE)})))
	norm_likelihood <- sum_likelihood-rowSums(sum_likelihood)[row(sum_likelihood)]
	max_likelihood <- apply(norm_likelihood, 1, function(x){which(x == max(x,na.rm=TRUE))})
	subtype_names <- colnames(subtype.weights)
	max_subtype <- do.call(cbind,lapply(max_likelihood,function(x){max_likelihood[max_likelihood == x] <- subtype_names[x]}))
	return(max_subtype)
}

# Jaccard Distance
`compute.distance` <- function(subtype.classif, test, binary_barcode, subtype_features,subtype.weights, dist.type){
    subtype.classif_type <- lapply(1:length(levels(subtype.classif)), function(x){subtype.classif[subtype.classif == x]})
    binary_barcode_type <- lapply(1:length(subtype.classif_type),function(x){binary_barcode[subtype_features[,x],]})
    test <- lapply(1:length(subtype.classif_type), function(x){test[subtype_features[,x]]})
    train <- lapply(1:length(subtype.classif_type),function(x){cbind(test[[x]],binary_barcode_type[[x]])})
    results <- lapply(1:length(subtype.classif_type),function(x){lapply(2:ncol(train[[x]]),function(y){dist.type(t(train[[x]][,c(1,y)]))})})
    # results <- lapply(1:length(subtype.classif_type),function(x){lapply(2:ncol(train[[x]]),function(y){dist.type(t(train[[x]][,c(1,y)]))})})
    final <- lapply(1:ncol(subtype.weights),function(x){(unlist(results[[x]]))*subtype.weights[,x]}) 
    maximum <- lapply(1:ncol(subtype.weights),function(x){max(unlist(final[x]))})
    final_maximum <- max(unlist(maximum))
    final_result <- which(final_maximum == maximum)
return(final_result)
}

# Similarity value is the set intersection.
dist.standard <- function( adj ) {
	w <- adj %*% t( adj )
	## normalize
	w <- w / max( w )
	
	return( w[1,2] )
}

# Jaccard Distance
# Similarity is the set intersection over the union
dist.jaccard <- function( adj ) {
	w <- vegdist( adj, method="jaccard", upper=TRUE, diag=TRUE )
	w <- as.matrix( w )
    # convert to similarity
    w <- 1 - w
    
    return( w[1,2] )
}

## Pearson's correlation
## Similarity between two items is relative to the mean of all other
## pair-wise similarities.
dist.pearson <- function( adj, cutoff=0 ) {
	adj <- t( adj )
	w <- cor( adj, method="pearson" )
	
	## remove edges less than cutoff
	w[ w <= cutoff ] <- 0
	
	return( w[1,2])
}

dist.hamming <- function( adj ) {
	w <- vegdist( adj, method="hamming", upper=TRUE, diag=TRUE )
	w <- as.matrix( w )
    # normalize
    # w <- w / max( unlist(w) )
    # convert to similarity
    w <- 1 - w
    
    return( w[1,2] )
}

# dist.hamming <- function( adj ) {
#     w <- hamming.distance( adj )
#     ## normalize
#     # w <- w / max( unlist(w) )
#     ## convert to similarity
#     w <- 1 - w
#     
#     return( w[1,2] )
# }


#############
## libraries
#############

#library(devtools)
#install_github("mRMRe", username="bhaibeka", branch="master")
#system("chmod -R 775 /stockage/Laboratoires/HAI/Rlib/mRMRe")
library(mRMRe)

library(genefu)
library(jetset)

library(glmnet)
library(multicore)
library(parallel)
library(vegan)


#############
#############

# training set
# let's use EXPO for training our subtype classification model
# load EXPO dataset
load(file.path("../data","EXPO.RData"))
# use jetset to get one probe per gene
# jetset results in csv
jetset_csvs <- c("jetset.scores.hgu95av2_1.2.0.csv", "jetset.scores.hgu133a_1.2.0.csv", "jetset.scores.hgu133plus2_1.2.0.csv", "jetset.scores.u133x3p_1.2.0.csv")
names(jetset_csvs) <- c("affy.u95", "affy.u133a", "affy.u133plus2", "affy.3x")
# read jetset csv
jetset_annots <- NULL
for(i in 1:length(jetset_csvs)) {
    dd <- read.csv(file.path("jetset_csv", jetset_csvs[i]), stringsAsFactors=FALSE)
    rownames(dd) <- dd[ ,"probeset"]
    jetset_annots <- c(jetset_annots, list(dd))
}
names(jetset_annots) <- names(jetset_csvs)
# jetset mapping
aa <- jetset_annots[["affy.u133plus2"]]
nn <- intersect(colnames(data), rownames(aa))
if(is.null(nn) || length(nn) < 1) { stop(sprintf("bad platform for dataset %s!", ddn[k])) }
## consider only the best probes
nn <- nn[aa[nn, "best"]]
aa <- aa[nn, ,drop=FALSE]
annot <- data.frame("probe"=nn, "gene.symbol"=aa[nn,"symbol"], "EntrezGene.ID"=aa[nn,"EntrezID"])
rownames(annot) <- nn
data <- data[ ,nn,drop=FALSE]
# rename probeset with gene id
colnames(data) <- rownames(annot) <- paste("geneid", as.character(annot[ ,"EntrezGene.ID"]), sep="_")


# fit SCMGENE on EXPO
pdf("scmgene_expo.pdf")
modgene <- list("ESR1"=cbind("probe"="geneid_2099", "EntrezGene.ID"=2099, "coefficient"=1), "ERBB2"=cbind("probe"="geneid_2064", "EntrezGene.ID"=2064, "coefficient"=1), "AURKA"=cbind("probe"="geneid_6790", "EntrezGene.ID"=6790, "coefficient"=1))
myscmgene <- subtype.cluster(module.ESR1=modgene$ESR1, module.ERBB2=modgene$ERBB2, module.AURKA=modgene$AURKA, data=data, annot=annot, do.mapping=FALSE, do.scale=TRUE, plot=TRUE)
dev.off()

# probabilities to belong to each of the 4 breast cancer molecular subtypes
# ER-/HER2-, HER2+, ER+/HER2- High Prolif, and ER+/HER2- Low Prolif
subtype.weights <- myscmgene$subtype.proba2
subtype.classif <- factor(myscmgene$subtype2, levels=c("ER-/HER2-", "HER2+", "ER+/HER2- High Prolif", "ER+/HER2- Low Prolif"))

# transform the continuous gene expressions in binary values (gene x > gene y)
# the set of binary values are referred to as barcode here after
# build an adjacency matrix that represent the pairs of genes one needs to consider
aa <- matrix(0, ncol=ncol(data), nrow=ncol(data), dimnames=list(colnames(data), colnames(data)))

# test with randomly selected 500 genes
set.seed(1)
aa[sample(1:100,replace = FALSE), sample(1:100,replace = FALSE)] <- 1

# binary barcode
binary_barcode <- compute.barcode(dataset=data, adjmatrix=aa)
binary_barcode[binary_barcode == 0] <- NA  # Made binary code equal to zero NA for maximum likelihood:  may want to use pseudocounts here (pseudocounts did not work)

# probability barcode for each subtype
subtype_barcode <- do.call(cbind,lapply(1:ncol(subtype.weights), function(x) { compute.barcode.proba(dataset = data, theweight = subtype.weights[,x], adjmatrix = aa)}))

# feature selection
subtype_features <- top.features(proba = subtype_barcode, n = 3, top = 50, m = rowMeans)

# maximum likelihood
maximum_likelihood <- compute.maxlikelihood(barcode = binary_barcode, subtype = subtype_barcode,features = subtype_features, n = 3, weights = subtype.weights)

# comparing subtypes: Maximum Likelihood
comp1 <- ifelse(maximum_likelihood == subtype.classif,1,0)
counts1 <- table(comp1)
accuracy1 <- (counts1[2]/(counts1[1]+counts1[2]))

# linear model to real expression values
# change levels from "ER-/HER2-", "HER2+", "ER+/HER2- High Prolif", "ER+/HER2- Low Prolif" to #1-4
for (i in 1:length(levels(subtype.classif))){	# use lapply here instead of for loop 
	levels(subtype.classif)[i] <- i
}

# glmnet fit:  Fit a generalized linear model via penalized maximum likelihood
fit1 <- glmnet(data, subtype.classif, family = "multinomial")

# cluster: ones and zeros for quality check -- need to do this! 

# distance
binary_barcode <- compute.barcode(dataset=data, adjmatrix=aa)
subtype.classif_type <- lapply(1:length(levels(subtype.classif)), function(x){subtype.classif[subtype.classif == x]})

# use testing set as training set
test <- binary_barcode[,1:353]
final <-  lapply(1:ncol(test), function(x){compute.distance(subtype.classif,test[,x],binary_barcode[,-test[,x]],subtype_features,subtype.weights[rownames(subtype.weights) %in% colnames(binary_barcode[,-test[,x]]),], dist.pearson)})

#comparing subtypes: distance
expected <- as.integer(lapply(1:ncol(test),function(x){as.integer(subtype.classif[names(subtype.classif) == colnames(binary_barcode)[x]])}))
comp2 <- lapply(1:ncol(test), function(x){ifelse(unlist(final[x]) == expected[x],1,0)})
counts2 <- table(unlist(comp2))
accuracy2 <- (counts2[2]/(counts2[1]+counts2[2]))

	