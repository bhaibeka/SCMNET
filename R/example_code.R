rm(list=ls())

#############
## functions
#############
# barcode and likelihood

`compute.barcode` <- function(dataset, adjmatrix)
{
	edges <- which(adjmatrix == 1)
	barcode <- do.call(rbind, lapply(edges, function(edge, nn) {
						j <- edge%%nn #modulus
						if(j == 0)
							j <- nn         
						i <- ceiling(edge/nn)
						c(i, j,dataset[, i] > dataset[, j])
					}, nn=ncol(adjmatrix)))
	barcode <- barcode[barcode[ ,1] != barcode[ ,2], , drop=FALSE]
	rownames(barcode) <- paste(rownames(adjmatrix)[barcode[ ,1]], rownames(adjmatrix)[barcode[ ,2]], sep="_")
	barcode <- barcode[ ,-c(1,2), drop=FALSE]
	return(barcode)
}

`compute.barcode.proba` <- function(dataset, weights, adjmatrix)
{
	edges <- which(adjmatrix == 1)
	barcode <- do.call(rbind, lapply(edges, function(edge, nn, weights) {
						j <- edge%%nn
						if(j == 0)
							j <- nn        
						i <- ceiling(edge/nn)
						c(i, j, weighted.mean(dataset[, i] > dataset[, j], w=weights, na.rm=TRUE))
					}, nn=ncol(adjmatrix)))
	barcode <- barcode[barcode[ ,1] != barcode[ ,2], , drop=FALSE]
	rownames(barcode) <- paste(rownames(adjmatrix)[barcode[ ,1]], rownames(adjmatrix)[barcode[ ,2]], sep="_")
	# barcode <- barcode[ ,-c(1,2), drop=FALSE]
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
	
		
#############
## libraries
#############

#library(devtools)
#install_github("mRMRe", username="bhaibeka", branch="master")
#system("chmod -R 775 /stockage/Laboratoires/HAI/Rlib/mRMRe")
library(mRMRe)

library(genefu)
library(jetset)


#############
#############

# training set
# let's use EXPO for training our subtype classification model
# load EXPO dataset
load(file.path("data","EXPO.RData"))
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

# test with on;y the first 100 genes
aa[1:100, 1:100] <- 1

# barcode
binary_barcode <- compute.barcode(dataset=data, adjmatrix=aa)
probability_barcode <- compute.barcode.proba(dataset=data, adjmatrix=aa)

#Select a feature set of gene pairs