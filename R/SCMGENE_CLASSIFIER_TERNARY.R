rm(list=ls())

#############
## functions
#############
# barcode and likelihood

`compute.barcode` <- function(dataset, adjmatrix)
{	
	edges <- which(adjmatrix == 1) ## a vector of indices
	barcode1 <- do.call(rbind, lapply(edges, function(edge, nn, theweights) {
						j <- edge%%nn ## index 1:number of columns in adjmatrix
						if(j == 0)
							j <- nn   ## if j == 0, replace 0 by 49
						
						i <- ceiling(edge/nn)  ## column index of genes
                        c(i,j, dataset[,i] > dataset[,j] + sd(as.numeric(dataset)))  ## select columns of dataset, j is an element of 1:nn
					}, nn=ncol(adjmatrix)))
    barcode1 <- barcode1[barcode1[ ,1] != barcode1[ ,2], , drop=FALSE]  ## removes diagonal values gene1_gene1 gene1_gene1
	rownames(barcode1) <- paste(rownames(adjmatrix)[barcode1[ ,1]], rownames(adjmatrix)[barcode1[ ,2]], sep="-")  ## rownames gene1_gene2
    # barcode1 <- barcode1[ ,-c(1,2), drop=FALSE]
	barcode2 <- do.call(rbind, lapply(edges, function(edge, nn, theweights) {
						j <- edge%%nn ## index 1:number of columns in adjmatrix
						if(j == 0)
							j <- nn   ## if j == 0, replace 0 by 49
						
						i <- ceiling(edge/nn)  ## column index of genes
                        c(i,j,dataset[,i] < dataset[,j] - sd(as.numeric(dataset)))  ## select columns of dataset, j is an element of 1:nn
					}, nn=ncol(adjmatrix)))
    barcode2 <- barcode2[barcode2[ ,1] != barcode2[ ,2], , drop=FALSE]  ## removes diagonal values gene1_gene1 gene1_gene1
	rownames(barcode2) <- paste(rownames(adjmatrix)[barcode2[ ,1]], rownames(adjmatrix)[barcode2[ ,2]], sep="-") ## removes diagonal values gene1_gene1 gene1_gene1
    barcode <- barcode1 - barcode2 
    barcode <- barcode[ ,-c(1,2), drop=FALSE]
	return(barcode)
}

`compute.barcode.proba` <- function(dataset, theweight, adjmatrix)
{
	edges <- which(adjmatrix == 1) ## a vector of indices
	barcode1 <- do.call(rbind, lapply(edges, function(edge, nn, theweights) {
						j <- edge%%nn ## index 1:number of columns in adjmatrix
						if(j == 0)
							j <- nn   ## if j == 0, replace 0 by 49
						
						i <- ceiling(edge/nn)  ## column index of genes
                        c(i,j, dataset[,i] > dataset[,j] + sd(as.numeric(dataset)))  ## select columns of dataset, j is an element of 1:nn
					}, nn=ncol(adjmatrix)))
    barcode1 <- barcode1[barcode1[ ,1] != barcode1[ ,2], , drop=FALSE]  ## removes diagonal values gene1_gene1 gene1_gene1
	rownames(barcode1) <- paste(rownames(adjmatrix)[barcode1[ ,1]], rownames(adjmatrix)[barcode1[ ,2]], sep="-")  ## rownames gene1_gene2
    # barcode1 <- barcode1[ ,-c(1,2), drop=FALSE]
	barcode2 <- do.call(rbind, lapply(edges, function(edge, nn, theweights) {
						j <- edge%%nn ## index 1:number of columns in adjmatrix
						if(j == 0)
							j <- nn   ## if j == 0, replace 0 by 49
						
						i <- ceiling(edge/nn)  ## column index of genes
                        c(i,j,dataset[,i] < dataset[,j] - sd(as.numeric(dataset)))  ## select columns of dataset, j is an element of 1:nn
					}, nn=ncol(adjmatrix)))
    barcode2 <- barcode2[barcode2[ ,1] != barcode2[ ,2], , drop=FALSE]  ## removes diagonal values gene1_gene1 gene1_gene1
	rownames(barcode2) <- paste(rownames(adjmatrix)[barcode2[ ,1]], rownames(adjmatrix)[barcode2[ ,2]], sep="-") ## removes diagonal values gene1_gene1 gene1_gene1
    # barcode2 <- barcode2[ ,-c(1,2), drop=FALSE]
    barcode <- barcode1 - barcode2 
    barcode <- barcode[ ,-c(1,2), drop=FALSE]
    barcode <- apply(barcode, 1, function(x) {weighted.mean(x,w = theweight)})
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

# MRMR feature selection
#    features <- feature.selection(proba = subtype_barcode, barcode = binary_barcode, n = 3, top = 5, redundant = 0, w = 1)
`feature.selection` <- function(proba, barcode, n, top,redundant, w){
    redundancy <- list()
    features <- list()
    # subtype_freq <- proba[,seq(n,ncol(proba),n)]
    salience <- lapply(1:ncol(proba), function(x){apply(proba[,x] - proba[,-x] ,1,min)})
    # salience <- lapply(1:ncol(proba), function(x){abs(salience[[x]])})
        if(redundant == 1){
            # first genepair: MRMR depends highly on first gene pair
            gene_pair <- lapply(1:ncol(proba), function(x){sort(unlist(salience[[x]]), decreasing = TRUE)[1]})
            # features: make first gene pair first feature
            features[[1]] <- do.call(cbind, lapply(1:ncol(proba), function(x){names(gene_pair[[x]])}))
            i = 2
            while(i <= top){
                # redundancy[[i-1]] <- do.call(cbind,lapply(1:ncol(subtype_freq), function(x){lapply(1:nrow(barcode), function(y){(subtype_freq[,x][rownames(subtype_freq) == names(gene_pair[[x]])])*(1-vegdist(barcode[c(which(rownames(barcode) == names(gene_pair[[x]])),y),],methods = "jaccard"))})}))
                redundancy[[i-1]] <- do.call(cbind,lapply(1:ncol(proba), function(x){lapply(1:nrow(barcode), function(y){(vegdist(barcode[c(which(rownames(barcode) == names(gene_pair[[x]])),y),],methods = "jaccard", binary = FALSE))})}))
                # tmp1_redundancy: makes redundancy a matrix
                tmp1_redundancy <- abind(redundancy,along = 3)
                # avg_redundancy: calculate average redundancy between gene_pairs
                avg_redundancy <- apply(tmp1_redundancy,c(1,2),mean) 
                # tmp2_redundancy: makes redundancy a list
                tmp2_redundancy <- split(avg_redundancy,col(avg_redundancy))
                # tmp3_redundancy: set names 
                tmp3_redundancy <- lapply(1:ncol(proba),function(x){setNames(tmp2_redundancy[[x]], rownames(barcode))})
                # tmp4_redundancy: remove gene pairs in redundancy
                tmp4_redundancy <- lapply(1:ncol(proba), function(x){tmp3_redundancy[[x]][!names(tmp3_redundancy[[x]]) %in% lapply(features,`[[`,x)]})
                # tmp1_saliency: remove gene pairs in saliency
                tmp1_salience <- lapply(1:ncol(proba), function(x){salience[[x]][!names(salience[[x]]) %in% lapply(features,`[[`,x)]})
                # top_genes: saliency - redundancy
                top_genes <- lapply(1:ncol(proba), function(x){which.max(w*unlist(tmp1_salience[[x]])+unlist(tmp4_redundancy[[x]]))})
                # ifelse(i == 2, top_genes <- lapply(1:ncol(subtype_freq), function(x){which.max(unlist(tmp4_redundancy[[x]]))}),top_genes <- lapply(1:ncol(subtype_freq), function(x){which.max(unlist(tmp1_salience[[x]])-unlist(tmp4_redundancy[[x]]))}))
                features[[i]] <- do.call(cbind, lapply(1:ncol(proba), function(x){names(top_genes[[x]])}))
                gene_pair <- top_genes
                i = i + 1
            }
            features <- do.call(rbind,features)
    }
    else{
        features <- do.call(cbind,lapply(1:ncol(proba), function(x){names(sort(unlist(salience[[x]]), decreasing = TRUE))}[1:top]))
    }
    return(features)
}

# Finds standard binary vectors using the frequencies of the features selected in feature.selection: if frequency is > 0.5 it is 1, else 0
`standard.vectors` <- function(proba, features, n,threshold, dataset){
    # subtype_freq <- proba[,seq(n,ncol(proba),n)]
    # frequencies: finds the frequencies for each feature
    frequencies <- lapply(1:ncol(proba), function(x){proba[,x][rownames(proba) %in% features[,x]]})
    # standard_vector: if frequency is > 0.5 it is 1, else 0
    standard_vector <- lapply(1:ncol(proba), function(x){ifelse(frequencies[[x]] + sd(as.numeric(dataset))/mean(dataset) < 0, -1, ifelse(frequencies[[x]] - sd(as.numeric(dataset))/mean(dataset) > 0, 1, 0))})
    return(standard_vector)   
}

# Finds the jaccard distance (1 being close, 0 being far away) between the standard vectors for each subtype and the samples, and the finds the max and assigns the sample to that subtype
`jaccard.distance` <- function(proba, barcode, features,subtype.weights,n,standard){
    # subtype_freq <- proba[,seq(n,ncol(proba),n)]
    # samples:  selects the feature values for each sample
    samples <-lapply(1:ncol(proba), function(x){barcode[rownames(barcode) %in% features[,x],]})
    # tmp_samples: binds the standard vector to the feature values for each sample to do the jaccard
    # tmp_samples <- lapply(1:ncol(subtype_freq),function(y){cbind(standard[[x]],samples[[x]])})
    tmp_samples <- lapply(1:ncol(proba),function(x){lapply(1:ncol(proba), function(y){cbind(standard[[x]],samples[[y]])})})
    # jaccard: find the jaccard distance between standard vector and sample
    # jaccard <- do.call(cbind,lapply(1:ncol(subtype_freq), function(x){lapply(1:ncol(tmp_samples[[x]]), function(y){1-vegdist(t(tmp_samples[[x]][,c(1,y)]),methods = "jaccard")})}))
    jaccard <- lapply(1:ncol(proba),function(z){do.call(cbind,lapply(1:ncol(proba), function(x){lapply(1:ncol(tmp_samples[[z]][[x]]), function(y){vegdist(t(tmp_samples[[z]][[x]][,c(1,y)]),methods = "jaccard")})}))})
    # tmp1_jaccard: make jaccard matrix numeric
    tmp1_jaccard <- lapply(1:ncol(proba),function(x){apply(jaccard[[x]],2,as.numeric)})
    # remove the first row of jaccard as it is only a comparison of the standard vectors to themselves so jaccard distance equals 0
    tmp2_jaccard <- lapply(1:ncol(proba),function(x){tmp1_jaccard[[x]][-1,]})
    #assign colnames to jaccard
    # colnames(tmp2_jaccard) <- c(1,2,3,4)
    #max_likelihood:  maximum likelihood
    max_likelihood <- do.call(cbind,lapply(1:ncol(proba),function(x){tmp2_jaccard[[x]][,x]/rowSums(tmp2_jaccard[[x]])}))
    # maximum <- lapply(1:ncol(subtype_freq), function(x){max_likelihood[,x]*subtype.weights[,x]})
    colnames(max_likelihood) <- c(1,2,3,4)
    #predicted: predicted subtypes
    predicted <- colnames(max_likelihood)[apply(max_likelihood,1,which.min)]
    return(predicted)
}

`maximum_likelihood` <- function(tmp2_jaccard,binary_barcode,subtype_barcode)
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

library(pbapply)
library(abind)


#############
#############

# training set
# let's use EXPO for training our subtype classification model
# load EXPO dataset
load(file.path("../data","METABRIC_GE_ENTREZ.RData"))
load(file.path("../data","METABRIC_GE_SUBTYPING.RData"))
data <- tumor.ge # METABRIC dataset ONLY
colnames(data) <- sub("geneid.","",colnames(data))
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
aa <- jetset_annots[["affy.u133plus2"]] #EXPO is U133plus2
pam50 <- read.table("../data/pam50_genes.csv", sep = ",", header = TRUE) #read in PAM50 genes
data <- data[,colnames(data) %in% pam50$EntrezGene] #consider only PAM50 genes
colnames(data) <- aa$probeset[match(colnames(data),aa$EntrezID)] #match EntrezID to probeset
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
# data <- data[,colnames(data) %in% c("geneid_2099", "geneid_2064", "geneid_6790")] #using only SCMGENE genes


# fit SCMGENE on EXPO
# pdf("scmgene_metabric.pdf")
# modgene <- list("ESR1"=cbind("probe"="geneid_2099", "EntrezGene.ID"=2099, "coefficient"=1), "ERBB2"=cbind("probe"="geneid_2064", "EntrezGene.ID"=2064, "coefficient"=1), "AURKA"=cbind("probe"="geneid_6790", "EntrezGene.ID"=6790, "coefficient"=1))
# myscmgene <- subtype.cluster(module.ESR1=modgene$ESR1, module.ERBB2=modgene$ERBB2, module.AURKA=modgene$AURKA, data=data, annot=annot, do.mapping=FALSE, do.scale=TRUE, plot=TRUE)
# dev.off()

# probabilities to belong to each of the 4 breast cancer molecular subtypes
# ER-/HER2-, HER2+, ER+/HER2- High Prolif, and ER+/HER2- Low Prolif
subtype.weights <- sbt.proba[[1]]
subtype.classif <- factor(sbt[,1], levels=c("ER-/HER2-", "HER2+", "ER+/HER2- High Prolif", "ER+/HER2- Low Prolif"))

# transform the continuous gene expressions in binary values (gene x > gene y)
# the set of binary values are referred to as barcode here after
# build an adjacency matrix that represent the pairs of genes one needs to consider
aa <- matrix(1, ncol=ncol(data), nrow=ncol(data), dimnames=list(colnames(data), colnames(data)))

# test with randomly selected 500 genes
# set.seed(1)
# aa[sample(1:100,replace = FALSE), sample(1:100,replace = FALSE)] <- 1

# binary barcode
binary_barcode <- compute.barcode(dataset=data, adjmatrix=aa)
# binary_barcode <- binary_barcode[!apply(binary_barcode==0,1,all),]  # Delete all rows having all zeros
# binary_barcode <- binary_barcode[!apply(binary_barcode==1,1,all),]  # Delete all rows having all ones

dim(binary_barcode)   #print dim of binary barcode

# probability barcode for each subtype
subtype_barcode <- do.call(cbind,lapply(1:ncol(subtype.weights), function(x) { compute.barcode.proba(dataset = data, theweight = subtype.weights[,x], adjmatrix = aa)}))
subtype_barcode <- subtype_barcode[rownames(subtype_barcode) %in% rownames(binary_barcode),]

dim(subtype_barcode) #print dim of subtype_barcode (number of rows should be the same as binary)

# subtype_barcode <- subtype_barcode[rownames(subtype_barcode) %in% rownames(binary_barcode),]

# feature selection
features <- feature.selection(proba = subtype_barcode, barcode = binary_barcode, n = 3, top = 5, redundant = 0, w = 1)

# standard vectors
standard <- standard.vectors(proba = subtype_barcode,features = features,n = 3, threshold = 0.5)

# jaccard distance
predicted <- jaccard.distance(proba = subtype_barcode, barcode = binary_barcode, features = features,n = 3,standard = standard)
# 
# # linear model to real expression values
# # change levels from "ER-/HER2-", "HER2+", "ER+/HER2- High Prolif", "ER+/HER2- Low Prolif" to #1-4
# for (i in 1:length(levels(subtype.classif))){    # use lapply here instead of for loop 
#     levels(subtype.classif)[i] <- i
# }
# 
# # accuracy
expected <- as.integer(lapply(1:ncol(binary_barcode),function(x){as.integer(subtype.classif[names(subtype.classif) == colnames(binary_barcode)[x]])}))
comp2 <- lapply(1:ncol(binary_barcode), function(x){ifelse(predicted[x] == expected[x],1,0)})
counts2 <- table(unlist(comp2))
accuracy2 <- (counts2[2]/(counts2[1]+counts2[2]))


# nfold <- 10
# nr <- ncol(binary_barcode)
# if(nfold > 1) k <- floor(nr/nfold) else {
#     k <- 1 
#     nfold <- nr
# }
# smpl <- sample(nr)
# 
# expected <- list()
# predicted <- list()
# accuracy2 <- list()
# counts2 <- list()
# for(i in 1:nfold){
#     if(i == nfold)
#         s.ix <- smpl[c(((i-1)*k+1):nr)]
#     else(i != nfold)
#         s.ix <- smpl[c(((i-1)*k+1):(i*k))]
#     train <- binary_barcode[,s.ix]
#     validation <- binary_barcode[,-s.ix]
#     features <- feature.selection(proba = subtype_barcode, barcode = train, n = 3, top = 5, redundant = 0, w = 1)
#     standard <- standard.vectors(proba = subtype_barcode,features = features,n = 3, threshold = 0.5)
#     predicted[[i]] <- jaccard.distance(proba = subtype_barcode, barcode = validation, features = features,n = 3,standard = standard)
#     expected[[i]] <- as.integer(lapply(1:ncol(validation),function(x){as.integer(subtype.classif[names(subtype.classif) == colnames(validation)[x]])}))
#     comp2 <- lapply(1:ncol(validation), function(x){ifelse(predicted[[i]][x] == expected[[i]][x],1,0)})
#     counts2 <- table(unlist(comp2))
#     accuracy2[[i]] <- (counts2[2]/(counts2[1]+counts2[2]))
# }

    
# foreach(i=1:10) %dopar% {
#     accuracy2 <- list()
#     train <- binary_barcode[,folds$subsets[folds$which != i]]
#     validation <- binary_barcode[,folds$subsets[folds$which == i]]
#     features <- feature.selection(proba = subtype_barcode, barcode = train, n = 3, top = 5, redundant = 1, w = 1)
#     standard <- standard.vectors(proba = subtype_barcode,features = features,n = 3, threshold = 0.5)
#     predicted <- jaccard.distance(proba = subtype_barcode, barcode = validation, features = features,n = 3,standard = standard)
#     expected <- as.integer(lapply(1:ncol(validation),function(x){as.integer(subtype.classif[names(subtype.classif) == colnames(validation)[x]])}))
#     comp2 <- lapply(1:ncol(validation), function(x){ifelse(predicted[x] == expected[x],1,0)})
#     counts2 <- table(unlist(comp2))
#     accuracy2 <- (counts2[2]/(counts2[1]+counts2[2]))}

# final <- lapply(1:10,function(x){diag(prop.table(table(predicted[[x]],expected[[x]]),2))})
# all.matrix <- abind(final,along = 2)
# themean <- apply(all.matrix,1,mean)*100
# names(themean) <- c("Basal-Like", "HER2-enriched", "Luminal A", "Luminal B")
# thesd <- melt(apply(all.matrix,1,sd)*100)
# df <- melt(themean)
# colnames(df) <- c("mean")
# df$sd <- thesd$value
# df$order <- c(4,3,1,2)
# df <- df[order(df$order),]
# ggplot(df,aes(x = factor(rownames(df)), y = mean, fill = factor(rownames(df))))+geom_bar(position = position_dodge(), stat = "identity") + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd))+ scale_fill_gray() + scale_x_discrete(limits=c("Luminal A","Luminal B","HER2-enriched", "Basal-Like")) + ylab("Accuracy (%)") + xlab("Breast Cancer Subtypes") + theme(legend.position = "none") + scale_y_continuous(limits=c(0,100))
# ggplot(df,aes(x = factor(rownames(df), levels=unique(as.character(rownames(df))) ), y = mean, fill = factor(rownames(df), levels=unique(as.character(rownames(df))))))+geom_bar(position = position_dodge(), stat = "identity") + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd)) + ylab("Accuracy (%)") + xlab("Breast Cancer Subtypes") + theme(legend.position = "none") + scale_y_continuous(limits=c(0,100)) + scale_fill_grey()
