SCMNET
======

SCMNET

Procedure
=========

1) Use METABRIC dataset as training set

2) Use the SCMGENE classification (see subtyping)

3) Consider only the 50 genes in PAM50, compute simple rules based on gene pairs (is gene i more expressed than gene j, binary value 0/1)

4) Feature selection, for each subtype s
  4.1) measure of discriminative value: for each gene pair, frequency of values 1 in subtype s - max(frequency of values 1 in subtype q \ s),the larger the value, the more discriminative
	4.2) measure of redundancy: Jaccard index between two gene pairs, the larger the index, the more redundant
	4.3) apply mRMR (minimum redundancy, maximum relevance) feature selection the select the top n (10, 20, 30, 40, 50, 100) gene pairs which are the most discriminative but the less redundant
	4.4) for each gene pair selected for subtype s, compute the frequency of values 1 in subtype s , if this frequency is > 0.5, the gene pair should ideally be equal to 1 in subtype s, 0 otherwise -> store this "typical" vector r for each subtype

5) Classification model:
	5.1) for a new sample/patient/tumor, compute the gene pairs and build the vector v of gene pairs selected in 4 for each subtype s
	5.2) compute the Jaccard index between v1 and r1, v2 and r2, v3 and r3, and v4 and r4 (assuming you have four subtypes as in our case)
	5.3) classify the new sample based on the maximum Jaccard index as computed in 5.2
