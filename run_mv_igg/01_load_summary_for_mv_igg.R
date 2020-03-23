###############################################################################

##########################     Load summary      ##############################

###############################################################################

library(MultiABEL)

# Load list of indendent snps 

setwd("/home/common/projects/Multivariate_analysis_IgG/results/Rdata")
load ("20170206_indep_snps_for_IgGMV.RData")

# Load summary statistics  

setwd("/mnt/polyomica/projects/for_igg/data/")


files.name <- list.files()[-c(1)]
files.name <- files.name[c(1,12,17:23, 2:11, 13:16)]

sst = load.summary(files=files.name,
	cor.pheno = NULL, indep.snps = indep.snps, est.var = FALSE,
	columnNames = c("rs_id", "ea", "ra", "eaf", "beta", "se", "n"), fixedN = NULL)

save(file='/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/23_igg_sst.Rdata', list=c("sst"))