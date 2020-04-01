##############################################################################################

##########################     Load summary statistics for the   #############################

##########################     replication cohort  (3K)          #############################

#############################################################################################

library(MultiABEL)

# Load list of indendent snps 
setwd("/home/common/projects/Multivariate_analysis_IgG/results/Rdata")
load ("20170206_indep_snps_for_IgGMV.RData")

# Load summary statistics for all 3K GWASES (because one of the traits is N-glycosylation total)
setwd("/mnt/polyomica/projects/MV_igg_replication/MANOVA/")
filenames_3k <- paste0(paste0("3K_meta_IGP", seq(1,23)), ".csv")
files.name <- c(paste0("../3k_meta_23traits/meta_results/", filenames_3k))

sst = load.summary(files=files.name,
                   cor.pheno = NULL, indep.snps = indep.snps, est.var = FALSE,
                   columnNames = c("rs_id", "ea", "ra", "eaf_ma", "beta_ma", "se_ma", "n_ma"), fixedN = NULL)
save(file='/mnt/polyomica/projects/MV_igg_replication/MANOVA/Rdata/20200306_SUM_STATs_3k_iggs.Rdata', list=c("sst"))

