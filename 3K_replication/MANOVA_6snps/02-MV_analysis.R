###############################################################################

######################     Run MV analysis for N-glycosylation    ###############################

###############################################################################

library(MultiABEL)

## in the loading sst part all IGP1-23 were loaded in the right order - checked

# load sum stats
load('/mnt/polyomica/projects/MV_igg_replication/MANOVA/Rdata/20200306_SUM_STATs_3k_iggs.Rdata')

multi_N_glycosylation <- MultiSummary(sst, index=c(1:23), type = "outbred")
save(file='/mnt/polyomica/projects/MV_igg_replication/MANOVA/Rdata/20200325_multi_N_glycosylation_3K.Rdata', list=c("multi_N_glycosylation"))

#N-glycosylation

Ngly <- multi_N_glycosylation$scan
gly_coefs <-  multi_N_glycosylation$coef
rm(multi_N_glycosylation)

## snps sign for N-glycosylation in discovery

gly_snps <- c("rs1372288", "rs12635457", "rs479844", "rs4561508")

gly.repl <- Ngly[Ngly$marker %in% gly_snps,]

gly.repl <- gly.repl[,c(1:6)]
fwrite(gly.repl, "./data/20200325_N_glycosylation_repl_3K_CORRECTED.csv") 

## corresponding coefficients

ind <- which(Ngly$marker %in%  gly_snps)
gly_coefs.repl <- gly_coefs[ind,]

gly_coefs.repl.df <- as.data.frame(gly_coefs.repl)
a <- gsub("../3k_meta_23traits/meta_results/", "", names(gly_coefs.repl.df))
names(gly_coefs.repl.df) <- a
gly_coefs.repl.df <- cbind(Ngly[ind,"marker"],gly_coefs.repl.df) 
fwrite(gly_coefs.repl.df, "./data/20200325_N_glycosylation_coefs_repl.csv")


###############################################################################

######################     Run MV analysis for GALACTOSYLATION     ###############################

###############################################################################

## previously used 6:18 & 19 :23 -> 19 is an unidnentified structure, remove -> recalculate the trait
multi_N_galactosylation <- MultiSummary(sst, index=c(6:18, 20:23), type = "outbred")
save(file='/mnt/polyomica/projects/MV_igg_replication/MANOVA/Rdata/20200325_multi_N_galactosylation_3K.Rdata', list=c("multi_N_galactosylation"))

gal <- multi_N_galactosylation$scan
gal_coefs <-  multi_N_galactosylation$coef
rm(multi_N_galactosylation)

## snps sign for galactosylation in discovery
gal_snps <- c("rs12049042")

gal.repl <- gal[gal$marker %in% gal_snps,]

gal.repl <- gal.repl[,c(1:6)]
fwrite(gal.repl, "./data/20200325_galactosylation_repl_3K_CORRECTED.csv") 

## corresponding coefficients

ind <- which(gal$marker %in%  gal_snps)
gal_coefs.repl <- gal_coefs[ind,]
gal_coefs.repl.df <- as.data.frame(gal_coefs.repl)
a <- gsub("../3k_meta_23traits/meta_results/", "", names(gal_coefs.repl.df))
names(gal_coefs.repl.df) <- a
gal_coefs.repl.df <- cbind(gal[ind,"marker"],gal_coefs.repl.df) 
fwrite(gal_coefs.repl.df, "./data/20200325_galactosylation_coefs_repl.csv")

###############################################################################

######################     Run MV analysis for BISECTION    ###############################

###############################################################################

multi_N_bisection <- MultiSummary(sst, index=c(5, 9, 10, 12, 14, 18, 23), type = "outbred")
save(file='/mnt/polyomica/projects/MV_igg_replication/MANOVA/Rdata/20200325_multi_N_bisection_3K.Rdata', list=c("multi_N_bisection"))

### original sum stats list based on unmodified 3K GWASes

b <- multi_N_bisection$scan
b_coefs <-  multi_N_bisection$coef
rm(multi_N_bisection)

## snps sign for galactosylation in discovery
b_snps <- c("rs11895615")

b.repl <-b[b$marker %in% b_snps,]

b.repl <- b.repl[,c(1:6)]
fwrite(b.repl, "./data/20200325_bisection_repl_3K_CORRECTED.csv") 

## corresponding coefficients

ind <- which(b$marker %in%  b_snps)
b_coefs.repl <- b_coefs[ind,]
b_coefs.repl.df <- as.data.frame(b_coefs.repl)
a <- gsub("../3k_meta_23traits/meta_results/", "", names(b_coefs.repl.df))
names(b_coefs.repl.df) <- a
b_coefs.repl.df <- cbind(b[ind,"marker"],b_coefs.repl.df) 
fwrite(b_coefs.repl.df, "./data/20200325_bisection_coefs_repl.csv")



