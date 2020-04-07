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


###############################################################################

######################     Run MV analysis for FUCOSYLATION    ###############################

###############################################################################


multi_fucosylation <- MultiSummary(sst, index=c(1, 3, 5, 7:10, 13:15, 17, 18, 22, 23), type = "outbred")
save(file='/mnt/polyomica/projects/MV_igg_replication/MANOVA/Rdata/20200407_multi_fucosylation_3K.Rdata', list=c("multi_fucosylation"))


#fucosylation

fucose <- multi_fucosylation$scan
fucose_coefs <-  multi_fucosylation$coef
rm(multi_fucosylation)

## snps sign for fucosylation in discovery

fucose_snps <- c("rs6964421", "rs540719", "rs7257072")

fucose.repl <- fucose[fucose$marker %in% fucose_snps,]

fucose.repl <- fucose.repl[,c(1:6)]
fwrite(fucose.repl, "./data/20200407_fucosylation_repl_3K.csv") 

## corresponding coefficients

ind <- which(fucose$marker %in%  fucose_snps)
fucose_coefs.repl <- fucose_coefs[ind,]

fucose_coefs.repl.df <- as.data.frame(fucose_coefs.repl)
a <- gsub("../3k_meta_23traits/meta_results/", "", names(fucose_coefs.repl.df))
names(fucose_coefs.repl.df) <- a
fucose_coefs.repl.df <- cbind(fucose[ind,"marker"],fucose_coefs.repl.df) 
fwrite(fucose_coefs.repl.df, "./data/20200407_fucosylation_coefs_repl.csv")
###############################################################################

######################     Run MV analysis for SYALYLATION    ###############################

###############################################################################

multi_sialylation <- MultiSummary(sst, index=c(15:18, 20:23), type = "outbred")
save(file='/mnt/polyomica/projects/MV_igg_replication/MANOVA/Rdata/20200407_multi_sialylation_3K.Rdata', list=c("multi_sialylation"))

### original sum stats list based on unmodified 3K GWASes

s <- multi_sialylation$scan
s_coefs <-  multi_sialylation$coef
rm(multi_sialylation)

## snps sign for galactosylation in discovery
s_snps <- c("rs2745851")

s.repl <-s[s$marker %in% s_snps,]

s.repl <- s.repl[,c(1:6)]
fwrite(s.repl, "./data/20200407_sialylation_repl_3K.csv") 

## corresponding coefficients

ind <- which(s$marker %in%  s_snps)
s_coefs.repl <- s_coefs[ind,]
s_coefs.repl.df <- as.data.frame(s_coefs.repl)
a <- gsub("../3k_meta_23traits/meta_results/", "", names(s_coefs.repl.df))
names(s_coefs.repl.df) <- a
s_coefs.repl.df <- cbind(s[ind,"marker"],s_coefs.repl.df) 
fwrite(s_coefs.repl.df, "./data/20200407_sialylation_coefs_repl.csv")
