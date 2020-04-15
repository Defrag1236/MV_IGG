
library(data.table)
library(dplyr)
library(tidyr)

##### functions

## create linear combinations of gwases
GWAS_linear_combination=function(a,beta,se1,vary1=1,covm,N){
  vary2=sum(covm*(a%o%a))
  b=(beta%*%a)
  
  vary1_varg_n=se1^2+(beta[,1]^2)/N 
  
  se2=vary1_varg_n*(vary2/vary1)
  
  se2=se2-b^2/N
  
  se=sqrt(se2)
  
  b=b/sqrt(vary2)
  se=se/sqrt(vary2)
  
  out=cbind(b=b,se=se)
  colnames(out)=c("b","se")
  out=as.data.frame(out)
  return(out)
}




## load phenotypic correlations -> same for any trait
## obtain from the sst data
load('/mnt/polyomica/projects/MV_igg_replication/MANOVA/Rdata/20200306_SUM_STATs_3k_iggs.Rdata')

cormatrix <- sst$cor.pheno
a1 <- attributes(sst$cor.pheno)$dimnames[[2]]
a1 <- gsub("../3k_meta_23traits/meta_results/3K_meta_", "", a1)
attributes(sst$cor.pheno)$dimnames[[1]]<- a1
attributes(sst$cor.pheno)$dimnames[[2]]<- a1
cormatrix <- sst$cor.pheno


###############################################################################

########################    EXAMPLE using      ################################

#######################   N-glycosylation trait  ############################## 

###############################################################################

## assign trait name
comp_trait <- "Nglycosylation"

## which IGPs have to be selected to construct this trait
## for N-glycosylation that is all 23 IGPs
igps <- paste0("IGP", c(1:23))

  
  igps_csv <- paste0("3K_meta_", paste0(igps,".csv"))
  frames.names <- paste0("/mnt/polyomica/projects/MV_igg_replication/3k_meta_23traits/meta_results/", igps_csv)
  for (i in 1:length(frames.names)) assign(frames.names[[i]],fread(frames.names[[i]]))
  gwases <- lapply(frames.names, get)
  ## add to each gwas a column containing its id 
  names(gwases) <- paste0(igps,".csv")
  gwases<- lapply(seq_along(gwases), function(i) {gwases[[i]]$gwas_id <-names(gwases[i]) 
  gwases[[i]]})
  
  ## combine all snp effects into one dataframe
  betas.select <- function(mydata){
    betas <- subset(mydata, select = c( rs_id, beta_ma))
    gwas.id <- unique(mydata$gwas_id)
    colnames(betas)[2] <- paste("beta", gwas.id, sep ="_")
    return(betas)
  }
  
  list.betas  <- lapply(gwases, betas.select)
  
  tmp <- Reduce(merge, list.betas)
  
  ### add a column with se for the first GWAS
  se_gwas1 <- subset(gwases[[1]], select = c(rs_id, se_ma))
  
  tmp_1 <- merge(tmp, se_gwas1)
  
  colnames(tmp_1)[ncol(tmp_1)] <- "se_gwas1"
  
  ## add a column with minimal sample size across all GWASes for each snp 
  size.list <- lapply(gwases, subset, select = c(rs_id, n_ma))
  names(size.list) <- seq(1:length(size.list))
  size.list_2 <- lapply(seq_along(size.list), function(i) {colnames(size.list[[i]])[2] <- paste("n", names(size.list)[i]) 
  size.list[[i]]})
  size.list.m <-Reduce(merge, size.list_2)
  size.list.m$N_min <- do.call(pmin, size.list.m[,2:ncol(size.list.m)])
  
  betas.se.N <- merge(tmp_1, size.list.m[, c("rs_id", "N_min")])




### leave the correlations just between the gwases that comprise the trait
### for N-glycosylation that would be all IGPs 1:23

cor.matrix.filtered <- cormatrix[c(1:23), c(1:23)]

### correction betas need to be picked up from the discovery cohort
### each trait has an accordingly named Rdata file in 
### '/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/'

load('/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_N_glycosylation.Rdata')

coef <- as.data.frame(multi_N_glycosylation$coef)
### select the snps for which we want to replicate the association with the trait
### snps for N-glycosylation are listed below

snps <- c("rs1372288", "rs4561508", "rs12635457", "rs479844")

for (snp in snps) {
ind <- which(multi_N_glycosylation$scan$marker == snp)
coef.f <- coef[ind,grepl("est", names(coef))]
a1 <- as.numeric(coef.f)
beta.matrix <- as.matrix(betas.se.N[,2:(length(a1)+1)])

trait_lin_comb <- GWAS_linear_combination(a=a1, beta = beta.matrix, se1 = betas.se.N$se_gwas1, covm=cor.matrix.filtered , N = betas.se.N$N_min)


ind=which(betas.se.N$rs_id == snp)

res <- trait_lin_comb[ind, ]
library(dplyr)
library(tidyr)
res = mutate(res, Z=b/se,Pval=pchisq(Z^2,df=1,low=F))
res$rs_id <- snp
res$trait <- comp_trait

fwrite(res, file=paste0("./data/",paste(snp,paste0(comp_trait,"_linear_combination_3K.csv"), sep="_")))
}
