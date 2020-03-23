###############################################################################

######################     Run clumping for all tables     ####################

###############################################################################

# Merge UV traits in one big table 

setwd("/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/clumping/clumping_uv")
lf <- list.files()
out <- read.table(lf[1], head=T, stringsAsFactors=F)
cn <- c("SNP","Chr","Pos","P","trait")
colnames(out) <- cn

 for (f in lf[-1]) {
    x <- read.table(f, head=T, stringsAsFactors=F)
      if (!(is.na(x[1,1]))) {
      x <- subset (x, x[,4]<5e-8/21, select=1:5)
      }
        if (!(is.na(x[1,1]))) {
        colnames(x)=cn
        out <- rbind(out, x)
        }
  }

write.table(file="/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/clumping/general_clumping/UV_traits_77.txt", x=out, row.names=F, col.names=T, quote=FALSE, sep="\t")

# Merge MV traits in one big table 

setwd("/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/clumping/clumping_mv")
lf1 <- list.files()
out1 <- read.table(lf1[1], head=T, stringsAsFactors=F)
out1 <- out1[,1:5]
cn <- c("SNP","Chr","Pos","P","trait")
colnames(out1) <- cn

 for (f in lf1[-1]) {
    y <- read.table(f, head=T, stringsAsFactors=F)
      if (!(is.na(y[1,1]))) {
      y <- subset (y, y[,4]<5e-8/30, select=1:5)
      }
        if (!(is.na(y[1,1]))) {
        colnames(y)=cn
        out1 <- rbind(out1, y)
        }
  }

write.table(file="/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/clumping/general_clumping/MV_traits_9.txt", x=out1, row.names=F, col.names=T, quote=FALSE, sep="\t")

# Do clumping for UV and MV and for general table

function_for_shlop_28_12_2017=function(locus_table,p_value="P",pos="POS",snp="SNP",
                                       delta=2.5e5,chr="CHR",thr=5e-8,trait=NULL){
 #locus_table=bt
    locus_table[,p_value]=as.numeric(locus_table[,p_value])
    
    if (!is.null(trait)){
        traits="traits"
        locus_table=cbind(locus_table,traits=locus_table[,trait])
        locus_table[,traits]=as.character(locus_table[,traits])
    }
    
    out=locus_table[0,]
    
    locus_table=locus_table[locus_table[,p_value]<=thr,]
    
    i=1
    if (nrow(locus_table)>0){
        locus_table[,pos]=as.numeric(locus_table[,pos])
        locus_table[,p_value]=as.numeric(locus_table[,p_value])
        Zx <-locus_table
        Zx=Zx[order(Zx[,p_value]),]
        #n_traits=1
        #Zx=cbind(Zx,n_traits)
        i=1
        while (nrow(Zx)>0){       
            ind=which((abs(Zx[i,pos]-Zx[,pos])<=delta)&(Zx[i,chr]==Zx[,chr]))
                      
            if (!is.null(trait)){
                Zx[i,traits]=paste(unique(Zx[ind,trait]),collapse = ";")
            }
            
            out=rbind(out,Zx[i,])
            Zx=Zx[-ind,]            
        }
        rownames(out)=as.character(out[,snp])
    }
    
    if (!is.null(trait)){
        j=1
        out=cbind(out,Ntraits=1)
        out[,"Ntraits"]=as.numeric(out[,"Ntraits"])
        for (j in 1:nrow(out)){
            trs=unique(unlist(strsplit(out[j,traits],split = ";")))
            out[j,traits]=paste(trs,collapse = ";")
            out[j,"Ntraits"]=length(trs)
        }
    }
    
 return(out)
}
general <- rbind (out,out1)
write.table(file="/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/clumping/general_clumping/All_traits_77+9.txt", x=general, row.names=F, col.names=T, quote=FALSE, sep="\t")


clump_uv <- function_for_shlop_28_12_2017(out, p_value="P", pos="Pos", snp="SNP", delta=5e5, chr="Chr", thr=5e-8/21, trait="trait")
clump_mv <- function_for_shlop_28_12_2017(out1, p_value="P", pos="Pos", snp="SNP", delta=5e5, chr="Chr", thr=5e-8/30, trait="trait")
clump_g <- function_for_shlop_28_12_2017(general, p_value="P", pos="Pos", snp="SNP", delta=5e5, chr="Chr", thr=5e-8/21, trait="trait")


write.table(file="/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/clumping/general_clumping/Clumping_77_uv_traits.txt", x=clump_uv, row.names=F, col.names=T, quote=FALSE, sep="\t")
write.table(file="/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/clumping/general_clumping/Clumping_9_mv_traits.txt", x=clump_mv, row.names=F, col.names=T, quote=FALSE, sep="\t")
write.table(file="/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/clumping/general_clumping/Clumping_77+9_traits.txt", x=clump_g, row.names=F, col.names=T, quote=FALSE, sep="\t")
