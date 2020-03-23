###############################################################################

######################     Run clumping     ###################################

###############################################################################


# Function_for_shlop

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
# Cbind mv_results with m37 file (chr+pos) and do clumping

idnames_mv <- c('/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_N_glycosylation.Rdata',
	'/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_monogalactosylation.Rdata',
 '/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_digalactosylation.Rdata',
 '/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_galactosylation.Rdata',
 '/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_monosialylation.Rdata',
 '/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_disialylation.Rdata',
 '/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_sialylation.Rdata',
 '/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_fucosylation.Rdata',
 '/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_bisectingGlcNAc.Rdata')

names_mv <- c("multi_N_glycosylation", "multi_monogalactosylation", "multi_digalactosylation", "multi_galactosylation",
	"multi_monosialylation", "multi_disialylation", "multi_sialylation", "multi_fucosylation", "multi_bisectingGlcNAc")

names_for_save_mv <- c("multi_N_glycosylation.txt", "multi_monogalactosylation.txt", "multi_digalactosylation.txt", "multi_galactosylation.txt",
	"multi_monosialylation.txt", "multi_disialylation.txt", "multi_sialylation.txt", "multi_fucosylation.txt", "multi_bisectingGlcNAc.txt")

m37 <- read.table(file="/mnt/polyomica/projects/for_igg/data/IGP1.csv", head=T, stringsAsFactors=F, sep="\t")
m37 <- m37[,c(2,4:5)]

setwd("/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/clumping/delta_250kb/clumping_mv/")

for (n in 1:9) {

	load(idnames_mv[n])
	eval(expr=parse(text=paste("MV_P=",names_mv[n])))
	trait <- names_mv[n]
	eval(expr=parse(text=paste("rm(",names_mv[n],")")))

	y=MV_P$scan
	maf=pmin(1-y[,"freq"],y[,"freq"])
	ind=which(maf>=0.01 & y[,"n"]>3000)
	y=y[ind,]

	cycle_chr_pos <- y[(y[,1] %in% m37[,1]),]
	ind <- match(rownames(cycle_chr_pos),m37[,1])
	Chr <- m37[ind,2]
	Pos <- m37[ind,3]
	cycle_chr_pos <- cbind(cycle_chr_pos, Chr, Pos, trait)
	cycle_chr_pos <- cycle_chr_pos[,c("marker", "Chr", "Pos", "p", "trait")]
	
	clumped <- function_for_shlop_28_12_2017(cycle_chr_pos,p_value="p",pos="Pos",snp="marker",delta=2.5e5,chr="Chr",thr=5e-8/(9+21),trait="trait")


	write.table(file=names_for_save_mv[n],x=clumped,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
	}

# Cbind uv_results with m37 file (chr+pos) and do clumping

idnames_uv <- NULL
names_uv <- NULL
names_for_save_uv <- NULL



for (n in c(1:77)) {

    idnames_uv <-c(idnames_uv, paste ("IGP", n, ".csv", sep=""))

    names_uv <- c(names_uv, paste("IGP", n, sep=""))

    names_for_save_uv <- c(names_for_save_uv, paste("IGP", n, ".txt", sep=""))

    }

library(data.table)
setwd("/mnt/polyomica/projects/for_igg/data/")
for (n in 1:77) {

    y <- fread(idnames_uv[n], head=T, stringsAsFactors=F,data.table=F)
    trait <- names_uv[n]

    maf=pmin(1-y[,"eaf"],y[,"eaf"])
    ind=which(maf>=0.01 & y[,"n"]>3000)
    y=y[ind,]
    
    cycle_chr_pos <- y[(y[,2] %in% m37[,1]),]
    ind <- match(cycle_chr_pos[,2],m37[,1])
    Chr <- m37[ind,2]
    Pos <- m37[ind,3]
    cycle_chr_pos <- cbind(cycle_chr_pos, Chr, Pos, trait)
    cycle_chr_pos <- cycle_chr_pos[,c("rs_id", "Chr", "Pos", "p", "trait")]
    colnames(cycle_chr_pos) <- c("marker", "Chr", "Pos", "p", "trait")

    clumped <- function_for_shlop_28_12_2017(cycle_chr_pos, p_value="p", pos="Pos", snp="marker", delta=2.5e5, chr="Chr", thr=5e-8/21, trait=NULL)


    z <- paste("/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/clumping/delta_250kb/clumping_uv/", names_for_save_uv[n], sep="")
    write.table(file=z, x=clumped,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
    }