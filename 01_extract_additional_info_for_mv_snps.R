### extract additional onfo for mv snps ###


# load data

mv_snps <-  read.table("/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/clumping/general_clumping/Clumping_9_mv_traits.txt", head=T, stringsAsFactors=F)

library(data.table)

IGP1 <- fread(file="/mnt/polyomica/projects/for_igg/data/IGP1.csv", head=T, stringsAsFactors=F, sep="\t", data.table=F)



# make table for mv

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


extracted_data <- list()

for (n in 1:9) {

	load(idnames_mv[n])
	eval(expr=parse(text=paste("MV_P=",names_mv[n])))
	x <- names_mv[n]
	eval(expr=parse(text=paste("rm(",names_mv[n],")")))

		for (i in (1:nrow(mv_snps))) {
			
			if (x==mv_snps$trait[i]) {

				extracted_data[[i]] <- MV_P$scan[(which(mv_snps$SNP[i]==MV_P$scan[,1])),1:6]


				print(paste(n, i, sep=":")) 
					
			} 
		
		}

	print(n)
	}

extr_info <-  do.call("rbind", extracted_data)


alleles_mv <- matrix(ncol=2, nrow=nrow(mv_snps))

for ( n in (1:nrow(mv_snps))) {

	alleles_mv[n,] <- as.matrix(IGP1[(which(mv_snps$SNP[n]==IGP1[,2])),6:7])

	}


mv_final_table <- cbind(mv_snps, extr_info[,c(2:3,5,6)], alleles_mv)


write.table(mv_final_table, "/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/clumping/general_clumping/MV_9_traits_clumping_with_beta_se_n_freq_alleles.txt", col.names=T, row.names=F, quote=F)