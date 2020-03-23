### extract additional onfo for mv snps ###


# load data


load("/home/common/projects/Multivariate_analysis_IgG/results/Rdata/20170208_sst_23traits.Rdata")

mv_snps <-  read.table("/home/common/projects/Multivariate_analysis_IgG/results/clumping/clumping_merge/Clumping_9_mv_traits.txt", head=T, stringsAsFactors=F)



# make table for mv

idnames_mv <- c('/home/common/projects/Multivariate_analysis_IgG/results/Rdata/20170208_multi_N_glycosylation.Rdata',
	'/home/common/projects/Multivariate_analysis_IgG/results/Rdata/20170208_multi_monogalactosylation.Rdata',
 '/home/common/projects/Multivariate_analysis_IgG/results/Rdata/20170208_multi_digalactosylation.Rdata',
 '/home/common/projects/Multivariate_analysis_IgG/results/Rdata/20170208_multi_galactosylation.Rdata',
 '/home/common/projects/Multivariate_analysis_IgG/results/Rdata/20170208_multi_monosialylation.Rdata',
 '/home/common/projects/Multivariate_analysis_IgG/results/Rdata/20170208_multi_disialylation.Rdata',
 '/home/common/projects/Multivariate_analysis_IgG/results/Rdata/20170208_multi_sialylation.Rdata',
 '/home/common/projects/Multivariate_analysis_IgG/results/Rdata/20170208_multi_fucosylation.Rdata',
 '/home/common/projects/Multivariate_analysis_IgG/results/Rdata/20170208_multi_bisectingGlcNAc.Rdata')

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

	alleles_mv[n,] <- as.matrix(sst$alleles[(which(mv_snps$SNP[n]==rownames(sst$alleles))),1:2])

	}


mv_final_table <- cbind(mv_snps, extr_info[,c(2:3,5,6)], alleles_mv)
colnames(mv_final_table)[12:13] <- c("A1", "A2")

write.table(mv_final_table, "/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/MV_9_traits_clumping_with_beta_se_n_freq_alleles.txt", col.names=T, row.names=F, quote=F)