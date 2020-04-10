### do table with additional results of mv for replicated snps ###

# extract info 

idnames_mv <- c("/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_N_glycosylation.Rdata",
				"/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_bisectingGlcNAc.Rdata",
				"/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_fucosylation.Rdata",
				"/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_sialylation.Rdata")


names_mv <- c("multi_N_glycosylation", "multi_bisectingGlcNAc", "multi_fucosylation", "multi_sialylation")

snps <- matrix(ncol=2, nrow=7) 


snps[,1] <- c("rs11895615", "rs1372288", "rs12635457", "rs479844", "rs4561508", "rs7257072", "rs2745851")
snps[,2] <- c("multi_bisectingGlcNAc", "multi_N_glycosylation", "multi_N_glycosylation", "multi_N_glycosylation", "multi_N_glycosylation", "multi_fucosylation", "multi_sialylation")

colnames(snps) <- c("SNP", "trait")
snps <- as.data.frame(snps)

extracted_data <- vector("list", nrow(snps))
names(extracted_data) <- paste(snps[,1], "scan", sep="_")


for (n in 1:4) {

	load(idnames_mv[n])
	eval(expr=parse(text=paste("MV_P=",names_mv[n])))
	x <- names_mv[n]
	eval(expr=parse(text=paste("rm(",names_mv[n],")")))

		for (i in (1:nrow(snps))) {
			
			if (x==snps$trait[i]) {

				extracted_data[[i]] <- MV_P$scan[(which(snps$SNP[i]==MV_P$scan[,1])),-c(2:6)]


				print(paste(n, i, sep=":")) 
					
			} 
		
		}

	print(n)	
	}


extracted_data_coef <- vector("list", nrow(snps))
names(extracted_data_coef) <- paste(snps[,1], "coef", sep="_")

for (n in 1:4) {

	load(idnames_mv[n])
	eval(expr=parse(text=paste("MV_P=",names_mv[n])))
	x <- names_mv[n]
	eval(expr=parse(text=paste("rm(",names_mv[n],")")))

	coef <- MV_P$coef
	rownames(coef) <- MV_P$scan[,1]

		for (i in (1:nrow(snps))) {
			
			if (x==snps$trait[i]) {

				extracted_data_coef[[i]] <- coef[(which(snps$SNP[i]==rownames(coef))),]


				print(paste(n, i, sep=":")) 
					
			} 
		
		}

	print(n)	
	}


# save results 

for (n in (1:7)) {

	path1 <- paste("/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/coefficients_and_other_info/", names(extracted_data)[n], ".txt", sep="")
	path2 <- paste("/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/coefficients_and_other_info/", names(extracted_data_coef)[n], ".txt", sep="")

	write.table(extracted_data[[n]], path1, col.names=T, row.names=F, quote=F)
	write.table(t(as.matrix(extracted_data_coef[[n]])), path2, col.names=T, row.names=F, quote=F)

	}