### extract additional info for general clumping data ###

# load data

clumped <- read.table("/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/clumping/general_clumping/Clumping_77+9_traits.txt", head=T, stringsAsFactors=F)

library(data.table)

IGP2 <- fread(file="/mnt/polyomica/projects/for_igg/data/IGP2.csv", head=T, stringsAsFactors=F, sep="\t", data.table=F)
IGP29 <- fread(file="/mnt/polyomica/projects/for_igg/data/IGP29.csv", head=T, stringsAsFactors=F, sep="\t", data.table=F)
IGP45 <- fread(file="/mnt/polyomica/projects/for_igg/data/IGP45.csv", head=T, stringsAsFactors=F, sep="\t", data.table=F)
IGP59,  <- fread(file="/mnt/polyomica/projects/for_igg/data/IGP59.csv", head=T, stringsAsFactors=F, sep="\t", data.table=F)
IGP66 <- fread(file="/mnt/polyomica/projects/for_igg/data/IGP66.csv", head=T, stringsAsFactors=F, sep="\t", data.table=F)
IGP70 <- fread(file="/mnt/polyomica/projects/for_igg/data/IGP70.csv", head=T, stringsAsFactors=F, sep="\t", data.table=F)
IGP74 <- fread(file="/mnt/polyomica/projects/for_igg/data/IGP74.csv", head=T, stringsAsFactors=F, sep="\t", data.table=F)
IGP77 <- fread(file="/mnt/polyomica/projects/for_igg/data/IGP77.csv", head=T, stringsAsFactors=F, sep="\t", data.table=F)


# extract info

idnames_mv <- c("/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_N_glycosylation.Rdata",
				"/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_galactosylation.Rdata",
				"/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_bisectingGlcNAc.Rdata",
				"/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_sialylation.Rdata",
				"/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_fucosylation.Rdata")


names_mv <- c("multi_N_glycosylation", "multi_galactosylation", "multi_sialylation", "multi_fucosylation", "multi_bisectingGlcNAc")


extracted_data <- list()

for (n in 1:5) {

	load(idnames_mv[n])
	eval(expr=parse(text=paste("MV_P=",names_mv[n])))
	x <- names_mv[n]
	eval(expr=parse(text=paste("rm(",names_mv[n],")")))

		for (i in (1:nrow(clumped))) {
			
			if (x==clumped$trait[i]) {

				extracted_data[[i]] <- MV_P$scan[(which(clumped$SNP[i]==MV_P$scan[,1])),c(1:3, 5:6)]


				print(paste(n, i, sep=":")) 
					
			} 
		
		}

	print(n)
	}


igp_list <- list(IGP2, IGP29, IGP45, IGP59, IGP66, IGP70, IGP74, IGP77)
names(igp_list) <- c("IGP2", "IGP29", "IGP45", "IGP59", "IGP66", "IGP70", "IGP74", "IGP77")



for (n in (1:8)) {

	x <- names(igp_list)[n]
	y <- igp_list[[n]]

		for (i in (1:nrow(clumped))) {
			
			if (x==clumped$trait[i]) {

				extracted_data[[i]] <- y[(which(clumped$SNP[i]==y[,2])),c(2,8,13,10:11)]
				colnames(extracted_data[[i]]) <- colnames(extracted_data[[1]]) 

				print(paste(n, i, sep=":")) 
					
			} 

		}
}


extr_info <-  do.call("rbind", extracted_data)

alleles <- matrix(ncol=2, nrow=nrow(clumped))

for ( n in (1:nrow(clumped))) {

	alleles[n,] <- as.matrix(igp_list[[1]][(which(clumped$SNP[n]==igp_list[[1]][,2])),6:7])

	}


final_table <- cbind(clumped, extr_info[,c(2:5)], alleles)
colnames(final_table)[12:13] <- c("ea", "ra")

write.table(final_table, "/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/general_clumping_result_with_beta_se_n_freq_alleles.txt", col.names=T, row.names=F, quote=F)