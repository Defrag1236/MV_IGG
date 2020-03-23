### extraqct gwas information ###

# load data

load("/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_N_glycosylation.Rdata")
load("/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_galactosylation.Rdata")
load("/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_bisectingGlcNAc.Rdata")

igp59 <- read.table(file="/mnt/polyomica/projects/for_igg/data/IGP59.csv", head=T, stringsAsFactors=F, sep="\t")
igp14 <- read.table(file="/mnt/polyomica/projects/for_igg/data/IGP14.csv", head=T, stringsAsFactors=F, sep="\t")


# filter data

multi_N_clear <- multi_N_glycosylation$scan[,c(1,4)]
multi_gal_clear <- multi_galactosylation$scan[,c(1,4)]
multi_bis_clear <- multi_bisectingGlcNAc$scan[,c(1,4)]
igp59_clear <- igp59[,c(2,12)]
igp14_clear <- igp14[,c(2,12)]

colnames(igp59_clear) <- colnames(multi_N_clear)
colnames(igp14_clear) <- colnames(multi_N_clear)

# save data

write.table(multi_N_clear, "/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/for_zoom_locus/raw_gwas/multi_N_glycosylation.txt", row.names=F, quote=F)
write.table(multi_gal_clear, "/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/for_zoom_locus/raw_gwas/multi_galactosylation.txt", row.names=F, quote=F)
write.table(multi_bis_clear, "/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/for_zoom_locus/raw_gwas/multi_bisectingGlcNAc.txt", row.names=F, quote=F)
write.table(igp59_clear, "/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/for_zoom_locus/raw_gwas/IGP59.txt", row.names=F, quote=F)
write.table(igp14_clear, "/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/for_zoom_locus/raw_gwas/IGP14.txt", row.names=F, quote=F)