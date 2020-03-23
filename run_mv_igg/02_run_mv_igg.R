###############################################################################

######################     Run MV analysis     ###############################

###############################################################################

library(MultiABEL)
load('/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/23_igg_sst.Rdata')

# Run MV for each group of traits and save results in .Rdata

	multi_N_glycosylation <- MultiSummary(sst, index=c(1:23), type = "outbred")
	save(file='/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_N_glycosylation.Rdata', list=c("multi_N_glycosylation"))
	rm(multi_N_glycosylation)

	multi_monogalactosylation <- MultiSummary (sst, index=c(6:10, 15), type = "outbred")
	save(file='/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_monogalactosylation.Rdata', list=c("multi_monogalactosylation"))
	rm(multi_monogalactosylation)

	multi_digalactosylation <- MultiSummary (sst, index=c(11:14, 16:18, 20:23), type = "outbred")
	save(file='/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_digalactosylation.Rdata', list=c("multi_digalactosylation"))
	rm(multi_digalactosylation)

	multi_galactosylation <- MultiSummary (sst, index=c(6:18, 20:23), type = "outbred")
	save(file='/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_galactosylation.Rdata', list=c("multi_galactosylation"))
	rm(multi_galactosylation)

	multi_monosialylation <- MultiSummary (sst, index=c(15:18), type = "outbred")
	save(file='/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_monosialylation.Rdata', list=c("multi_monosialylation"))
	rm(multi_monosialylation)

	multi_disialylation <- MultiSummary (sst, index=c(20:23), type = "outbred")
	save(file='/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_disialylation.Rdata', list=c("multi_disialylation"))
	rm(multi_disialylation)

	multi_sialylation <- MultiSummary (sst, index=c(15:18, 20:23), type = "outbred")
	save(file='/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_sialylation.Rdata', list=c("multi_sialylation"))
	rm(multi_sialylation)

	multi_fucosylation <- MultiSummary (sst, index=c(1, 3, 5, 7:10, 13:15, 17, 18, 22, 23), type = "outbred")
	save(file='/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_fucosylation.Rdata', list=c("multi_fucosylation"))
	rm(multi_fucosylation)

	multi_bisectingGlcNAc <- MultiSummary (sst, index=c(5, 9, 10, 12, 14, 18, 23), type = "outbred")
	save(file='/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/multi_bisectingGlcNAc.Rdata', list=c("multi_bisectingGlcNAc"))
	rm(multi_bisectingGlcNAc)