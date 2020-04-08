### do qq plots for mv traits ###

# do plots

mv_result <- list.files("/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/")
mv_result <- mv_result[-c(1,2)]

library(qqman)

setwd("/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/")

for (n in (1:9)) {

	load(mv_result[n])
	name_mv <- gsub(".Rdata", "", mv_result[n])
	eval(expr=parse(text=paste("MV_P=",name_mv)))
	eval(expr=parse(text=paste("rm(",name_mv,")")))

	name_for_save <- paste("/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/plots/", name_mv, ".jpg", sep="")
	title_for_plot <- paste ("Q-Q plot of ", name_mv, " p-values", sep="")

	jpeg(name_for_save) 

	qq(MV_P$scan$p, main = title_for_plot, col = "blue4")

	dev.off()

	}

