###############################################################################

#############     Multi-Traits Correlation Replication     ###################

###############################################################################

library(MultiABEL)

## load discovery gwases - object sst
load('/home/common/projects/Multivariate_analysis_IgG/2020_work_for_paper/results/Rdata/23_igg_sst.Rdata')
sst_d <- sst  

## check trait names
a <- attributes(sst_d$cor.pheno)$dimnames[[2]]


# [1] "IGP1.csv"  "IGP2.csv"  "IGP3.csv"  "IGP4.csv"  "IGP5.csv"  "IGP6.csv" 
# [7] "IGP7.csv"  "IGP8.csv"  "IGP9.csv"  "IGP10.csv" "IGP11.csv" "IGP12.csv"
# [13] "IGP13.csv" "IGP14.csv" "IGP15.csv" "IGP16.csv" "IGP17.csv" "IGP18.csv"
# [19] "IGP19.csv" "IGP20.csv" "IGP21.csv" "IGP22.csv" "IGP23.csv"


#load replication
load('/mnt/polyomica/projects/MV_igg_replication/MANOVA/Rdata/20200306_SUM_STATs_3k_iggs.Rdata')
sst_r <- sst

# rename replication traits to be consistent with discovery

a1 <- attributes(sst_r$cor.pheno)$dimnames[[2]]

# [1] "../3k_meta_23traits/meta_results/3K_meta_IGP1.csv" 
# [2] "../3k_meta_23traits/meta_results/3K_meta_IGP2.csv" 
# [3] "../3k_meta_23traits/meta_results/3K_meta_IGP3.csv" 
# [4] "../3k_meta_23traits/meta_results/3K_meta_IGP4.csv" 
# [5] "../3k_meta_23traits/meta_results/3K_meta_IGP5.csv" 
# [6] "../3k_meta_23traits/meta_results/3K_meta_IGP6.csv" 
# [7] "../3k_meta_23traits/meta_results/3K_meta_IGP7.csv" 
# [8] "../3k_meta_23traits/meta_results/3K_meta_IGP8.csv" 
# [9] "../3k_meta_23traits/meta_results/3K_meta_IGP9.csv" 
# [10] "../3k_meta_23traits/meta_results/3K_meta_IGP10.csv"
# [11] "../3k_meta_23traits/meta_results/3K_meta_IGP11.csv"
# [12] "../3k_meta_23traits/meta_results/3K_meta_IGP12.csv"
# [13] "../3k_meta_23traits/meta_results/3K_meta_IGP13.csv"
# [14] "../3k_meta_23traits/meta_results/3K_meta_IGP14.csv"
# [15] "../3k_meta_23traits/meta_results/3K_meta_IGP15.csv"
# [16] "../3k_meta_23traits/meta_results/3K_meta_IGP16.csv"
# [17] "../3k_meta_23traits/meta_results/3K_meta_IGP17.csv"
# [18] "../3k_meta_23traits/meta_results/3K_meta_IGP18.csv"
# [19] "../3k_meta_23traits/meta_results/3K_meta_IGP19.csv"
# [20] "../3k_meta_23traits/meta_results/3K_meta_IGP20.csv"
# [21] "../3k_meta_23traits/meta_results/3K_meta_IGP21.csv"
# [22] "../3k_meta_23traits/meta_results/3K_meta_IGP22.csv"
# [23] "../3k_meta_23traits/meta_results/3K_meta_IGP23.csv"

a1 <- gsub("../3k_meta_23traits/meta_results/3K_meta_", "", a1)
attributes(sst_r$cor.pheno)$dimnames[[1]]<- a1
attributes(sst_r$cor.pheno)$dimnames[[2]]<- a1

a2 <- attributes(sst_r$gwa)$dimnames[[2]]
a2.corr <-  gsub("../3k_meta_23traits/meta_results/3K_meta_", "", a2)
attributes(sst_r$gwa)$dimnames[[2]] <- a2.corr

### join with the alleles discovery
gwas_d <- as.data.frame(sst_d$gwa)
A1_d <- sst_d$alleles$A1 
A2_d <-  sst_d$alleles$A2 
gwas_d <- cbind(gwas_d, A1_d, A2_d)
colnames(gwas_d)[(dim(gwas_d)[2]-1):dim(gwas_d)[2]] <- c("A1", "A2") 



### join with the alleles replication
gwas_r <- as.data.frame(sst_r$gwa)
A1_r <- sst_r$alleles$A1 
A2_r <-  sst_r$alleles$A2 
gwas_r <- cbind(gwas_r, A1_r, A2_r)
colnames(gwas_r)[(dim(gwas_r)[2]-1):dim(gwas_r)[2]] <- c("A1", "A2") 

###############################################################################

#############     Multi-Traits Correlation Replication  PLOTS    ###################

###############################################################################


require(ggplot2)
require(cowplot)


#############     Multi-Traits Correlation Replication BISECTION     ###################



set.seed(2451)
traits <- paste0(paste0("IGP", c(5, 9, 10, 12, 14, 18, 23)), ".csv")

res_bisect <- MV.cor.test(marker = "rs11895615", gwa.1 = gwas_d, gwa.2 = gwas_r, R.1 = sst_d$cor.pheno,
            R.2 = sst_r$cor.pheno, traits = traits, plot=TRUE)

ci_bisect <- as.data.frame(res_bisect$res)

# $res
# correlation     ci.left    ci.right 
# 0.8095238   0.3333333   1.0000000 

pdf("./figs/Bisection_rs11895615.pdf", height=15, width=15)
df.plot <- res_bisect$df.plot
p1 <- ggplot()+ 
  geom_point(data=df.plot, mapping=aes(x=rank.1, y=rank.2, color=traits), size=2) + 
  geom_point(data=df.plot, mapping=aes(x=rank.1, y=rank.2, color=traits, size = se.beta), alpha = 0.2) + 
  stat_smooth(data=df.plot, mapping=aes(x=rank.1, y=rank.2), method = "lm", se=FALSE, color="black", size=0.3, fullrange = TRUE) + 
  coord_cartesian(xlim = c(0.5, 8), ylim = c(0.5, 8)) + xlim(0,200) + 
  scale_size_continuous(range = c(3, 10)) +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14,face="bold"), 
        strip.text.x = element_text(size = 16))+ 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.title.y=element_blank(), legend.position = c(0.8,0.3), 
        legend.background=element_rect(colour='NA', fill='transparent'), legend.key=element_blank(), 
        legend.title=element_text(size=14), 
        legend.text=element_text(size=12), legend.key.size = unit(1.4, 'lines')) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)), size = FALSE) +
  scale_colour_discrete(name = "Traits")

p2 <- ggplot(data=df.plot, aes(x=rank.1,y=mean.conc)) +
  coord_cartesian(xlim = c(0.5, 8), ylim = c(0, 18)) + 
  geom_bar(stat = "identity", aes(fill=traits), width = 0.4) + theme(legend.position="none") + theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) + geom_errorbar(aes(ymin = mean.conc - sd.conc,ymax = mean.conc + sd.conc), width = 0.1)  + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.title.y=element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

plot_grid(p1,p2,ncol=1,align = "v", rel_heights = c(2,1))
dev.off()

#############     Multi-Traits Correlation Replication N-glycosylation     ##################
traits <- paste0(paste0("IGP", c(1:23)), ".csv")

set.seed(86)

res_Nglyc_rs1372288 <- MV.cor.test(marker = "rs1372288", gwa.1 = gwas_d, gwa.2 = gwas_r, R.1 = sst_d$cor.pheno,
                          R.2 = sst_r$cor.pheno, traits = traits, plot=TRUE)

ci_Nglyc_rs1372288 <- as.data.frame(res_Nglyc_rs1372288$res)


pdf("./figs/Nglycosylation_rs1372288.pdf", height=15, width=15)
df.plot <- res_Nglyc_rs1372288$df.plot
p1 <- ggplot()+ 
  geom_point(data=df.plot, mapping=aes(x=rank.1, y=rank.2, color=traits), size=2) + 
  geom_point(data=df.plot, mapping=aes(x=rank.1, y=rank.2, color=traits, size = se.beta), alpha = 0.2) + 
  stat_smooth(data=df.plot, mapping=aes(x=rank.1, y=rank.2), method = "lm", se=FALSE, color="black", size=0.3, fullrange = TRUE) + 
  coord_cartesian(xlim = c(0.5, 24), ylim = c(0.5, 29)) + xlim(0,200) + 
  scale_size_continuous(range = c(3, 10)) +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14,face="bold"), 
        strip.text.x = element_text(size = 16))+ 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.title.y=element_blank(), legend.position = c(0.8,0.3), 
        legend.background=element_rect(colour='NA', fill='transparent'), legend.key=element_blank(), 
        legend.title=element_text(size=14), 
        legend.text=element_text(size=12), legend.key.size = unit(1.4, 'lines')) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)), size = FALSE) +
  scale_colour_discrete(name = "Traits")

p2 <- ggplot(data=df.plot, aes(x=rank.1,y=mean.conc)) +
  coord_cartesian(xlim = c(0.5, 24), ylim = c(0, 24)) + 
  geom_bar(stat = "identity", aes(fill=traits), width = 0.4) + theme(legend.position="none") + theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) + geom_errorbar(aes(ymin = mean.conc - sd.conc,ymax = mean.conc + sd.conc), width = 0.1)  + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.title.y=element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

plot_grid(p1,p2,ncol=1,align = "v", rel_heights = c(2,1))
dev.off()

### rs4561508
set.seed(736)

res_Nglyc_rs4561508 <- MV.cor.test(marker = "rs4561508", gwa.1 = gwas_d, gwa.2 = gwas_r, R.1 = sst_d$cor.pheno,
                                   R.2 = sst_r$cor.pheno, traits = traits, plot=TRUE)

ci_Nglyc_rs4561508 <- as.data.frame(res_Nglyc_rs4561508$res)


pdf("./figs/Nglycosylation_rs4561508.pdf", height=15, width=15)
df.plot <- res_Nglyc_rs4561508$df.plot
p1 <- ggplot()+ 
  geom_point(data=df.plot, mapping=aes(x=rank.1, y=rank.2, color=traits), size=2) + 
  geom_point(data=df.plot, mapping=aes(x=rank.1, y=rank.2, color=traits, size = se.beta), alpha = 0.2) + 
  stat_smooth(data=df.plot, mapping=aes(x=rank.1, y=rank.2), method = "lm", se=FALSE, color="black", size=0.3, fullrange = TRUE) + 
  coord_cartesian(xlim = c(0.5, 24), ylim = c(0.5, 29)) + xlim(0,200) + 
  scale_size_continuous(range = c(3, 10)) +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14,face="bold"), 
        strip.text.x = element_text(size = 16))+ 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.title.y=element_blank(), legend.position = c(0.8,0.3), 
        legend.background=element_rect(colour='NA', fill='transparent'), legend.key=element_blank(), 
        legend.title=element_text(size=14), 
        legend.text=element_text(size=12), legend.key.size = unit(1.4, 'lines')) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)), size = FALSE) +
  scale_colour_discrete(name = "Traits")

p2 <- ggplot(data=df.plot, aes(x=rank.1,y=mean.conc)) +
  coord_cartesian(xlim = c(0.5, 24), ylim = c(0, 24)) + 
  geom_bar(stat = "identity", aes(fill=traits), width = 0.4) + theme(legend.position="none") + theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) + geom_errorbar(aes(ymin = mean.conc - sd.conc,ymax = mean.conc + sd.conc), width = 0.1)  + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.title.y=element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

plot_grid(p1,p2,ncol=1,align = "v", rel_heights = c(2,1))
dev.off()

### rs12635457
set.seed(5923)

res_Nglyc_rs12635457 <- MV.cor.test(marker = "rs12635457", gwa.1 = gwas_d, gwa.2 = gwas_r, R.1 = sst_d$cor.pheno,
                                   R.2 = sst_r$cor.pheno, traits = traits, plot=TRUE)

ci_Nglyc_rs12635457 <- as.data.frame(res_Nglyc_rs12635457$res)


pdf("./figs/Nglycosylation_rs12635457.pdf", height=15, width=15)
df.plot <- res_Nglyc_rs12635457$df.plot
p1 <- ggplot()+ 
  geom_point(data=df.plot, mapping=aes(x=rank.1, y=rank.2, color=traits), size=2) + 
  geom_point(data=df.plot, mapping=aes(x=rank.1, y=rank.2, color=traits, size = se.beta), alpha = 0.2) + 
  stat_smooth(data=df.plot, mapping=aes(x=rank.1, y=rank.2), method = "lm", se=FALSE, color="black", size=0.3, fullrange = TRUE) + 
  coord_cartesian(xlim = c(0.5, 30), ylim = c(0.5, 24)) + xlim(0,200) + 
  scale_size_continuous(range = c(3, 10)) +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14,face="bold"), 
        strip.text.x = element_text(size = 16))+ 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.title.y=element_blank(), legend.position = c(0.8,0.3), 
        legend.background=element_rect(colour='NA', fill='transparent'), legend.key=element_blank(), 
        legend.title=element_text(size=14), 
        legend.text=element_text(size=12), legend.key.size = unit(1.4, 'lines')) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)), size = FALSE) +
  scale_colour_discrete(name = "Traits")

p2 <- ggplot(data=df.plot, aes(x=rank.1,y=mean.conc)) +
  coord_cartesian(xlim = c(0.5, 24), ylim = c(0, 24)) + 
  geom_bar(stat = "identity", aes(fill=traits), width = 0.4) + theme(legend.position="none") + theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) + geom_errorbar(aes(ymin = mean.conc - sd.conc,ymax = mean.conc + sd.conc), width = 0.1)  + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.title.y=element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

plot_grid(p1,p2,ncol=1,align = "v", rel_heights = c(2,1))
dev.off()


### rs479844
set.seed(253)

res_Nglyc_rs479844 <- MV.cor.test(marker = "rs479844", gwa.1 = gwas_d, gwa.2 = gwas_r, R.1 = sst_d$cor.pheno,
                                    R.2 = sst_r$cor.pheno, traits = traits, plot=TRUE)

ci_Nglyc_rs479844 <- as.data.frame(res_Nglyc_rs479844$res)


pdf("./figs/Nglycosylation_rs479844.pdf", height=15, width=15)
df.plot <- res_Nglyc_rs479844$df.plot
p1 <- ggplot()+ 
  geom_point(data=df.plot, mapping=aes(x=rank.1, y=rank.2, color=traits), size=2) + 
  geom_point(data=df.plot, mapping=aes(x=rank.1, y=rank.2, color=traits, size = se.beta), alpha = 0.2) + 
  stat_smooth(data=df.plot, mapping=aes(x=rank.1, y=rank.2), method = "lm", se=FALSE, color="black", size=0.3, fullrange = TRUE) + 
  coord_cartesian(xlim = c(0.5, 30), ylim = c(0.5, 24)) + xlim(0,200) + 
  scale_size_continuous(range = c(3, 10)) +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14,face="bold"), 
        strip.text.x = element_text(size = 16))+ 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.title.y=element_blank(), legend.position = c(0.8,0.3), 
        legend.background=element_rect(colour='NA', fill='transparent'), legend.key=element_blank(), 
        legend.title=element_text(size=14), 
        legend.text=element_text(size=12), legend.key.size = unit(1.4, 'lines')) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)), size = FALSE) +
  scale_colour_discrete(name = "Traits")

p2 <- ggplot(data=df.plot, aes(x=rank.1,y=mean.conc)) +
  coord_cartesian(xlim = c(0.5, 24), ylim = c(0, 24)) + 
  geom_bar(stat = "identity", aes(fill=traits), width = 0.4) + theme(legend.position="none") + theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) + geom_errorbar(aes(ymin = mean.conc - sd.conc,ymax = mean.conc + sd.conc), width = 0.1)  + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.title.y=element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

plot_grid(p1,p2,ncol=1,align = "v", rel_heights = c(2,1))
dev.off()


#############     Multi-Traits Correlation Replication Galactosylation     ###################



set.seed(49)
traits <- paste0(paste0("IGP", c(6:18, 20:23)), ".csv")

res_gal <- MV.cor.test(marker = "rs12049042", gwa.1 = gwas_d, gwa.2 = gwas_r, R.1 = sst_d$cor.pheno,
                          R.2 = sst_r$cor.pheno, traits = traits, plot=TRUE)

ci_gal <- as.data.frame(res_gal$res)


pdf("./figs/Galactosylation_rs12049042.pdf", height=15, width=15)
df.plot <- res_gal$df.plot
p1 <- ggplot()+ 
  geom_point(data=df.plot, mapping=aes(x=rank.1, y=rank.2, color=traits), size=2) + 
  geom_point(data=df.plot, mapping=aes(x=rank.1, y=rank.2, color=traits, size = se.beta), alpha = 0.2) + 
  stat_smooth(data=df.plot, mapping=aes(x=rank.1, y=rank.2), method = "lm", se=FALSE, color="black", size=0.3, fullrange = TRUE) + 
  coord_cartesian(xlim = c(0.5, 18), ylim = c(0.5, 18)) + xlim(0,200) + 
  scale_size_continuous(range = c(3, 10)) +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14,face="bold"), 
        strip.text.x = element_text(size = 16))+ 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.title.y=element_blank(), legend.position = c(0.8,0.3), 
        legend.background=element_rect(colour='NA', fill='transparent'), legend.key=element_blank(), 
        legend.title=element_text(size=14), 
        legend.text=element_text(size=12), legend.key.size = unit(1.4, 'lines')) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)), size = FALSE) +
  scale_colour_discrete(name = "Traits")

p2 <- ggplot(data=df.plot, aes(x=rank.1,y=mean.conc)) +
  coord_cartesian(xlim = c(0.5, 18), ylim = c(0, 18)) + 
  geom_bar(stat = "identity", aes(fill=traits), width = 0.4) + theme(legend.position="none") + theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) + geom_errorbar(aes(ymin = mean.conc - sd.conc,ymax = mean.conc + sd.conc), width = 0.1)  + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.title.y=element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

plot_grid(p1,p2,ncol=1,align = "v", rel_heights = c(2,1))
dev.off()


### correlations with 95% CIs
ci_list <- lapply(ls(pattern="ci_*"), get)

ci <- Reduce(cbind, ci_list)

ci$parameter <- rownames(ci)

library(tidyr)
library(dplyr)

## transpose
ci.t<- ci %>% 
  gather(trait, value, 1:6) %>%
  spread(parameter, value) 

fwrite(ci.t, file="./data/MT_correlation_results.csv") 
