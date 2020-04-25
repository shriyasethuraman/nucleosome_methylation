setwd("/home/shriyas/cifs-lab/Shriya/Samples/manuscript_prep/nucleosome_methylation/working/methyl_mutants_DMRs/All_final_DMRs/drm2_DMRs/")
source('~/cifs-lab/RIP_manuscript/Revised_Figures/Fig7/Fig_7A/hp-profile-heatmap_modGB.r')
source('~/cifs-lab/Shriya/Codes/window_generator_mod.R')
library(boot)
library(ggplot2)


setwd("/home/shriyas/cifs-lab/Shriya/Samples/manuscript_prep/nucleosome_methylation/working/methyl_mutants_DMRs/All_final_DMRs/drm2_DMRs/")

table1_1 <- read.csv("col0_MNase_drm2_CHH_DMR_counts_2.csv", header = F, na.strings = "NA")
table2_1 <- read.csv("drm2_MNase_drm2_CHH_DMR_counts_2.csv", header = F, na.strings = "NA")

values1 <- ((table2_1[,-1]+0.1)/0.25)/((table1_1[,-1]+0.1)/0.99)
output1 <- window_generator(values1, 50, 1, 1)

# png("ratio-mean drm2 over col0 MNase 500bases around drm2-CHH DMRs.png", width=1000, height=600)
# matplot(colMeans(values1, na.rm=TRUE),col="red", type="l", lty=c(1,1,1,1,2,2,2,2), ylab="drm2/col0 RPM",lwd=2, main="ratio-mean drm2/col0 drm2 CHH-DMRs")
# abline(v=500)
# abline(v=600)
# dev.off()
# 
# 
# minim<-apply(table1_1[,-1]/0.99, 1, min)
# maxim<-apply(table1_1[,-1]/0.99, 1, max)
# x1 <-((table1_1[,-1]/0.99)-minim)/(maxim-minim)
# minim<-apply(table2_1[,-1]/0.25, 1, min)
# maxim<-apply(table2_1[,-1]/0.25, 1, max)
# x2<-((table2_1[,-1]/0.25)-minim)/(maxim-minim)
# 
# values1 <- ((x2[,-1]+0.1))/((x1[,-1]+0.1))
# 
# png("max_min ratio-mean drm2 over col0 MNase 500bases around drm2-CHH DMRs.png", width=1000, height=600)
# matplot(colMeans(values1, na.rm=TRUE),col="red", type="l", lty=c(1,1,1,1,2,2,2,2), ylab="drm2/col0 RPM",lwd=2, main="max_min ratio-mean drm2/col0 rep2 drm2 CHH-DMRs")
# abline(v=500)
# abline(v=600)
# dev.off()
# 
# 
# mean1 <- colMeans(x1[,-1], na.rm=T)
# mean2 <- colMeans(x2[,-1], na.rm=T)
# values1 <- (mean2+0.1)/(mean1+0.1)
# png("max_min mean-ratio rep1 drm2 over col0 MNase 500bases around drm2-CHH DMRs.png", width=1000, height=600)
# matplot(values1,col="red", type="l", lty=c(1,1,1,1,2,2,2,2), ylab="drm2/col0 RPM",lwd=2, main="max_min mean-ratio rep1 drm2/col0 drm2 CHH-DMRs")
# abline(v=500)
# abline(v=600)
# dev.off()
# 
# #values1 <- calculate_values("col0_MNase_drm2_CHH_DMR_counts_2.csv","drm2_MNase_drm2_CHH_DMR_counts_2.csv",0.99,0.25,2)


tab1 <- output1
ci_output <- matrix(0, nrow=(dim(tab1)[2])-2, ncol=3)
my.mean = function(x, index){ return(mean(x[index], na.rm = TRUE))}
for(i in seq(2,(dim(tab1)[2])-1)){
  b <- boot(tab1[,i+1], my.mean, R=1000)
  #plot(b)
  ci_output[i-1,1] <- (boot.ci(b, type = "norm"))$normal[2]
  ci_output[i-1,2] <- (boot.ci(b, type = "norm"))$normal[3]
  ci_output[i-1,3] <- (boot.ci(b, type = "norm"))$t0
}

colnames(ci_output) <- c("down","up","mean")

png("bootstrap 50bp_win_sl_1bp drm2 over col0 MNase 500bases around drm2-CHH DMRs.png", width=1000, height=600)
pdf("bootstrap 50bp_win_sl_1bp drm2 over col0 MNase 500bases around drm2-CHH DMRs.pdf", width=10, height=6)
ggplot(data.frame(ci_output), aes(1:1046)) + 
  xlab("position") + 
  ylab("drm2/col0 MNase") +
  geom_line(aes(y=mean), colour="blue") + 
  geom_ribbon(aes(ymin=down, ymax=up), alpha=0.5) +
  geom_vline(xintercept = 472) +
  geom_vline(xintercept = 572) +
  ggtitle("bootstrapped 50bp_win_sl_1bp drm2/col0 MNase profile at drm2 CHH-DMRs and calculated the means and std deviations for 1000permutations")
dev.off()




table1_1 <- read.csv("col0_rep2_MNase_drm2_d25_CHH_DMR_counts_2.csv", header = F, na.strings = "NA")
table2_1 <- read.csv("drm2_rep2_MNase_drm2_d25_CHH_DMR_counts_2.csv", header = F, na.strings = "NA")

values1 <- ((table2_1[,-1]+0.1)/3.12)/((table1_1[,-1]+0.1)/1.42)

output1 <- window_generator(values1, 50, 1, 1)
matplot(cbind(colMeans(output1, na.rm=TRUE)),col="red", type="l", lty=c(1,1,1,1,2,2,2,2), ylab="drm2/col0 RPM",lwd=2, main="drm2/col0 drm2 CHH-DMRs")
abline(v=450)
abline(v=550)

tab1 <- output1
ci_output <- matrix(0, nrow=(dim(tab1)[2])-2, ncol=3)
my.mean = function(x, index){ return(mean(x[index], na.rm = TRUE))}
for(i in seq(2,(dim(tab1)[2])-1)){
  b <- boot(tab1[,i+1], my.mean, R=1000)
  #plot(b)
  ci_output[i-1,1] <- (boot.ci(b, type = "norm"))$normal[2]
  ci_output[i-1,2] <- (boot.ci(b, type = "norm"))$normal[3]
  ci_output[i-1,3] <- (boot.ci(b, type = "norm"))$t0
}

colnames(ci_output) <- c("down","up","mean")

png("bootstrap 50bp_win_sl_1bp rep2 drm2 over col0 MNase 500bases around drm2-CHH DMRs.png", width=1000, height=600)
pdf("bootstrap 50bp_win_sl_1bp rep2 drm2 over col0 MNase 500bases around drm2-CHH DMRs.pdf", width=10, height=6)
ggplot(data.frame(ci_output), aes(1:1046)) + 
  xlab("position") + 
  ylab("drm2/col0 MNase") +
  geom_line(aes(y=mean), colour="blue") + 
  geom_ribbon(aes(ymin=down, ymax=up), alpha=0.5) +
  geom_vline(xintercept = 472) +
  geom_vline(xintercept = 572) +
  ggtitle("bootstrapped 50bp_win_sl_1bp rep2 drm2/col0 MNase profile at drm2 CHH-DMRs and calculated the means and std deviations for 1000permutations")
dev.off()




table1_1 <- read.csv("col0_MNase_drm2_CHH_DMR_counts_2.csv", header = F, na.strings = "NA")
table2_1 <- read.csv("drm2_MNase_drm2_CHH_DMR_counts_2.csv", header = F, na.strings = "NA")

table1_2 <- read.csv("col0_rep2_MNase_drm2_d25_CHH_DMR_counts_2.csv", header = F, na.strings = "NA")
table2_2 <- read.csv("drm2_rep2_MNase_drm2_d25_CHH_DMR_counts_2.csv", header = F, na.strings = "NA")

table1 <- table1_1[,-1] + table1_2[,-1]
table2 <- table2_1[,-1] + table2_2[,-1]
values1 <- ((table2+0.1)/3.37)/((table1+0.1)/2.41)

output1 <- window_generator(values1, 50, 1, 1)
matplot(cbind(colMeans(output1, na.rm=TRUE)),col="red", type="l", lty=c(1,1,1,1,2,2,2,2), ylab="drm2/col0 RPM",lwd=2, main="drm2/col0 drm2 CHH-DMRs")
abline(v=450)
abline(v=550)

tab1 <- output1
ci_output <- matrix(0, nrow=(dim(tab1)[2])-2, ncol=3)
my.mean = function(x, index){ return(mean(x[index], na.rm = TRUE))}
for(i in seq(2,(dim(tab1)[2])-1)){
  b <- boot(tab1[,i+1], my.mean, R=1000)
  #plot(b)
  ci_output[i-1,1] <- (boot.ci(b, type = "norm"))$normal[2]
  ci_output[i-1,2] <- (boot.ci(b, type = "norm"))$normal[3]
  ci_output[i-1,3] <- (boot.ci(b, type = "norm"))$t0
}

colnames(ci_output) <- c("down","up","mean")

png("bootstrap 50bp_win_sl_1bp rep1+2 drm2 over col0 MNase 500bases around drm2-CHH DMRs.png", width=1000, height=600)
pdf("bootstrap 50bp_win_sl_1bp rep1+2 drm2 over col0 MNase 500bases around drm2-CHH DMRs.pdf", width=10, height=6)
ggplot(data.frame(ci_output), aes(1:1046)) + 
  xlab("position") + 
  ylab("drm2/col0 MNase") +
  geom_line(aes(y=mean), colour="blue") + 
  geom_ribbon(aes(ymin=down, ymax=up), alpha=0.5) +
  geom_vline(xintercept = 472) +
  geom_vline(xintercept = 572) +
  ggtitle("bootstrapped 50bp_win_sl_1bp rep1+2 drm2/col0 MNase profile at drm2 CHH-DMRs and calculated the means and std deviations for 1000permutations")
dev.off()





setwd("/home/shriyas/cifs-lab/Shriya/Samples/manuscript_prep/nucleosome_methylation/working/methyl_mutants_DMRs/All_final_DMRs/drm2_DMRs/")

table1 <- read.csv("WT_CHH_dataset.csv", header = F)
table2 <- read.csv("drm2_CHH_dataset.csv", header = F)

xx <- (table2[,-1]-table1[,-1])

output1 = window_generator(xx, 10, 1)

tab1 <- output1
ci_output <- matrix(0, nrow=(dim(tab1)[2])-2, ncol=3)
my.mean = function(x, index){ return(mean(x[index], na.rm = TRUE))}
for(i in seq(2,(dim(tab1)[2])-1)){
  b <- boot(tab1[,i+1], my.mean, R=1000)
  #plot(b)
  ci_output[i-1,1] <- (boot.ci(b, type = "norm"))$normal[2]
  ci_output[i-1,2] <- (boot.ci(b, type = "norm"))$normal[3]
  ci_output[i-1,3] <- (boot.ci(b, type = "norm"))$t0
}

colnames(ci_output) <- c("down","up","mean")

png("bootstrap 10bp-1bp-slide drm2-col0 CHHmethylation at nucleotide 500bases around drm2-CHH DMRs.png", width=1000, height=600)
ggplot(data.frame(ci_output), aes(1:1086)) + 
  xlab("position") +
  ylab("Average CHHmethylation of all permutations of nucleosomes") +
  geom_line(aes(y=mean), colour="blue") + 
  geom_ribbon(aes(ymin=down, ymax=up), alpha=0.5) +
  geom_vline(xintercept = 495) +
  geom_vline(xintercept = 595) +
  ggtitle("bootstrapped drm2-col0 CHHme 10bp-window-1bp-slide around drm2-CHH DMRs and calculated the means and std dev for 1000permutations")
dev.off()

output1 = window_generator(xx, 50, 1)

tab1 <- output1
ci_output <- matrix(0, nrow=(dim(tab1)[2])-2, ncol=3)
my.mean = function(x, index){ return(mean(x[index], na.rm = TRUE))}
for(i in seq(2,(dim(tab1)[2])-1)){
  b <- boot(tab1[,i+1], my.mean, R=1000)
  #plot(b)
  ci_output[i-1,1] <- (boot.ci(b, type = "norm"))$normal[2]
  ci_output[i-1,2] <- (boot.ci(b, type = "norm"))$normal[3]
  ci_output[i-1,3] <- (boot.ci(b, type = "norm"))$t0
}

colnames(ci_output) <- c("down","up","mean")

png("bootstrap 50bp-1bp-slide drm2-col0 CHHmethylation at nucleotide 500bases around drm2-CHH DMRs.png", width=1000, height=600)
pdf("bootstrap 50bp-1bp-slide drm2-col0 CHHmethylation at nucleotide 500bases around drm2-CHH DMRs.pdf", width=10, height=6)
ggplot(data.frame(ci_output), aes(1:1046)) + 
  xlab("position") +
  ylab("Average CHHmethylation of all permutations of nucleosomes") +
  geom_line(aes(y=mean), colour="blue") + 
  geom_ribbon(aes(ymin=down, ymax=up), alpha=0.5) +
  geom_vline(xintercept = 475) +
  geom_vline(xintercept = 575) +
  ggtitle("bootstrapped drm2-col0 CHHme 50bp-window-1bp-slide around drm2-CHH DMRs and calculated the means and std dev for 1000permutations")
dev.off()

