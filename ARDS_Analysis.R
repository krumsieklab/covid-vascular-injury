#### Initialize ----

# clean workspace
rm(list = ls())
# set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(openxlsx)
library(dplyr)
library(magrittr)
library(maplet)

library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)

library(survival)
library(survminer)

# source helper functions
source("HelperFunctions.R")

# check if directory to save results exists, otherwise create
if(!dir.exists(paste0(Sys.Date(),sep=""))) {
  dir.create(paste0(Sys.Date(),sep=""))
}

#### Load Data ----

file_data <- "ARDS.xlsx"

D0 <- 
  # load proteomics data
  maplet::mt_load_xls(file=file_data,sheet="assay",samples_in_rows=F,id_col="feature_id") %>%
  # load sample annotations
  maplet::mt_anno_xls(file=file_data,sheet="clin",anno_type="samples",anno_id_col="sample_id",data_id_col="sample") %>%
  {.}
colnames(rowData(D0)) <- "ids"
D0 %<>%
  # load protein annotations
  maplet::mt_anno_xls(file=file_data,sheet="rowData",anno_type="features",anno_id_col="feature_id",data_id_col="ids") %>%
  # flag data as log-transformed
  maplet::mt_load_flag_logged()

#### load protein subsets ----

# load protein set defined in at risk cohort
prot <- load_protein_list(filename=sprintf("%s/Protlist.Rdata",Sys.Date()))

#### Figure 4A ----

D1 <- D0 %>%
  # scale data
  maplet::mt_pre_trans_scale() %>%
  # select samples that either were sampled within 10d or are not COVID
  maplet::mt_modify_filter_samples(filter= sampled_within_10d=="Yes" | vector != "COVID19") %>%
  # filter out COVID pneumonia
  maplet::mt_modify_filter_samples(filter= (vector=="COVID19" | disease == "ARDS") & disease != "PNA") %>%
  # subselect proteins
  maplet::mt_modify_filter_features(filter= Assay %in% prot) %>%
  {.}

# get proteomics data frame
X <- assay(D1)
colnames(X) <- D1$subject_id
rownames(X) <- rowData(D1)$Assay

# get clinical annotations
anno <- D1 %>% colData %>% as.data.frame
rownames(anno) <- colnames(X)
# define average protein abundance
anno$averageProt <- colSums(X[,match(rownames(anno),colnames(X))])/length(prot)

# define annotation colors
mypal <- brewer.pal(name = "Paired", n = 12)
ann_colors = list(
  Group = c(`Covid-19 ARDS`=mypal[8],`Bacterial sepsis ARDS`=mypal[12], `Influenza ARDS`=mypal[7]),
  Platelet = colorRampPalette(c("#F0EC19","#0AAD77"))(length(anno$Platelet)),
  Mortality = c(Yes="#F05349", No="#4577EA"),
  ARDS = c(`Acute Lung Injury`=brewer.pal(name = "Purples", n = 9)[2], Mild=brewer.pal(name = "Purples", n = 9)[3],
           Moderate=brewer.pal(name = "Purples", n = 9)[5], Severe=brewer.pal(name = "Purples", n = 9)[7]),
  LogANG2 = colorRampPalette(c("#ffcccc","#660066"))(length(anno$LogANG2)),
  LogRIPK3 = colorRampPalette(c("#d9d9d9","#0d0d0d"))(length(anno$LogRIPK3)),
  age = colorRampPalette(c("#a5b6d2","#0f2c64") )(nrow(anno)),
  averageProt = colorRampPalette(c("#f2f2f2","#404040"))(54)
)

# define heatmap breaks
breaksList = seq(-max(abs(min(assay(D1))), abs(max(assay(D1)))), 
                 max(abs(min(assay(D1))), abs(max(assay(D1)))), by=0.1)

# Modify ordering of the clusters using clustering callback option
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv, agglo.FUN = function(x){max(x)})
  as.hclust(dend)
}

# heatmap
out <- pheatmap(X,annotation_col = anno %>% dplyr::select(Group,Mortality,Platelet,
                                                          # ARDS,
                                                          LogANG2, age,averageProt),
                # cellwidth = 10, cellheight=10,
                clustering_method = "ward.D2",
                cutree_cols = 2, cutree_rows = 2,
                annotation_colors = ann_colors,
                clustering_callback = callback,
                border_color = NA,
                color = colorRampPalette(c("#4575B4","white","#D73027"))(length(breaksList)),
                breaks = breaksList)

# save heatmap to file
pdf(sprintf("%s/Figure_4A.pdf",Sys.Date()), width = 12, height = 8)
out
dev.off()

#### Compute Cluster differential analysis ----

# extract clusters and get annotations for differential analysis
z <- data.frame(Cluster=sort(cutree(out$tree_col, k=2))) %>% 
  tibble::rownames_to_column(var="subject_id") %>% 
  dplyr::left_join(anno) %>%
  dplyr::select(Cluster, subject_id, Group, Mortality, Platelet, 
                # ARDS,
                intubation_pf_ratio,LogRIPK3,LogANG2,age,averageProt) %>%
  dplyr::mutate(Cluster2=ifelse(Cluster %in% c(2),1,0)) %>%
  dplyr::mutate(death10=ifelse(Mortality == "Yes",1,0)) %>%
  dplyr::mutate(Cluster2=as.factor(Cluster2))

# perform differential tests
da <- data.frame(# Fisher's test for Death
  mortality = table(z %>% dplyr::select(death10, Cluster2)) %>% fisher.test(alternative = "greater") %>% .$p.value,
  # Wilcoxon test for continuous variables
  platelet = wilcox.test(x=z$Platelet[z$Cluster==1], y=z$Platelet[z$Cluster==2])$p.value,
  pf = wilcox.test(x=z$intubation_pf_ratio[z$Cluster==1], y=z$intubation_pf_ratio[z$Cluster==2])$p.value,
  RIPK3 = wilcox.test(x=z$LogRIPK3[z$Cluster==1], y=z$LogRIPK3[z$Cluster==2])$p.value,
  ANG2 = wilcox.test(x=z$LogANG2[z$Cluster==1], y=z$LogANG2[z$Cluster==2])$p.value,
  age = wilcox.test(x=z$age[z$Cluster==1], y=z$age[z$Cluster==2])$p.value,
  avProt = wilcox.test(x=z$averageProt[z$Cluster==1], y=z$averageProt[z$Cluster==2])$p.value)

#### Figure 4B -----

dt_surv <- z %>%
  dplyr::left_join(D1 %>% colData %>% as_data_frame %>% dplyr::select(subject_id, event,time_diff), 
                   by="subject_id") %>%
  dplyr::filter(time_diff>=0)

fit <- survfit(formula = Surv(time = time_diff, event = event) ~ Cluster, data = dt_surv)
pval <- survminer::surv_pvalue(fit = fit, data = dt_surv, method = "survdiff")
# Visualize with survminer
p_death_surv <- survminer::ggsurvplot(fit = fit, data = dt_surv, risk.table = T,
                                      title=sprintf("P=%.3f",pval$pval), conf.int = F,
                                      xlim = c(0, 60), break.time.by=10,
                                      xlab = "Time (days from blood draw)",
                                      palette = c("black","gray50")) 

# save plots to files
pdf(sprintf("%s/Figure_4B.pdf", Sys.Date()), width=5, height=6, onefile = F)
p_death_surv
dev.off()

#### Figure 4C ----

# plot differential log(ANGPT2) in the two clusters
p_ANG2 <- z %>%
  ggplot(aes(x=as.factor(Cluster), y=LogANG2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=Group),position=position_jitter(0.2)) +
  scale_colour_manual(values = ann_colors$Group) +
  theme_bw() +
  xlab("Cluster") +
  ylab("Log(ANGPT2)") +
  ggtitle(sprintf("Log(ANGPT2) p-value: %.2e", da$ANG2))

# save plots to files
pdf(sprintf("%s/Figure_4C.pdf", Sys.Date()), width=5, height=6, onefile = F)
p_ANG2
dev.off()

#### Figure 5 ----

# plot differential log(RIPK3) in the two clusters
p_RIPK3 <- z %>%
  ggplot(aes(x=as.factor(Cluster), y=LogRIPK3)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=Group),position=position_jitter(0.2)) +
  scale_colour_manual(values = ann_colors$Group) +
  theme_bw() +
  xlab("Cluster") +
  ylab("Log(RIP3K)") +
  ggtitle(sprintf("Log(RIP3K) p-value: %.2e", da$RIPK3))

# plot correlation between log(ANGPT2) and log(RIPK3) in the two clusters
zsub <- z %>%
  dplyr::mutate(Cluster=as.factor(Cluster))
cc <- cor.test(zsub$LogANG2,zsub$LogRIPK3,method = "pearson",alternative = "two.sided")
p_cor <- zsub %>%
  ggplot(aes(x=LogANG2, y=LogRIPK3, color=Cluster)) +
  geom_point(aes(color=Group)) +
  scale_colour_manual(values = ann_colors$Group) +
  geom_smooth(method='lm', formula= y~x, se = TRUE, size=0.5,
              aes(group=1), color = "black") +
  theme_bw() +
  ylim(c(min(z$LogRIPK3, na.rm=T), max(z$LogRIPK3, na.rm=T))) +
  xlab("Log(ANGPT2)") +
  ylab("Log(RIPK3)") +
  ggtitle(sprintf("Log(ANGPT2) vs Log(RIPK3): cor %.2f, p-value %.2e", cc$estimate, cc$p.value))

pdf(sprintf("%s/Figure_5.pdf", Sys.Date()), width=12, height=5, onefile = F)
ggarrange(p_RIPK3, p_cor,
          ncol=2, nrow=1,
          labels = "AUTO",
          common.legend = T)
dev.off()

##### Supplementary Figure 5 ----

# plot differential death in the two clusters
p_death <- z %>% 
  dplyr::select(Cluster, Mortality) %>% table() %>%
  as.data.frame %>%
  ggplot(aes(fill=Mortality, y=Freq, x=Cluster)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = ann_colors$Mortality) +
  theme_bw() +
  ylab("Number of patients") +
  ggtitle(sprintf("Mortality p-value: %.2e", da$mortality))

# plot differential P:F ratio in the two clusters
p_PF <- z %>%
  ggplot(aes(x=as.factor(Cluster), y=intubation_pf_ratio)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=Group),position=position_jitter(0.2)) +
  scale_colour_manual(values = ann_colors$Group) +
  theme_bw() +
  xlab("Cluster") +
  ylab("P:F Ratio at intubation") +
  ggtitle(sprintf("P:F Ratio p-value: %.2e", da$pf))

# plot differential Platelet in the two clusters
p_platelet <- z %>%
  ggplot(aes(x=as.factor(Cluster), y=Platelet)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=Group),position=position_jitter(0.2)) +
  scale_colour_manual(values = ann_colors$Group) +
  theme_bw() +
  xlab("Cluster") +
  ylab("Platelet") +
  ggtitle(sprintf("Platelet p-value: %.2e", da$platelet))

# plot differential age in the two clusters
p_age <- z %>%
  ggplot(aes(x=as.factor(Cluster), y=age)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=Group),position=position_jitter(0.2)) +
  scale_colour_manual(values = ann_colors$Group) +
  theme_bw() +
  xlab("Cluster") +
  ylab("age") +
  ggtitle(sprintf("Age p-value: %.2e", da$age))

# plot differential averageProtein in the two clusters
p_avProt <- z %>%
  ggplot(aes(x=as.factor(Cluster), y=averageProt)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=Group),position=position_jitter(0.2)) +
  scale_colour_manual(values = ann_colors$Group) +
  theme_bw() +
  xlab("Cluster") +
  ylab("averageProt") +
  ggtitle(sprintf("averageProt p-value: %.2e", da$avProt))


# save plots to files
pdf(sprintf("%s/SupplementaryFigure_5.pdf", Sys.Date()), width=15, height=5, onefile = F)
ggarrange(p_platelet,p_avProt, p_age,
          ncol=3,nrow=1,
          labels = "AUTO",
          common.legend = T)
dev.off()

#### Supplementary Figure 4 (Heatmap) ----

D2 <- D0 %>%
  # subselect only covid
  maplet::mt_modify_filter_samples(filter = vector=="COVID19") %>%
  # scale data
  maplet::mt_pre_trans_scale() %>%
  # select samples that either were sampled within 10d or are not COVID
  maplet::mt_modify_filter_samples(filter= sampled_within_10d == "Yes" | vector != "COVID19") %>%
  # filter out COVID pneumonia
  maplet::mt_modify_filter_samples(filter= (vector=="COVID19" | disease == "ARDS") & disease != "PNA") %>%
  # subselect proteins
  maplet::mt_modify_filter_features(filter= name %in% prot) %>%
  {.}

# get proteomics data frame
X2 <- assay(D2)
colnames(X2) <- D2$subject_id
rownames(X2) <- rowData(D2)$name

# get clinical annotations
anno2 <- D2 %>% colData %>% as.data.frame
rownames(anno2) <- colnames(X2)
# define average protein abundance
anno2$averageProt <- colSums(X2[,match(rownames(anno2),colnames(X2))])/length(prot)

# define annotation colors
mypal <- brewer.pal(name = "Paired", n = 12)
ann_colors = list(
  Group = c(`Covid-19 ARDS`=mypal[8],`Bacterial sepsis ARDS`=mypal[12], `Influenza ARDS`=mypal[7]),
  Platelet = colorRampPalette(c("#F0EC19","#0AAD77"))(length(anno2$Platelet)),
  Mortality = c(Yes="#F05349", No="#4577EA"),
  ARDS = c(`Acute Lung Injury`=brewer.pal(name = "Purples", n = 9)[2], Mild=brewer.pal(name = "Purples", n = 9)[3],
           Moderate=brewer.pal(name = "Purples", n = 9)[5], Severe=brewer.pal(name = "Purples", n = 9)[7]),
  LogANG2 = colorRampPalette(c("#ffcccc","#660066"))(length(anno2$LogANG2)),
  LogRIPK3 = colorRampPalette(c("#d9d9d9","#0d0d0d"))(length(anno2$LogRIPK3)),
  age = colorRampPalette(c("#a5b6d2","#0f2c64") )(nrow(anno2)),
  averageProt = colorRampPalette(c("#f2f2f2","#404040"))(54)
)

# define heatmap breaks
breaksList2 = seq(-max(abs(min(assay(D2))), abs(max(assay(D2)))), 
                  max(abs(min(assay(D2))), abs(max(assay(D2)))), by=0.1)

# Modify ordering of the clusters using clustering callback option
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv, agglo.FUN = function(x){max(x)})
  as.hclust(dend)
}

# heatmap
out2 <- pheatmap(X2,annotation_col = anno2 %>% dplyr::select(Group,Mortality,Platelet,
                                                             # ARDS,
                                                             LogANG2, age,averageProt),
                 # cellwidth = 10, cellheight=10,
                 clustering_method = "ward.D2",
                 cutree_cols = 2, cutree_rows = 2,
                 annotation_colors = ann_colors,
                 clustering_callback = callback,
                 border_color = NA,
                 color = colorRampPalette(c("#4575B4","white","#D73027"))(length(breaksList2)),
                 breaks = breaksList2)

# save heatmap to file
pdf(sprintf("%s/SupplementaryFigure_4A.pdf",Sys.Date()), width = 12, height = 8)
out2
dev.off()

#### Supplementary Figure 4 (Boxplots) ----

# extract clusters and get annotations for differential analysis
z2 <- data.frame(Cluster=sort(cutree(out2$tree_col, k=2))) %>% 
  tibble::rownames_to_column(var="subject_id") %>% 
  dplyr::left_join(anno2) %>%
  dplyr::select(Cluster, subject_id, Group, Mortality, Platelet, 
                # ARDS,
                intubation_pf_ratio,LogRIPK3,LogANG2,age,averageProt) %>%
  dplyr::mutate(Cluster2=ifelse(Cluster %in% c(2),1,0)) %>%
  dplyr::mutate(death10=ifelse(Mortality == "Yes",1,0)) %>%
  dplyr::mutate(Cluster2=as.factor(Cluster2))

# perform differential tests
da <- data.frame(# Fisher's test for Death
  mortality = table(z2 %>% dplyr::select(death10, Cluster2)) %>% fisher.test(alternative = "greater") %>% .$p.value,
  # Wilcoxon test for continuous variables
  platelet = wilcox.test(x=z2$Platelet[z2$Cluster==1], y=z2$Platelet[z2$Cluster==2])$p.value,
  pf = wilcox.test(x=z2$intubation_pf_ratio[z2$Cluster==1], y=z2$intubation_pf_ratio[z2$Cluster==2])$p.value,
  RIPK3 = wilcox.test(x=z2$LogRIPK3[z2$Cluster==1], y=z2$LogRIPK3[z2$Cluster==2])$p.value,
  ANG2 = wilcox.test(x=z2$LogANG2[z2$Cluster==1], y=z2$LogANG2[z2$Cluster==2])$p.value,
  age = wilcox.test(x=z2$age[z2$Cluster==1], y=z2$age[z2$Cluster==2])$p.value,
  avProt = wilcox.test(x=z2$averageProt[z2$Cluster==1], y=z2$averageProt[z2$Cluster==2])$p.value)

dt_surv2 <- z2 %>%
  dplyr::left_join(D2 %>% colData %>% as_data_frame %>% dplyr::select(subject_id, event,time_diff), 
                   by="subject_id") %>%
  dplyr::filter(time_diff>=0)

fit <- survfit(formula = Surv(time = time_diff, event = event) ~ Cluster, data = dt_surv2)
pval <- survminer::surv_pvalue(fit = fit, data = dt_surv2, method = "survdiff")
# Visualize with survminer
p_death_surv <- survminer::ggsurvplot(fit = fit, data = dt_surv2, risk.table = T,
                                      title=sprintf("P=%.3f",pval$pval), conf.int = F,
                                      xlim = c(0, 60), break.time.by=10,
                                      xlab = "Time (days from blood draw)",
                                      palette = c("black","gray50")) 

# plot differential death in the two clusters
p_death <- z2 %>% 
  dplyr::select(Cluster, Mortality) %>% 
  table() %>% as.data.frame() %>%
  ggplot(aes(fill=Mortality, y=Freq, x=Cluster)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = ann_colors$Mortality) +
  theme_bw() +
  ylab("Number of patients") +
  ggtitle(sprintf("Mortality p-value: %.2e", da$mortality))

# plot differential P:F ratio in the two clusters
p_PF <- z2 %>%
  ggplot(aes(x=as.factor(Cluster), y=intubation_pf_ratio)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=Group),position=position_jitter(0.2)) +
  scale_colour_manual(values = ann_colors$Group) +
  theme_bw() +
  xlab("Cluster") +
  ylab("P:F Ratio at intubation") +
  ggtitle(sprintf("P:F Ratio p-value: %.2e", da$pf))

# plot differential log(ANGPT2) in the two clusters
p_ANG2 <- z2 %>%
  ggplot(aes(x=as.factor(Cluster), y=LogANG2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=Group),position=position_jitter(0.2)) +
  scale_colour_manual(values = ann_colors$Group) +
  theme_bw() +
  xlab("Cluster") +
  ylab("Log(ANGPT2)") +
  ggtitle(sprintf("Log(ANGPT2) p-value: %.2e", da$ANG2))

# plot differential Platelet in the two clusters
p_platelet <- z2 %>%
  ggplot(aes(x=as.factor(Cluster), y=Platelet)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=Group),position=position_jitter(0.2)) +
  scale_colour_manual(values = ann_colors$Group) +
  theme_bw() +
  xlab("Cluster") +
  ylab("Platelet") +
  ggtitle(sprintf("Platelet p-value: %.2e", da$platelet))

# plot differential log(RIPK3) in the two clusters
p_RIPK3 <- z2 %>%
  ggplot(aes(x=as.factor(Cluster), y=LogRIPK3)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=Group),position=position_jitter(0.2)) +
  scale_colour_manual(values = ann_colors$Group) +
  theme_bw() +
  xlab("Cluster") +
  ylab("Log(RIP3K)") +
  ggtitle(sprintf("Log(RIP3K) p-value: %.2e", da$RIPK3))

# plot differential age in the two clusters
p_age <- z2 %>%
  ggplot(aes(x=as.factor(Cluster), y=age)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=Group),position=position_jitter(0.2)) +
  scale_colour_manual(values = ann_colors$Group) +
  theme_bw() +
  xlab("Cluster") +
  ylab("age") +
  ggtitle(sprintf("Age p-value: %.2e", da$age))

# plot differential averageProtein in the two clusters
p_avProt <- z2 %>%
  ggplot(aes(x=as.factor(Cluster), y=averageProt)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=Group),position=position_jitter(0.2)) +
  scale_colour_manual(values = ann_colors$Group) +
  theme_bw() +
  xlab("Cluster") +
  ylab("averageProt") +
  ggtitle(sprintf("averageProt p-value: %.2e", da$avProt))

# plot correlation between log(ANGPT2) and log(RIPK3) in the two clusters
zsub <- z2 %>%
  dplyr::mutate(Cluster=as.factor(Cluster))
cc2 <- cor.test(zsub$LogANG2,zsub$LogRIPK3,method = "pearson",alternative = "two.sided")
p_cor <- zsub %>%
  ggplot(aes(x=LogANG2, y=LogRIPK3, color=Cluster)) +
  geom_point(aes(color=Group)) +
  scale_colour_manual(values = ann_colors$Group) +
  geom_smooth(method='lm', formula= y~x, se = TRUE, size=0.5,
              aes(group=1), color = "black") +
  theme_bw() +
  ylim(c(min(z2$LogRIPK3, na.rm=T), max(z2$LogRIPK3, na.rm=T))) +
  xlab("Log(ANGPT2)") +
  ylab("Log(RIPK3)") +
  ggtitle(sprintf("Log(ANGPT2) vs Log(RIPK3): cor %.2f, p-value %.2e", cc2$estimate, cc2$p.value))

# save plots to files
pdf(sprintf("%s/SupplementaryFigure_4B.pdf", Sys.Date()), width=5, height=6, onefile = F)
p_death_surv
dev.off()

# save plots to files
pdf(sprintf("%s/SupplementaryFigure_4C-E.pdf", Sys.Date()), width=15, height=5, onefile = F)
ggarrange(p_platelet,p_age,p_ANG2, 
          ncol=3,nrow=1,
          labels = c("C","D","E"),
          common.legend = T)
dev.off()