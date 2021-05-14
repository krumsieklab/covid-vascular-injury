#### Initialize ----

# clean workspace
rm(list = ls())
# set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load libraries
library(dplyr)
library(tidyverse)
library(magrittr)
library(maplet)

library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(pheatmap)

# source helper functions
source("HelperFunctions.R")

# check if directory to save results exists, otherwise create
if(!dir.exists(paste0(Sys.Date(),sep=""))) {
  dir.create(paste0(Sys.Date(),sep=""))
}

#### Load Data ----

file_data <- "Recovery.xlsx"

D <- 
  # load proteomics data
  maplet::mt_load_xls(file=file_data,sheet="data",samples_in_rows=F) %>%
  # load sample annotations
  maplet::mt_anno_xls(file=file_data,sheet="sampleanno",anno_type="samples",anno_id_col="sample",data_id_col="sample") %>%
  # load protein annotations
  maplet::mt_anno_xls(file=file_data,sheet="proteinanno",anno_type="features",anno_id_col="OlinkID",data_id_col="name") %>%
  # flag data as log-transformed
  maplet::mt_load_flag_logged()

#### Load protein set ----

# load protein set defined in at risk cohort
prot <- load_protein_list(filename=sprintf("%s/Protlist.Rdata",Sys.Date()))

#### Prepare df for plotting ----

D %<>%
  # remove third sample for one patient
  maplet::mt_modify_filter_samples(filter=Label %in% c("ICU","postICU")) %>%
  # scale data
  maplet::mt_pre_trans_scale() %>%
  # filter only protein set
  maplet::mt_modify_filter_features(filter= Assay %in% prot)

dt <- D %>% assay %>% as.data.frame %>% 
  # substitute protein names to OlinkIDs
  tibble::rownames_to_column(var="OlinkID") %>%
  dplyr::left_join(D %>% rowData %>% as.data.frame %>% dplyr::select(OlinkID,Assay)) %>%
  dplyr::select(-OlinkID) %>%
  tibble::column_to_rownames(var="Assay") %>%
  # transpose
  t %>% as.data.frame %>%
  # only select protein set
  dplyr::select(any_of(prot)) %>%
  # left join sample annotations
  tibble::rownames_to_column(var="SAMPLE_NAME") %>%
  dplyr::left_join(D %>% colData %>% as.data.frame) 

dt0 <- dt[dt$Label=="ICU",]
dt1 <- dt[dt$Label=="postICU",]
# make sure two dataset are in the same order
dt1 <- dt1[match(dt0$Patient_ID,dt1$Patient_ID),]

#### Figure 6A ----

anno <- dt %>% 
  dplyr::select(SAMPLE_NAME, Label, Patient_ID, death, age, platelet, LogANG2, FollowUp_30d, FollowUp_12m) %>%
  # dplyr::mutate(outcome=ifelse(!is.na(outcome),outcome,"unknown")) %>%
  tibble::column_to_rownames(var="SAMPLE_NAME")
# define average protein abundance
X <- dt %>% tibble::column_to_rownames("SAMPLE_NAME") %>% dplyr::select(any_of(prot))
anno$averageProt <- rowSums(X[match(rownames(anno),rownames(X)),])/length(prot)

ann_colors <- list(
  Label = c(ICU="white",postICU="black"),
  death = c(Yes="#F05349", No="#4577EA"),
  platelet = colorRampPalette(c("#F0EC19","#0AAD77"))(nrow(anno)),
  LogANG2 = colorRampPalette(c("#ffcccc","#660066"))(length(anno$LogANG2)),
  age = colorRampPalette(c("#a5b6d2","#0f2c64") )(nrow(anno)),
  # outcome = c(good="cyan4",bad="tomato4",unknown="gray"),
  averageProt = colorRampPalette(c("#f2f2f2","#404040"))(nrow(anno)),
  FollowUp_30d = colorRampPalette(c("wheat","firebrick"))(10),
  FollowUp_12m = colorRampPalette(c("wheat","firebrick"))(10)
)

# define heatmap breaks
breaksList = seq(-max(abs(min(assay(D))), abs(max(assay(D)))), 
                 max(abs(min(assay(D))), abs(max(assay(D)))), by=0.1)

# Modify ordering of the clusters using clustering callback option
# 1-max(x)+1/min(x)
callback4 = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv, agglo.FUN = function(x){1-max(x)+1/min(x)}) #1/(max(x)+3)+1/min(x)+mean(x)
  as.hclust(dend)
}

n_clust <- 2
hh1 <- pheatmap(dt %>% dplyr::filter(Label=="postICU") %>% 
                  tibble::column_to_rownames(var="SAMPLE_NAME") %>% 
                  dplyr::select(any_of(prot)) %>% 
                  t %>% as.data.frame, 
                annotation_col = anno %>% dplyr::filter(Label=="postICU") %>%
                  dplyr::select(platelet, LogANG2, age, FollowUp_12m, averageProt),
                annotation_colors = ann_colors,
                # cellwidth = 10, cellheight=10,
                clustering_method = "ward.D2", cutree_cols = n_clust, cutree_rows = 2,
                labels_col = anno %>% dplyr::filter(Label=="postICU") %>% dplyr::pull(Patient_ID),
                clustering_callback = callback4,
                border_color = NA,
                color = colorRampPalette(c("#4575B4","white","#D73027"))(length(breaksList)),
                breaks = breaksList)

pdf(sprintf("%s/Figure_6A.pdf", Sys.Date()), width = 9, height = 6, onefile = F)
print(hh1)
dev.off()

#### Cluster Differential Analysis ----

# extract clusters and get annotations for differential analysis
z1 <- data.frame(Cluster=sort(cutree(hh1$tree_col, k=n_clust))) %>% 
  tibble::rownames_to_column(var="SAMPLE_NAME") %>% 
  dplyr::left_join(anno %>% tibble::rownames_to_column("SAMPLE_NAME"), 
                   by="SAMPLE_NAME") %>%
  dplyr::mutate(color=ifelse(Cluster==1,"tomato","cadetblue"))

sortvar <- z1 %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarize(var=mean(averageProt)) %>%
  dplyr::arrange(var) %>%
  dplyr::mutate(Cluster_ord=rank(-var, ties.method = "min"))

z1 %<>% 
  dplyr::left_join(sortvar,by="Cluster") %>%
  dplyr::arrange(Cluster_ord)

# perform differential tests
da <- data.frame(
  # kendall correlation for follow-up scores
  fu_12m = cor.test(x=z1$FollowUp_12m, y=z1$Cluster_ord,method = "kendall", use="pairwise.complete.obs")$p.value,
  # Wilcoxon test for continuous variables
  platelet = wilcox.test(x=z1$platelet[z1$Cluster_ord==1], y=z1$platelet[z1$Cluster_ord==2])$p.value,
  LogANG2 = wilcox.test(x=z1$LogANG2[z1$Cluster_ord==1], y=z1$LogANG2[z1$Cluster_ord==2])$p.value,
  age = wilcox.test(x=z1$age[z1$Cluster_ord==1], y=z1$age[z1$Cluster_ord==2])$p.value,
  avProt = wilcox.test(x=z1$averageProt[z1$Cluster_ord==1], y=z1$averageProt[z1$Cluster_ord==2])$p.value)

#### Figure 6B ----

# plot follow-up
p_fu <- z1 %>%
  ggplot(aes(x=as.factor(Cluster_ord), y=FollowUp_12m)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  theme_bw() +
  xlab("Cluster") +
  ylab("Follow-Up Score (12 months)") +
  ggtitle(sprintf("Follow-up p-value: %.2e", da$fu_12m))

# save plots to files
pdf(sprintf("%s/Figure_6B.pdf", Sys.Date()), width=5, height=5, onefile = F)
p_fu
dev.off()

#### Figure 6C ----

dt_diff <- dt1 %>% dplyr::select(any_of(prot)) %>%
  subtract(dt0 %>% dplyr::select(any_of(prot))) %>%
  dplyr::mutate(name=dt0$Patient_ID) %>%
  tibble::column_to_rownames(var="name") %>%
  tibble::rownames_to_column("Patient_ID") %>%
  # use post ICU clustering
  dplyr::left_join(z1 %>% dplyr::select(Patient_ID,Cluster, color),
                   by="Patient_ID") %>%
  dplyr::mutate(Cluster=as.factor(Cluster)) %>%
  dplyr::left_join(anno %>% dplyr::select(Patient_ID,FollowUp_12m) %>% distinct(), by="Patient_ID") 

res_diff <- lapply(prot %>% {names(.)=.;.}, function(x){
  dat <- dt_diff %>%
    dplyr::mutate(var=!!sym(x)) %>%
    # dplyr::filter(Cluster!=3) %>%
    {.}
  fit <- lm(formula = var ~ Cluster , data = dat)
  list(protein=x, 
       statistic=summary(fit)$coefficients[2,"Estimate"],
       p.value=summary(fit)$coefficients[2,"Pr(>|t|)"])
}) %>% {do.call(rbind,.)} %>% as.data.frame %>%
  dplyr::mutate(protein=as.character(protein)) %>%
  dplyr::mutate(p.value=unlist(p.value)) %>%
  dplyr::mutate(p.adj=p.adjust(p.value,method = "BH")) %>%
  dplyr::arrange(p.value) 

# create boxplots
p <- lapply(res_diff$protein %>% {names(.)=.;.}, function(x){
  dt %>%
    ggplot(aes(x=Timepoint,y=!!sym(x))) +
    geom_boxplot(outlier.shape = NA) +
    geom_point() +
    geom_line(aes(group = Patient_ID, color=Patient_ID),size=1) +
    scale_color_manual(values = setNames(z1$color,z1$Patient_ID)) +
    theme_bw() +
    ggtitle(sprintf("%s adjusted p-value:%.2e",x,res_diff %>% filter(protein==x) %>% dplyr::pull(p.adj)))
})

# save boxplots to file
ncol=5
pdf(sprintf("%s/Figure_6C.pdf", Sys.Date()), width = 5*ncol, height = 5*ceiling(length(prot)/5), onefile = F)
ggarrange(plotlist = p,
          ncol=ncol, nrow=ceiling(length(prot)/ncol),
          common.legend = T)
dev.off()

#### Supplementary Figure 6 A-D -----

# plot differential platelet in the two clusters
p_platelet <- z1 %>%
  ggplot(aes(x=as.factor(Cluster_ord), y=platelet)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  theme_bw() +
  xlab("Cluster") +
  ylab("Platelet") +
  ggtitle(sprintf("Platelet p-value: %.2e", da$platelet))

# plot differential ANG2 in the two clusters
p_ang2 <- z1 %>%
  ggplot(aes(x=as.factor(Cluster_ord), y=LogANG2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  theme_bw() +
  xlab("Cluster") +
  ylab("Log(ANG2)") +
  ggtitle(sprintf("LogANG2 p-value: %.2e", da$LogANG2))

# plot differential age in the two clusters
p_age <- z1 %>%
  ggplot(aes(x=as.factor(Cluster_ord), y=age)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  theme_bw() +
  xlab("Cluster") +
  ylab("Age") +
  ggtitle(sprintf("Age p-value: %.2e", da$age))

# plot differential averageProt in the two clusters
p_avProt<- z1 %>%
  ggplot(aes(x=as.factor(Cluster_ord), y=averageProt)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  theme_bw() +
  xlab("Cluster") +
  ylab("Average Protein") +
  ggtitle(sprintf("Average Protein p-value: %.2e", da$avProt))

# save plots to files
pdf(sprintf("%s/SupplementaryFigure_6A-D.pdf", Sys.Date()), width=20, height=5, onefile = F)
ggarrange(p_platelet,p_age,p_avProt,p_ang2,
          ncol=4,nrow=1,
          common.legend = T)
dev.off()

#### Supplementary Figure 6 E-F ----

platelet_diff <- dt1 %>% dplyr::select(platelet) %>%
  subtract(dt0 %>% dplyr::select(platelet)) %>%
  dplyr::mutate(name=dt0$Patient_ID) %>%
  tibble::column_to_rownames(var="name") %>%
  tibble::rownames_to_column("Patient_ID") %>%
  # use post ICU clustering
  dplyr::left_join(z1 %>% dplyr::select(Patient_ID,Cluster, color),
                   by="Patient_ID") %>%
  dplyr::mutate(Cluster=as.factor(Cluster)) %>%
  dplyr::left_join(anno %>% dplyr::select(Patient_ID,FollowUp_12m) %>% distinct(), by="Patient_ID") 

res_pl <- lapply("platelet" %>% {names(.)=.;.}, function(x){
  dat <- platelet_diff %>%
    dplyr::mutate(var=!!sym(x)) %>%
    # dplyr::filter(Cluster!=3) %>%
    {.}
  fit <- lm(formula = var ~ Cluster , data = dat)
  list(statistic=summary(fit)$coefficients[2,"Estimate"],p.value=summary(fit)$coefficients[2,"Pr(>|t|)"])
}) %>% {do.call(rbind,.)} %>% as.data.frame

# create boxplots
p_pl <- dt %>%
    ggplot(aes(x=Timepoint,y=platelet)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point() +
    geom_line(aes(group = Patient_ID, color=Patient_ID),size=1) +
    scale_color_manual(values = setNames(z1$color,z1$Patient_ID)) +
    theme_bw() +
    ggtitle(sprintf("Platelet p-value:%.2e",res_pl$p.value))

ang2_diff <- dt1 %>% dplyr::select(LogANG2) %>%
  subtract(dt0 %>% dplyr::select(LogANG2)) %>%
  dplyr::mutate(name=dt0$Patient_ID) %>%
  tibble::column_to_rownames(var="name") %>%
  tibble::rownames_to_column("Patient_ID") %>%
  # use post ICU clustering
  dplyr::left_join(z1 %>% dplyr::select(Patient_ID,Cluster, color),
                   by="Patient_ID") %>%
  dplyr::mutate(Cluster=as.factor(Cluster)) %>%
  dplyr::left_join(anno %>% dplyr::select(Patient_ID,FollowUp_12m) %>% distinct(), by="Patient_ID")

res_ang2 <- lapply("LogANG2" %>% {names(.)=.;.}, function(x){
  dat <- ang2_diff %>%
    dplyr::mutate(var=!!sym(x)) %>%
    # dplyr::filter(Cluster!=3) %>%
    {.}
  fit <- lm(formula = var ~ Cluster , data = dat)
  list(statistic=summary(fit)$coefficients[2,"Estimate"],p.value=summary(fit)$coefficients[2,"Pr(>|t|)"])
}) %>% {do.call(rbind,.)} %>% as.data.frame

# create boxplots
p_a <- dt %>%
  ggplot(aes(x=Timepoint,y=LogANG2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  geom_line(aes(group = Patient_ID, color=Patient_ID),size=1) +
  scale_color_manual(values = setNames(z1$color,z1$Patient_ID)) +
  theme_bw() +
  ggtitle(sprintf("LogANG2 p-value:%.2e",res_ang2$p.value))

# save boxplots to file
pdf(sprintf("%s/SupplementaryFigure_6E-F.pdf", Sys.Date()), width = 10, height = 6, onefile = F)
ggarrange(p_pl,p_a,
          ncol=2,nrow=1,
          common.legend = T)
dev.off()
