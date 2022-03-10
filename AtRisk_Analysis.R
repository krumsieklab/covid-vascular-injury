#### Initialize ----

rm(list = ls())
# set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(magrittr)
library(tibble)
library(maplet)
library(readxl)
library(glue)
library(openxlsx)

library(RColorBrewer)
library(pheatmap)

library(ggplot2)
library(ggpubr)
library(ggrepel)

source("HelperFunctions.R")

# check if directory to save results exists, otherwise create
if(!dir.exists(paste0(Sys.Date(),sep=""))) {
  dir.create(paste0(Sys.Date(),sep=""))
}

#### Download Data ---

# download preprocessed proteomics data from figshare
# data files will be saved in data/
file_data <- "AtRisk.xlsx"
load.web.file(
  url="https://figshare.com/ndownloader/files/34347179",
  md5sum = "c705003b12c093507673c06e64165e12",
  outfile = file_data
)

file_data_other <- "OtherProteomics.xlsx"
load.web.file(
  url="https://figshare.com/ndownloader/files/34347242",
  md5sum = "ddbd0b10b1df1dfd3d9b062a6df16810",
  outfile = file_data_other
)

#### Load preprocessed data ----

Dyy <- 
  # load proteomics data
  maplet::mt_load_xls(file=sprintf("data/%s",file_data),sheet="data",samples_in_rows=T) %>%
  # load sample annotations
  maplet::mt_anno_xls(file=sprintf("data/%s",file_data),sheet="sampleanno",anno_type="samples",anno_id_col="sample_id_match",data_id_col="sample") %>%
  # load sample annotations
  maplet::mt_anno_xls(file=sprintf("data/%s",file_data),sheet="proteinanno",anno_type="features",anno_id_col="feature_id",data_id_col="name") %>%
  # flag as logged
  maplet::mt_load_flag_logged() %>%
  {.}

#### Run analysis on joint dataset ----

# define outcomes
outcomes <- data.frame(outcome = c("death","platelet"),
                       outcomeType = c("binary","numeric"))

# set significance threshold
p.adj.cut <- 0.1

# filter samples
Dzz <- Dyy %>% 
  # filter anchor samples from second batch
  maplet::mt_modify_filter_samples(filter = !(SET=="new" & anchor_autopsy=="Yes")) %>%
  # scale
  maplet::mt_pre_trans_scale() %>%
  # select only COVID samples
  maplet::mt_modify_filter_samples(filter = Status=="COVID") %>%
  {.}
  
# loop over outcomes
for (i in 1:nrow(outcomes)) {
  Dzz %<>% 
    # perform differential analysis
    p_analyze(outcome_info = outcomes[i,],
              compname = sprintf("%s", outcomes$outcome[[i]]))
}

#### Define Protein Set ----

res <- lapply(outcomes$outcome %>% {names(.)=.;.}, function(x){
  Dzz %>% 
    # extract statistical result table
    maplet::mtm_get_stat_by_name(name=x) %>%
    # merge protein names
    dplyr::left_join(Dzz %>% rowData %>% as.data.frame %>% dplyr::select(OlinkID, Assay), by=c("var"="OlinkID")) %>%
    # rename column
    dplyr::rename(name=Assay) %>%
    # sort by adjusted p-value
    dplyr::arrange(p.adj) %>%
    # filter significant results
    dplyr::filter(p.adj < p.adj.cut) %>%
    {.}
})

# select proteins that are significantly associate to both Mortality and Platelet 
sel1 <- intersect(res$death$name, res$platelet$name)
# add some known vascular injury proteins
sel2 <- c("ADAM-TS13","CD40-L","EGFR","TIE2","SELP","uPA","VEGFA","GP6","HO-1")
# merge the two lists
prot <- c(sel1, sel2) %>% unique

# save protein set to file
save(prot, file=sprintf("%s/Protlist.Rdata", Sys.Date()))

#### Print statistical results to file ----

res_full <- lapply(outcomes$outcome %>% {names(.)=.;.}, function(x){
  Dzz %>%
    maplet:::mtm_get_stat_by_name(name=x) %>%
    # merge protein names
    dplyr::left_join(Dzz %>% rowData %>% as.data.frame %>% dplyr::select(OlinkID,Assay), by=c("var"="OlinkID")) %>%
    dplyr::rename(name=Assay) %>%
    dplyr::arrange(p.adj) %>%
    {.}
})

wb = createWorkbook()
sheet = addWorksheet(wb, "Mortality")
writeData(wb, sheet=sheet, res_full$death %>% dplyr::select(name,estimate,std.error,statistic,df,fc,p.value,p.adj),
          rowNames = F, colNames = T)
sheet = addWorksheet(wb, "Platelet")
writeData(wb, sheet=sheet, res_full$platelet %>% dplyr::select(name,estimate,std.error,statistic,df,p.value,p.adj),
          rowNames = F, colNames = T)
saveWorkbook(wb, sprintf("%s/AtRisk_DifferentialAnalysisResults.xlsx",Sys.Date()), overwrite = TRUE)

#### Figure 2A ----

p_prset <- res_full$death %>%
  dplyr::mutate(death_association = -log10(p.adj)) %>%
  dplyr::mutate(sign=ifelse(p.adj<p.adj.cut,TRUE,FALSE)) %>%
  dplyr::select(var, name, death_association, sign) %>%
  dplyr::left_join(res_full$platelet %>%
                     dplyr::mutate(platelet_association = -log10(p.adj)) %>%
                     dplyr::mutate(sign=ifelse(p.adj<p.adj.cut,TRUE,FALSE)) %>%
                     dplyr::select(var, name, platelet_association, sign), 
                   by=c("var","name"), suffix=c("_death","_platelet")) %>%
  dplyr::mutate(sign=ifelse(sign_death==sign_platelet, sign_death, ifelse(name %in% sel2, TRUE, FALSE))) %>%
  dplyr::mutate(label=ifelse(sign==TRUE,name,ifelse(name=="TIE2",name,NA))) %>%
  dplyr::mutate(color=case_when(
    (sign_death==TRUE & sign_platelet==TRUE) ~ "Both",
    (sign_death==TRUE & sign_platelet==FALSE) ~ "Death only",
    (sign_death==FALSE & sign_platelet==TRUE) ~ "Platelets only",
    (sign_death==FALSE & sign_platelet==FALSE) ~ ifelse(name=="TIE2",name,"Neither")
  )) %>% 
  ggplot(aes(x=death_association, y=platelet_association, label=label, color=color)) +
  geom_point() +
  geom_text_repel(force=1.5) +
  # geom_vline(xintercept = log10(p.adj.cut), linetype="dotted") +
  geom_vline(xintercept = -log10(p.adj.cut), linetype="dotted") +
  # geom_hline(yintercept = log10(p.adj.cut), linetype="dotted") +
  geom_hline(yintercept = -log10(p.adj.cut), linetype="dotted") +
  xlim(c(0,4)) +
  ylim(c(0,4)) +  
  theme_bw() +
  xlab("Association with Death") +
  ylab("Association with Platelets") +
  scale_color_manual(values=c(Both="black", 
                              `Death only`="#F05349", 
                              `Platelets only`="#0AAD77", 
                              TIE2="gray40",
                              Neither="gray80"),
                     name="Significant\nassociations") +
  ggtitle("Statistical Results")

pdf(sprintf("%s/Figure_2A.pdf", Sys.Date()), width = 7, height = 5)
print(p_prset)
dev.off()

#### Figure 2B ----

# create dataframe for plotting
dt <- Dzz %>% assay %>% as.data.frame %>% t
colnames(dt) <- rowData(Dzz)$Assay
dt <- cbind.data.frame(dt,
                       Mortality=Dzz$death,
                       Status=Dzz$Status,
                       Patient=Dzz$Patient_ID,
                       Platelet = Dzz$platelet) %>%
  dplyr::filter(Status=="COVID") %>%
  dplyr::select(-Status)
colnames(dt) <- make.names(colnames(dt))

# get statistical results for Mortality and Platelet
res <- lapply(outcomes$outcome %>% {names(.)=.;.}, function(name) {
  r <- maplet:::mtm_get_stat_by_name(Dzz, name=name) %>%
    dplyr::left_join(Dzz %>% rowData %>% as.data.frame %>% dplyr::select(OlinkID,Assay), by=c("var"="OlinkID")) %>%
    dplyr::rename(labels=var) %>%
    dplyr::rename(var=Assay) %>%
    {.}
  pvals <- r[,which(colnames(r) %in% c("var","p.value", "statistic", "p.adj"))]
})

# mortality boxplot
g <- lapply(prot %>%  {names(.)=.;.}, function(p){
  pp <- make.names(p)
  dt %>%
    dplyr::select(Mortality, !!sym(pp), Patient) %>%
    dplyr::rename(Prot=!!sym(pp)) %>%
    dplyr::group_by(Patient) %>%
    # compute protein average across repeated measurements for plotting
    dplyr::summarise(Mortality=unique(Mortality), Patient=unique(Patient), ProtAverage=mean(Prot)) %>%
    dplyr::mutate(p.adj=res$death$p.adj[res$death$var==p]) %>%
    dplyr::arrange(p.adj) %>%
    ggplot(aes(x=Mortality, y=ProtAverage, fill=Mortality)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position=position_jitterdodge()) +
    scale_fill_manual(values = c("Yes"="#F05349","No"="#4577EA")) +
    theme_bw() +
    ylab(p) +
    ggtitle(sprintf("%s: p.adj %.2e", p, res$death$p.adj[res$death$var==p]))
})

# save to file
ncol=5
pdf(sprintf("%s/Figure_2B.pdf", Sys.Date()), width = 5*ncol, height = 5*ceiling(length(prot)/5), onefile = F)
print(ggarrange(plotlist = g,
          common.legend = T,
          ncol=ncol, nrow=ceiling(length(prot)/ncol)))
dev.off()

#### Figure 2C ----

D1 <- Dzz %>% 
  # filter only protein set
  maplet::mt_modify_filter_features(filter= Assay %in% prot) %>%
  {.}

# define heatmap breaks
breaksList = seq(-max(abs(min(assay(D1))), abs(max(assay(D1)))), 
                 max(abs(min(assay(D1))), abs(max(assay(D1)))), by=0.1)

# get proteomics data frame
dthh <- D1 %>% assay  %>% t %>% as.data.frame %>%
  dplyr::mutate(Patient = D1$Patient_ID) %>%
  dplyr::group_by(Patient) %>%
  # compute protein average across repeated measurements for plotting
  dplyr::summarise_all(mean) %>%
  dplyr::mutate(Patient=sprintf("Patient_%s",Patient)) %>%
  tibble::column_to_rownames("Patient") %>%
  t()

# get clinical annotations
anno <- D1 %>% colData %>% as.data.frame %>% 
  dplyr::select(Patient_ID,death,platelet,SET,age) %>%
  distinct() %>%
  # create autopsy variable
  dplyr::mutate(autopsy=ifelse(SET=="new","Yes","No")) %>%
  # patient 835 had no tissue staining
  dplyr::mutate(autopsy=ifelse(Patient_ID=="835","No",autopsy)) %>%
  dplyr::select(-SET)
rownames(anno) <- sprintf("Patient_%s",anno$Patient_ID)
# define average protein abundance
anno$averageProt <- colSums(dthh[,match(rownames(anno),colnames(dthh))])/length(prot)

# define annotation colors
ann_colors = list(
  death = c(Yes="#F05349", No="#4577EA"),
  platelet = colorRampPalette(c("#F0EC19","#0AAD77") )(nrow(anno)),
  age = colorRampPalette(c("#a5b6d2","#0f2c64") )(nrow(anno)),
  averageProt = colorRampPalette(c("#f2f2f2","#404040"))(nrow(anno)),
  autopsy = c(No="white",Yes="coral")
)

# Modify ordering of the clusters using clustering callback option
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv, agglo.FUN = function(x){min(x)+max(x)})
  as.hclust(dend)
}

# heatmap
n_clust <- 3
out2 <- pheatmap(mat = dthh, 
                 # cellwidth = 10, cellheight=10,
                 annotation_col = anno %>% dplyr::select(-Patient_ID, -autopsy),
                 annotation_colors = ann_colors,
                 clustering_method = "ward.D2", 
                 cutree_cols = n_clust, cutree_rows = 2,
                 # show_colnames = F,
                 labels_row = D1 %>% rowData %>% .$Assay,
                 clustering_callback = callback,
                 border_color = NA,
                 color = colorRampPalette(c("#4575B4","white","#D73027"))(length(breaksList)),
                 breaks = breaksList)

# save heatmap to file
pdf(sprintf("%s/Figure_2C.pdf", Sys.Date()), width = 12, height = 7)
print(out2)
dev.off()

#### Supplementary Figure 3 ---- 

# extract clusters and get annotations for Fisher's test
z <- data.frame(Cluster=sort(cutree(out2$tree_col, k=n_clust))) %>% 
  tibble::rownames_to_column(var="subject_id") %>% 
  dplyr::left_join(anno %>% dplyr::mutate(subject_id=sprintf("Patient_%s",Patient_ID))) %>%
  dplyr::select(Cluster, subject_id, averageProt, age) %>%
  # dplyr::mutate(Cluster=ifelse(Cluster==2,"High Vascular Injury", ifelse(Cluster==1,"Intermediate","Low Vascular Injury"))) %>%
  {.}

sortvar <- z %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarize(var=mean(averageProt)) %>%
  dplyr::arrange(var) %>%
  dplyr::mutate(Cluster_ord=rank(-var, ties.method = "min"))

z %<>% 
  dplyr::left_join(sortvar,by="Cluster") %>%
  dplyr::arrange(Cluster_ord)

# compute differential analysis between clusters for averageProt
da <- data.frame(avProt12 = wilcox.test(x=z$averageProt[z$Cluster_ord==1], y=z$averageProt[z$Cluster_ord==2])$p.value,
                 avProt23 = wilcox.test(x=z$averageProt[z$Cluster_ord==2], y=z$averageProt[z$Cluster_ord==3])$p.value,
                 avProt13 = wilcox.test(x=z$averageProt[z$Cluster_ord==1], y=z$averageProt[z$Cluster_ord==3])$p.value,
                 age12 = wilcox.test(x=z$age[z$Cluster_ord==1], y=z$age[z$Cluster_ord==2])$p.value,
                 age23 = wilcox.test(x=z$age[z$Cluster_ord==2], y=z$age[z$Cluster_ord==3])$p.value,
                 age13 = wilcox.test(x=z$age[z$Cluster_ord==1], y=z$age[z$Cluster_ord==3])$p.value)

# create boxplot
p_avProt <- z %>%
  ggplot(aes(x=as.factor(Cluster_ord), y=averageProt)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(width = 0.2)) +
  theme_bw() +
  xlab("Cluster") +
  ylab("averageProt") +
  ggtitle(sprintf("averageProt: 1-2: %.2e; 2-3: %.2e; 1-3: %.2e", da$avProt12,da$avProt23,da$avProt13))

# create boxplot
p_age <- z %>%
  ggplot(aes(x=as.factor(Cluster_ord), y=age)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(width = 0.2)) +
  theme_bw() +
  xlab("Cluster") +
  ylab("Age") +
  ggtitle(sprintf("age: 1-2: %.2e; 2-3: %.2e; 1-3: %.2e ", da$age12, da$age23, da$age13))

# save plot to file
pdf(sprintf("%s/SupplementaryFigure_3.pdf", Sys.Date()), width=12, height=6, onefile = F)
print(ggarrange(p_avProt,p_age,
          ncol=2, nrow=1,
          labels="AUTO",
          common.legend = T))
dev.off()

#### Figure 4B ----

# load autopsy quantification data
df <- read_excel(path="data/OtherProteomics.xlsx",sheet = "Autopsy_quantification",col_names = T) %>%
  dplyr::mutate(Group=factor(Group,levels = c("Low ANGPT2","High ANGPT2")))

# boxplot
p <- df %>%
  ggplot(aes(x=Group,y=log10(CD61))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitter(width = 0.2)) +
  theme_bw() +
  ggtitle(sprintf("CD61 p-value %.2e",wilcox.test(x=log10(df$CD61[df$Group=="High ANGPT2"]),y=log10(df$CD61[df$Group=="Low ANGPT2"]))$p.value))

# save to file
pdf(sprintf("%s/Figure_4B.pdf",Sys.Date()), width = 5, height = 5)
print(p)
dev.off()

#### Figure 4C ----

# define autopsy patients identifiers
autopsy_patients <- data.frame(Patient_ID=c(832,833,6),
                               subject_id=c("Patient_832","Patient_833","Patient_6"),
                               label=c("P1","P2","P3"))

# create dataframe for plotting 
df <- dthh %>% 
  as.data.frame %>%
  dplyr::select(any_of(autopsy_patients$subject_id)) %>%
  t %>% as.data.frame %>%
  tibble::rownames_to_column("subject_id") %>%
  reshape2::melt() %>%
  dplyr::left_join(autopsy_patients, by="subject_id") %>%
  dplyr::mutate(dummyvar=case_when(label=="P3"~3,
                                   label=="P2"~2,
                                   label=="P1"~1)) %>%
  dplyr::left_join(D1 %>% rowData %>% as.data.frame %>% dplyr::select(name, Assay), by=c("variable"="name")) %>%
  dplyr::rename(protein_name=Assay) %>%
  dplyr::mutate(text=ifelse(label=="P3",protein_name,""))

# create vector of deltas for scatterplot
b <- runif(nrow(df), -0.1, 0.1)

# plot boxplot
gg <- df %>%
  ggplot(aes(label=text)) +
  geom_boxplot(aes(x=dummyvar, y=value, group=subject_id), 
               outlier.shape = NA) +
  geom_point(aes(x=dummyvar + b, y=value)) +
  # geom_hline(yintercept = df$value %>% median, size=1) +
  geom_text_repel(aes(x=dummyvar + b, y=value)) +
  geom_line(aes(x=dummyvar + b, y=value, group=variable), 
            color="gray60", size=0.3) +
  scale_x_continuous(breaks = c(1,2,3), labels = c("P1", "P2","P3"))+
  xlab("Patient") +
  ylab("") +
  theme_bw() +
  ggtitle(sprintf("p-value 1-2: %.2e / p-value 2-3: %.2e / p-value 1-3: %.2e",
                  wilcox.test(x=df$value[df$dummyvar %in% c(1)],y=df$value[df$dummyvar %in% c(2)], paired = T) %>% .$p.value,
                  wilcox.test(x=df$value[df$dummyvar %in% c(2)],y=df$value[df$dummyvar %in% c(3)], paired = T) %>% .$p.value,
                  wilcox.test(x=df$value[df$dummyvar %in% c(1)],y=df$value[df$dummyvar %in% c(3)], paired = T) %>% .$p.value))

# save plot to file
pdf(sprintf("%s/Figure_4C.pdf", Sys.Date()), width = 8, height = 5)
print(gg)
dev.off()
