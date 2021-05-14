# Downloads files from the web and verifies their checksum. Will use local copy in current directory, if it exists
load.web.file <- function(
  url, md5sum, outfile, zipfile = F) {
  # check if local file exists
  if (file.exists(outfile)) {
    # verify checksum
    realsum <- tools::md5sum(outfile)[[1]]
    if (realsum != md5sum) stop(sprintf("Local file %s has wrong checksum: %s", outfile, realsum))
    # do not delete wrong file, it was already here before
    
  } else {
    if(zipfile){
      # download file
      temp <- tempfile()
      download.file(url,temp)
      unzip(zipfile = temp, files = outfile, exdir = ".")
    } else {
      # download file
      download.file(url, outfile)
    }
    # verify checksum
    realsum <- tools::md5sum(outfile)[[1]]
    if (realsum != md5sum) { 
      # delete wrong file
      unlink(outfile)
      stop(sprintf("Remote file %s has wrong checksum: %s", url, realsum))
    }
  }
}

# Modify ordering of the clusters using clustering callback option
callback2 = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv, agglo.FUN = function(x){min(x)} )
  as.hclust(dend)
}

# differential analysis function
p_analyze <- function(D, outcome_info, compname) {
  D %<>%
    # linear mixed-effect model
    maplet::mt_stats_univ_lm(formula = as.formula(glue("~{outcome_info$outcome} + Time + SET + (1|Patient_ID)")), 
                             stat_name = compname) %>%
    # multiple testing correction
    maplet::mt_post_multtest(stat_name = compname, method = "BH") %>%
    {.}
  
  if(outcome_info$outcomeType=="binary") {
    D %<>%
      # fold change
      maplet::mt_post_fold_change(stat_name = compname) %>%
      # volcano plot
      maplet::mt_plots_volcano(stat_name = compname, 
                               x = fc, 
                               feat_filter = p.adj < !!p.adj.cut,
                               colour = p.adj < !!p.adj.cut)
  }
  if(outcome_info$outcomeType=="numeric") {
    D %<>%
      # volcano plot
      maplet::mt_plots_volcano(stat_name = compname, 
                               x = statistic, 
                               feat_filter = p.adj < !!p.adj.cut,
                               colour = p.adj < !!p.adj.cut)
  }
  D
}

# load protein set
load_protein_list <- function(filename="data/Protlist.Rdata"){
  # protein list defined by GIM analysis
  if(file.exists(filename)) {
    warning(sprintf("Loading protein set saved in %s",filename))
    load(filename)
  } else{
    warning("Cannot find specified file.\nLoading precomputed protein set in data/Protlist.Rdata")
    load("data/Protlist.Rdata")
  }
  prot
}




