#Soren Charmsaz
#Last updated 03/26/2024
#R version 4.1.2 (2021-11-01)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: Mac OS Sonoma 14.5

rm(list = ls())
library(reshape2)
library(randomcoloR)
library(pals)
library(ggplot2)
library(Hmisc)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(corrplot)
library(readxl)
library(ggridges)
library(ggpubr)
library(raster)
library(matrixStats)
library(limma)

####READ and CLUSTER FUNCTIONS####
returnfcs <- function(FDR_cutoff=.05,
                      metaDataFile='~/',
                      panelDataFile='~/',
                      dataDirectory='~/',
                      shape_timepoint=NULL,
                      color_timepoint=NULL){
  #This function generates an fcs file, subtype_markers, colors and shapes for clustering 
  require(scales);require(readxl);require(dplyr);require(flowCore)
  ##directory and metadatafile checking
  if(!dir.exists(dataDirectory)) {stop('ERR: cannot find data directory')}
  if(!file.exists(metaDataFile)) {stop('ERR: cannot find metadata.xlsx or .csv file')}
  ##readin metadata and clean
  ifelse(grepl(metaDataFile,pattern='.xls'),md <- read_excel(metaDataFile),md <- read.csv(metaDataFile,header = TRUE))#must be in xl format or csv
  md$batch <- factor(md$Batch)
  md$timepoint <- factor(md$timepoint)
  md$run <- factor(md$Run)
  
  rownames(md) = md$sample_id;md$sample_id <- md$sample_id
  #Make sure all files in metadata present in datadirectory
  if(!all(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])){
    print(paste('ERR: not all filenames in metadata present in data folder - missing',md$file_name[!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])],'Subsetting...'))
    md <- md[-c(!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])),]
  }
  ##Define shapes for timepoint
  if(is.null(shape_timepoint)){shape_timepoint <- c(0:25)[1:length(levels(md$timepoint))]}#can specify as long as number is same
  if(length(shape_timepoint)!=length(levels(md$timepoint))){stop(paste0('ERR no. shapes specified is less than no. of timepoint (',length(levels(md$timepoint)),')'))}
  names(shape_timepoint) <- levels(md$timepoint)
  ## Define colors for the timepoint
  if(is.null(color_timepoint)){color_timepoint <- hue_pal()(length(levels(md$timepoint)))}#can specify as long as number is same
  if(length(color_timepoint)!=length(levels(md$timepoint))){stop(paste0('ERR no. shapes specified is less than no. of timepoint (',length(levels(md$timepoint)),')'))}
  ## read fcs
  fcs_raw <- read.flowSet(md$file_name, path = dataDirectory, transformation = FALSE, truncate_max_range = FALSE)
  #sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))
  panel <- read_excel(panelDataFile)
  head(data.frame(panel))
  ## Replace problematic characters
  panel$Metal <- gsub('-', '_', panel$Metal)
  panel_fcs <- pData(parameters(fcs_raw[[1]]))
  panel_fcs$desc <- gsub('-', '_', panel_fcs$desc)
  panel_fcs$desc[is.na(panel_fcs$desc)] <- paste0('NA_',which(is.na(panel_fcs$desc))) #was labelled 'todo'(mapping based on isotope for now just getting rid of NA keeping rownum) 
  # use panel$Antigen to fix description in panel_fcs
  # use metal+isotope as mapping between panel from xlsx and panel from the fcs files
  rownames(panel_fcs) = panel_fcs$name
  panel_fcs[paste0(panel$Metal,panel$Isotope,'Di'),2] <- panel$Antigen
  ## Replace paramater data in flowSet
  pData(parameters(fcs_raw[[1]])) <- panel_fcs
  ## Define variables indicating marker types
  subtype_markers <- panel$Antigen[panel$Subtype == 1]
  functional_markers <- panel$Antigen[panel$Functional == 1]
  if(!all(subtype_markers %in% panel_fcs$desc)){stop('ERR: Not all subtype_markers in panel_fcs$desc (isotopes)')}
  if(!all(functional_markers %in% panel_fcs$desc)){stop('ERR: Not all functional_markers in panel_fcs$desc (isotopes)')}
  ## arcsinh transformation and column subsetting
  fcs <- fsApply(fcs_raw, function(x, cofactor = 5){
    colnames(x) <- panel_fcs$desc
    expr <- exprs(x)
    expr <- asinh(expr[, union(subtype_markers,functional_markers)] / cofactor)
    expr[!is.finite(expr)] <- NA #convert inf to NA
    expr<-na.omit(expr) #remove NA
    exprs(x) <- expr
    x
  })
  sample_ids <- rep(md$sample_id, fsApply(fcs, nrow))
  return(list('fcs'=fcs,
              'subtype_markers'=subtype_markers,
              'functional_markers'=functional_markers,
              'shape_timepoint'=shape_timepoint,
              'color_timepoint'=color_timepoint,
              'sample_ids'=sample_ids,
              'meta_data'=md))
}

clusterfcs <- function(fcs=output$fcs,
                       subtype_markers = output$subtype_markers,
                       seed=1234,plottitle='consensus_plots',
                       numclusters=40){
  ## Cell population identification with FlowSOM and ConsensusClusterPlus
  require(dplyr);require(FlowSOM);require(ConsensusClusterPlus)
  set.seed(seed)
  som <- ReadInput(fcs, transform = FALSE, scale = FALSE) %>% BuildSOM(colsToUse = subtype_markers)
  ## Get the cell clustering into 100 SOM codes
  cell_clustering_som <- som$map$mapping[,1]
  ## Metaclustering into numclusters with ConsensusClusterPlus
  codes <- som$map$codes
  mc <- ConsensusClusterPlus(t(codes), maxK = numclusters, reps = 100,
                             pItem = 0.9, pFeature = 1, title = plottitle, 
                             plot = "png", clusterAlg = "hc", 
                             innerLinkage = "average", finalLinkage = "average",
                             distance = "euclidean", seed = 1234)
  
  ## Get cluster ids for each cell
  code_clustering <- mc[[numclusters]]$consensusClass#metaclusters consensus
  cell_clustering <- code_clustering[cell_clustering_som]#cell clustering from som
  return(list('code_clustering'=code_clustering,'cell_clustering'=cell_clustering,'metaclusters'=mc))
}

####CLUSTER HEATMAP FUNCTIONS ####
plot_clustering_heatmap_wrapper <- function(fcs, cell_clustering, nclusters=40,
                                            color_clusters='auto', cluster_merging = NULL, 
                                            subtype_markers,
                                            clusterMergeFile=NULL,
                                            fileName = 'clusteringheatmap.pdf'){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales)
  ## Will output the heatmap object and print it 
  if((color_clusters)=='auto'){color_clusters <- hue_pal()(nclusters)}
  #get expression
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  ## Scale expression of all markers to values between 0 and 1
  rng <- colQuantiles(expr, probs = c(0.01, 0.99))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,subtype_markers]
  ## Calculate the mean expression
  pdf(fileName, width=8, height=11) 
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  ## Calculate cluster frequencies
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  ## Sort the cell clusters with hierarchical clustering
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Colors for the heatmap
  color_heat <- colorRampPalette(brewer.pal(n = 9, name = "YlOrBr"))(100)
  #legend_breaks = seq(from = 0, to = 1, by = 0.2)
  labels_row <- paste0(expr01_mean$cell_clustering, " (", clustering_prop ,
                       "%)")
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  ## Annotation for the merged clusters
  
  if(!is.null(clusterMergeFile)){
    ifelse(grepl(clusterMergeFile,pattern='.xls'),cluster_merging <- read_excel(clusterMergeFile),cluster_merging <- read.csv(clusterMergeFile,header = TRUE))
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$Merged <- cluster_merging$new_cluster
    color_clusters2 <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters2) <- levels(cluster_merging$new_cluster)
    annotation_colors$Merged <- color_clusters2
  }
  
  p <- pheatmap(expr_heat, color = rev(colorRampPalette(brewer.pal(n = 11, name = "RdYlBu"))(100)), 
                cluster_cols = T,
                cluster_rows = T, 
                labels_row = labels_row,
                #scale="column",
                display_numbers = FALSE, number_color = "black",
                fontsize = 9, fontsize_number = 6,  
                #legend_breaks = legend_breaks,
                annotation_row = annotation_row, 
                annotation_colors = annotation_colors,
                cellwidth = 8,
                cellheight = 8
  )
  dev.off() 
  print('Colors:')
  print(color_clusters)
  print(p);return(p)
}

plot_clustering_heatmap_wrapper2 <- function(fcs, cell_clustering, nclusters=40,
                                             color_clusters='auto',
                                             subtype_markers,
                                             fileName = 'clusteringheatmap.pdf'){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales)
  ## Will output the heatmap object and print it 
  if((color_clusters)=='auto'){color_clusters <- hue_pal()(nclusters)}
  #get expression
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  ## Scale expression of all markers to values between 0 and 1
  rng <- colQuantiles(expr, probs = c(0.01, 0.99))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,subtype_markers]
  ## Calculate the mean expression##################################################
  pdf(fileName, width=8, height=11) 
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Colors for the heatmap
  
  labels_row <- paste0(expr01_mean$cell_clustering, " ")
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  p <- pheatmap(expr_heat, 
                color = c(rep(magma(100)[1],25),magma(100)[1:100]), 
                cluster_cols = T,
                cluster_rows = F, 
                labels_row = labels_row,
                #scale="column",
                display_numbers = F, 
                number_color = "black",
                fontsize = 9, fontsize_number = 6,  
                #legend_breaks = legend_breaks,
                annotation_row = annotation_row, 
                annotation_colors = annotation_colors,
                cellwidth = 8,
                cellheight = 8,
                border_color = "black",
                annotation_legend = F
  )
  dev.off() 
  print('Colors:')
  print(color_clusters)
  print(p);return(p)
}

####DIAGNOSTICS####
makeDiagnosticPlots = function(exprData, 
                               md = output$meta_data,
                               sample_ids = output$sample_ids,
                               fcs = output$fcs,
                               subtype_markers = output$subtype_markers,
                               color_conditions = clustercolors,
                               shape_conditions = c(1:13),
                               fileName = 'diagnostics.pdf', 
                               tit = '', 
                               fun = mean)
{
  pdf(file = fileName)
  
  # plot 1
  ggdf <- data.frame(sample_id = sample_ids, exprData)
  ggdf <- melt(ggdf, id.var = 'sample_id', value.name = 'expression', 
               variable.name = 'antigen')
  mm <- match(ggdf$sample_id, md$sample_id)
  ggdf$timepoint <- md$timepoint[mm]
  print(ggplot(ggdf, aes(x = expression, color = timepoint, group = sample_id)) + 
          geom_density() +
          facet_wrap(~ antigen, nrow = 4, scales = 'free') + theme_bw() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1), 
                strip.text = element_text(size = 7),
                axis.text = element_text(size = 5)) + 
          scale_color_manual(values = color_timepoint) )
  # plot 2
  
  ## Spot check - number of cells per sample
  cell_table <- table(sample_ids)
  ggdf <- data.frame(sample_id = names(cell_table), 
                     cell_counts = as.numeric(cell_table))
  mm <- match(ggdf$sample_id, md$sample_id)
  ggdf$timepoint <- md$timepoint[mm]
  print(ggplot(ggdf, aes(x = sample_id, y = cell_counts, fill = timepoint)) + 
          geom_bar(stat = 'identity') + 
          geom_text(aes(label = cell_counts), hjust = 0.5, vjust = -0.5, size = 2.5) + 
          theme_bw() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +  
          scale_fill_manual(values = color_timepoint, drop = FALSE) + 
          scale_x_discrete(drop = FALSE))
  
  # plot 3
  
  ## Define a function that calculates the Non-Redundancy Score per sample
  NRS <- function(x, ncomp = 3){
    pr <- prcomp(x, center = TRUE, scale. = FALSE)
    score <- rowSums(outer(rep(1, ncol(x)), pr$sdev[1:ncomp]^2) * 
                       abs(pr$rotation[,1:ncomp]))
    return(score)
  }
  
  ## Calculate the score
  ## May want to do the same with other markers
  nrs_sample <- fsApply(fcs[, subtype_markers], NRS, use.exprs = TRUE)
  rownames(nrs_sample) <- md$sample_id
  nrs <- colMeans(nrs_sample, na.rm = TRUE)
  
  ## Plot the NRS for ordered markers
  ## May be helpful to look at tissue instead of condition
  subtype_markers_ord <- names(sort(nrs, decreasing = TRUE))
  nrs_sample <- data.frame(nrs_sample)
  nrs_sample$sample_id <- rownames(nrs_sample)
  ggdf <- melt(nrs_sample, id.var = "sample_id",
               value.name = "nrs", variable.name = "antigen")
  ggdf$antigen <- factor(ggdf$antigen, levels = subtype_markers_ord)
  mm <- match(ggdf$sample_id, md$sample_id)
  ggdf$timepoint <- md$timepoint[mm]
  print(ggplot(ggdf, aes(x = antigen, y = nrs)) +
          geom_point(aes(color = timepoint), alpha = 0.9,
                     position = position_jitter(width = 0.3, height = 0)) +
          geom_boxplot(outlier.color = NA, fill = NA) +
          stat_summary(fun.y = "mean", geom = "point", shape = 21, fill = "white") +
          theme_bw() + ggtitle(tit)+ 
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) )
  # scale_color_manual(values = color_conditions)
  
  dev.off()
}


####CLUSTER HISTO####

plot_clustering_distr_wrapper <- function(expr = expr, 
                                          cell_clustering){
  # Calculate the median expression
  cell_clustering <- factor(cell_clustering)
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  # Sort the cell clusters with hierarchical clustering
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  # Calculate cluster frequencies
  freq_clust <- table(cell_clustering)
  freq_clust <- round(as.numeric(freq_clust)/sum(freq_clust)*100, 2)
  cell_clustering <- factor(cell_clustering,
                            labels = levels(cell_clustering))
  ### Data organized per cluster
  ggd <- melt(data.frame(cluster = cell_clustering, expr),
              id.vars = "cluster", value.name = "expression",
              variable.name = "antigen")
  ggd$antigen <- factor(ggd$antigen, levels = colnames(expr))
  ggd$reference <- "no"
  ### The reference data
  ggd_bg <- ggd
  ggd_bg$cluster <- "reference"
  ggd_bg$reference <- "yes"
  
  ggd_plot <- rbind(ggd, ggd_bg)
  ggd_plot$cluster <- factor(ggd_plot$cluster,
                             levels = c(levels(cell_clustering)[rev(cluster_rows$order)], "reference"))
  
  ggplot() +
    geom_density_ridges(data = ggd_plot, aes(x = expression, y = cluster,
                                             color = reference, fill = reference), alpha = 0.3) +
    facet_wrap( ~ antigen, scales = "free_x", nrow = 2) +
    theme_ridges() +
    theme(axis.text = element_text(size = 7),  
          strip.text = element_text(size = 7), legend.position = "none")
  
}
####UMAP####
do_umap <- function(fcs,subtype_markers,sample_ids,cell_clustering,metadata,
                    clusterMergeFile='~_merging.xlsx',
                    seed = 1234, ncells=2000,sample_subset=NULL){
  require(umap);require(flowCore);require(readxl)
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  ## Create vector to later find and skip duplicates
  dups <- duplicated(expr[, subtype_markers])
  dups <- which(!(dups))## Find and skip duplicates
  ifelse(grepl(clusterMergeFile,pattern='.xls'),cluster_merging <- read_excel(clusterMergeFile),cluster_merging <- read.csv(clusterMergeFile,header = TRUE))
  ## New clustering1m
  mm <- match(cell_clustering, cluster_merging$original_cluster)
  cell_clustering1m <- cluster_merging$new_cluster[mm]
  ## Create a data frame of sample_ids and cell_clustering1m
  dtf<-data.frame(ids=sample_ids,type=cell_clustering1m)
  ## Data subsampling: create indices by sample
  inds <- split(1:length(sample_ids), sample_ids) #to get original indexes belonging to each cluster
  samplenames <- names(inds) #create a name vector of the files
  custom.settings = umap.defaults
  custom.settings$seed = seed
  custom.settings$n.neighbors = neighbors
  ####umapindex generation####
  ifelse(is.null(sample_subset),
         umap_ncells <- pmin(table(sample_ids), ncells),
         umap_ncells <- pmin(table(sample_ids), ncells)[sample_subset]
  )
  if(!is.null(sample_subset)){inds <- inds[sample_subset]}
  umap_inds <- lapply(names(inds), function(i){
    s <- sample(inds[[i]], umap_ncells[i], replace = FALSE)
    intersect(s, dups)
  })
  set.seed(seed)
  umap_inds <- unlist(umap_inds)
  umap_out <- umap(expr[umap_inds, subtype_markers], config = custom.settings, method = 'naive')
  umapRes2D = data.frame(umap1 = umap_out$layout[, 1], umap2 = umap_out$layout[, 2], 
                         expr[umap_inds, subtype_markers],
                         sample_id = sample_ids[umap_inds], cell_clustering = factor(cell_clustering1m[umap_inds]), check.names = FALSE)
  return(umapRes2D)
}

plotUmap <- function(umapRes,seed=1234,neighbors=10,midpoint,color_clusters='auto',code_clustering,subtype_markers=NULL)
{require(umap);require(ggplot2);require(viridis);require(ggrepel)
  if((color_clusters)=='auto'){color_clusters <- hue_pal()(length(unique(code_clustering)))}
  custom.settings = umap.defaults
  custom.settings$seed = seed
  custom.settings$n.neighbors = neighbors
  ggp <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = cell_clustering)) +
    geom_point(size = 1) +
    theme_bw() +
    
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) +
    
    scale_color_manual(values = color_clusters, name="CLUSTERS") +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 2))
  
  print(ggp)
  #other options
  print(ggp + facet_wrap(~ timepoint, ncol = 3)+ggtitle('TIMEPOINTS'))
  print(ggp + facet_wrap(~ batch, ncol = 2)+ggtitle('BATCH'))
  
  
  ggp2 <- ggplot(umapRes,  aes(x = umap1, y = umap2)) +
    geom_point(size = 1) +
    theme_bw() +
    
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 2))
  print(ggp2 + facet_wrap(~ sample_id, ncol = 8)+ggtitle('SAMPLE'))
  
  
  
  
  ggp3 <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = sample_id)) +
    geom_point(size = 1) +
    theme_bw() +
    
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 2))
  
  print(ggp3)
  #can specify which markers to display
  if(!is.null(subtype_markers)){
    for(i in subtype_markers)
    {
      ggp <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = umapRes[,i])) +
        geom_point(size = 1) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        scale_color_gradient2(i, low="dark blue",mid="white",high="dark red", midpoint = mean(unlist(umapRes[,i])))
      print(ggp)
    }
  }
}


#======================
#     RUNNING DATA
#======================


####DATA LOADING AND CLUSTERING####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
workd<-getwd()

###Start here if loading output rds

output <- readRDS('backup_output.rds')

#set up factor levels
clusterlevels = c("TcEFF_I",
                  "TcEFF_II",
                  "TcEFF_III",
                  "TcEM",
                  "TcCM",
                  "TcN",
                  "ThCTL",
                  "Th2",
                  "Th2EM",
                  "Th2CM_I",
                  "Th2CM_II",
                  "Th17",
                  "Treg",
                  "ThN",
                  "DNT_I",
                  "DNT_II",
                  "DNT_III",
                  "DNT_IV",
                  "DNT_V",
                  "NK_I",
                  "NK_II",
                  "NK_III",
                  "B_I",
                  "B_II",
                  "Myeloid",
                  "UA")

samplevels <- c("468_0",
                "468_1",
                "963_0",
                "963_1",
                "615_0",
                "615_1",
                "645_0",
                "645_1",
                "732_0",
                "732_1",
                "438_0",
                "438_1",
                "133_0",
                "133_1",
                "682_0",
                "682_1",
                "134_0",
                "134_1",
                "163_0",
                "163_1",
                "582_0",
                "582_1",
                "937_0",
                "937_1",
                "523_0",
                "523_1",
                "523_2",
                "360_0",
                "360_1",
                "239_0",
                "239_1",
                "404_0",
                "404_1",
                "404_2",
                "941_0",
                "941_1",
                "789_0",
                "789_1",
                "600_0",
                "600_1",
                "290_0",
                "290_1",
                "504_0",
                "504_1",
                "064_0",
                "064_1",
                "485_0",
                "485_1",
                "718_0",
                "718_1",
                "718_2",
                "258_0",
                "258_1",
                "890_0_NA",
                "890_0",
                "890_1",
                "065_0",
                "065_1",
                "065_2",
                "617_0",
                "617_1",
                "617_2",
                "680_0",
                "680_1",
                "630_0",
                "630_1",
                "619_0",
                "619_1",
                "154_0",
                "154_1",
                "748_0",
                "748_1",
                "800_0",
                "800_1",
                "882_0",
                "882_1",
                "591_0",
                "591_1",
                "821_0",
                "821_1",
                "570_0",
                "570_1",
                "570_2",
                "108_0",
                "108_1",
                "188_0",
                "188_1",
                "409_0",
                "409_1",
                "409_2",
                "382_0",
                "382_1",
                "111_0",
                "111_1",
                "869_0",
                "869_1",
                "634_0",
                "634_1",
                "242_0",
                "242_1",
                "533_0",
                "533_1",
                "494_0",
                "494_1",
                "997_0",
                "997_1",
                "388_0",
                "388_1",
                "388_2",
                "636_0",
                "636_1",
                "636_2",
                "429_0",
                "429_1",
                "429_2",
                "446_0",
                "446_1",
                "450_0",
                "450_1",
                "527_0",
                "527_1",
                "758_0",
                "758_1",
                "758_2",
                "275_0",
                "275_1",
                "463_0",
                "463_1",
                "463_2",
                "701_0",
                "701_1",
                "959_0",
                "959_1",
                "412_0",
                "412_1",
                "020_0",
                "020_1",
                "453_0",
                "453_1",
                "961_0",
                "961_1",
                "009_0",
                "009_1",
                "543_0",
                "543_1",
                "543_2",
                "406_0",
                "406_1",
                "432_0",
                "432_1",
                "106_0",
                "106_1",
                "573_0",
                "573_1",
                "266_0",
                "266_1",
                "825_0",
                "825_1",
                "505_0",
                "505_1",
                "505_2",
                "750_0",
                "750_1",
                "530_0",
                "530_1",
                "605_0",
                "605_1",
                "875_0",
                "875_1",
                "589_0",
                "589_1",
                "541_0",
                "541_1",
                "488_0",
                "488_1",
                "387_0",
                "387_1",
                "464_0",
                "464_1",
                "464_2",
                "572_0",
                "572_1",
                "374_0",
                "374_1",
                "460_0",
                "460_1",
                "950_0",
                "950_1",
                "806_0",
                "806_1",
                "706_0",
                "706_1",
                "542_0",
                "542_1",
                "542_2",
                "499_0",
                "499_1",
                "499_2",
                "063_0",
                "063_1")

timelevels=c("0","1","2")

clustercolors <- as.character(c(cols25(n=25),alphabet(n=19)))

#metacluster heatmap
plot_clustering_heatmap_wrapper2(fcs=output$fcs,
                                 color_clusters = clustercolors,
                                 cell_clustering = factor(output$cell_clustering1m, levels=clusterlevels), 
                                 subtype_markers=output$subtype_markers,
                                 fileName = 'Extended_Aging_clusteringheatmap_final.pdf');dev.off()


#save output
#saveRDS(output, file = "backup_output.rds")

####DIFFERENTIAL PLOTS####
#set up count and prop matrices
counts_table <- table(output$cell_clustering1m, output$sample_ids)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)

write.csv(counts, file='Extended_Aging_counts.csv')
write.csv(props, file='Extended_Aging_props.csv')
