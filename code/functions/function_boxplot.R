# Custom function for boxplot

# Libraries

library(here)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(stats)
library(gridExtra)
library(grid)
library(ggbreak) 
library(ggsci)
library(viridis)

here::i_am("code/functions/function_boxplot.R")

# cytokine: two color input color list fxn

two.color.input_cyt_boxplot_fxn <- function(dataset, conc, factors, cyt_list, facet_names, pal_list) {
  color_pal <- paste(pal_list)
  x <- paste(factors)
  y <- paste(conc)
  dataset_sig <- subset(dataset, dataset$cytokines %in% cyt_list)
  dataset_sig <- dataset_sig %>% arrange(desc(`pt id`), desc(cytokines))
  
  plot <- ggplot(dataset_sig, aes(x=as.factor(.data[[x]]), y= as.numeric(.data[[y]]), 
                                  color = as.factor(.data[[x]]))) + 
    scale_color_manual(values = color_pal) + 
    geom_boxplot(size = 0.75, 
                 fill = 'white'
    ) +
    geom_jitter(alpha = 0.75, pch=19, size=1.5, position=position_jitterdodge(dodge.width=1, jitter.width=1)) +
    theme(plot.title = element_text(hjust = 0.5, size = 16), axis.text.y.left = element_text(angle = 0), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), panel.background = element_blank(), legend.position = 'none', axis.text.y = element_text(color="black", size=11, angle = 90, hjust = 0.5), axis.text.x = element_text(color="black",size=11, angle = 0), axis.line = element_line(colour = "black"), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
          strip.text = element_text(#face = "bold", 
            size = 12)) + 
    facet_wrap(~ `cytokines`, scales = "free", labeller = labeller(cytokines = facet_names))
  
  return(plot)
}

# cytokine: multi-color input list fxn

multi_color_boxplot_fxn <- function(dataset, conc, factors, cyt_list, facet_names) {
  colors <- c(pal_npg()(10))
  color_pal <- c(colors[2], colors[1], colors[3:10])
  x <- paste(factors)
  y <- paste(conc)
  dataset_sig <- subset(dataset, dataset$cytokines %in% cyt_list)
  dataset_sig <- dataset_sig %>% arrange(desc(`pt id`), desc(cytokines))
  
  plot <- ggplot(dataset_sig, aes(x=as.factor(.data[[x]]), y= as.numeric(.data[[y]]), 
                                  color = as.factor(.data[[x]]))) + 
    scale_color_manual(values = color_pal) + 
    geom_boxplot(size = 0.75, 
                 fill = 'white'
    ) +
    geom_jitter(alpha = 0.75, pch=19, size=1.5, position=position_jitterdodge(dodge.width=1, jitter.width=1.25)) +
    theme(plot.title = element_text(hjust = 0.5, size = 16), axis.text.y.left = element_text(angle = 0), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), panel.background = element_blank(), legend.position = 'none', axis.text.y = element_text(color="black", size=11, angle = 90, hjust = 0.5), axis.text.x = element_text(color="black",size=11, angle = 0), axis.line = element_line(colour = "black"), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
          strip.text = element_text(size = 12)) + 
    facet_wrap(~ `cytokines`, scales = "free", labeller = labeller(cytokines = facet_names))
  
  return(plot)
}

# cytokine: multi-color input with two grouping variables

multi_color_two.grouping.variables_boxplot_fxn <- function(dataset, conc, factors, secondary.factor, cyt_list, facet_names) {
  colors <- c(pal_npg()(10))
  color_pal <- c(colors[2], colors[1], colors[3:10])
  x <- paste(factors)
  y <- paste(conc)
  z <- paste(secondary.factor)
  dataset_sig <- subset(dataset, dataset$cytokines %in% cyt_list)
  dataset_sig <- dataset_sig %>% arrange(desc(`pt id`), desc(cytokines))
  
  plot <- ggplot(dataset_sig, aes(x=as.factor(.data[[x]]), y= as.numeric(.data[[y]]), 
                                  color = as.factor(.data[[z]]))) + 
    scale_color_manual(values = color_pal) + 
    geom_boxplot(size = 0.75, 
                 fill = 'white'
    ) +
    geom_jitter(alpha = 0.75, pch=19, size=1.5, position=position_jitterdodge(dodge.width=1, jitter.width=1.25)) +
    theme(plot.title = element_text(hjust = 0.5, size = 16), axis.text.y.left = element_text(angle = 0), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), panel.background = element_blank(), legend.position = 'none', axis.text.y = element_text(color="black", size=11, angle = 90, hjust = 0.5), axis.text.x = element_text(color="black",size=11, angle = 0), axis.line = element_line(colour = "black"), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
          strip.text = element_text(#face = "bold", 
            size = 12)) + 
    facet_wrap(~ `cytokines`, scales = "free", labeller = labeller(cytokines = facet_names))
  
  return(plot)
}

# cytokine: boxplot function that selects for significant cytokines (for comparisons of >=3)

boxplot_sig.cyt_kruskal_fxn <- function(dataset, conc, factors) {
  formula <- as.formula(paste0(conc, " ~ ", factors))
  colors <- c(pal_npg()(10))
  color_pal <- c(colors[2], colors[1], colors[3:10])
  list_comparisons_all <- compare_means(formula, group.by = 'cytokines', method = 'kruskal.test', data = dataset, p.adjust.method = 'fdr')
  x <- paste(factors)
  y <- paste(conc)
  sig_list_comparisons <- subset(list_comparisons_all,p <0.05)
  dataset_sig <- subset(dataset, dataset$cytokines %in% sig_list_comparisons$cytokines)
  dataset_sig <- dataset_sig %>% arrange(desc(`pt id`), desc(cytokines))
  
  baseline_comparison <- ggplot(dataset_sig, aes(x=as.factor(.data[[x]]), y= as.numeric(.data[[y]]), col = as.factor(.data[[x]]))) + 
    scale_color_manual(values = color_pal) +
    geom_boxplot(fill = 'white') +
    geom_jitter(alpha = 0.75, pch=19, size=1.5, position=position_jitterdodge(dodge.width=1, jitter.width=1.5)) + 
    theme(plot.title = element_text(hjust = 0.5, size = 16), axis.text.y.left = element_text(angle = 0), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), panel.background = element_blank(), legend.position = 'none', axis.text.y = element_text(color="black", size=11, angle = 90, hjust = 0.5), axis.text.x = element_text(color="black",size=11, angle = 0), axis.line = element_line(colour = "black"), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + 
    facet_wrap(~ `cytokines`, scales = "free")
  
  return(list(baseline_comparison,sig_list_comparisons))
}

# Two color boxplot function with cytof_list (input colors)

two.color.input_cytof_boxplot_fxn <- function(dataset, conc, factors, cytof_list, facet_names, pal_list) {
  color_pal <- paste(pal_list)
  x <- paste(factors)
  y <- paste(conc)
  dataset_sig <- subset(dataset, dataset$immune_cell %in% cytof_list)
  dataset_sig <- dataset_sig %>% arrange(desc(`pt id`), desc(immune_cell))
  
  plot <- ggplot(dataset_sig, aes(x=as.factor(.data[[x]]), y= as.numeric(.data[[y]]), 
                                  color = as.factor(.data[[x]]))) + 
    scale_color_manual(values = color_pal) + 
    geom_boxplot(size = 0.75, 
                 fill = 'white'
    ) +
    geom_jitter(alpha = 0.75, pch=19, size=1.5, position=position_jitterdodge(dodge.width=1, jitter.width=1)) +
    theme(plot.title = element_text(hjust = 0.5, size = 16), axis.text.y.left = element_text(angle = 0), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), panel.background = element_blank(), legend.position = 'none', axis.text.y = element_text(color="black", size=11, angle = 90, hjust = 0.5), axis.text.x = element_text(color="black",size=11, angle = 0), axis.line = element_line(colour = "black"), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
          strip.text = element_text(size = 12)) + 
    facet_wrap(~ `immune_cell`, scales = "free", labeller = labeller(immune_cell = facet_names))
  
  return(plot)
}

