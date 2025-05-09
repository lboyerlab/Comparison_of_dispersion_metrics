library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library(tidyr) # need this to pivot tables
library(tidyverse)
library(extraDistr)
library(DescTools) # contains Gini & Entropy functions
library(reshape2)
library(pheatmap)
library(BPSC) # for the Beta-Poisson
library(patchwork) # alternative to grid.arrange
library(ggplotify)
library(scModels)

#### set up #### 
#' Loads in CSVs from simulations for all 4 metrics for a given family of distributions as data.frames
#' @param dist String name of the distribution used to sample
#' @param str_mod String, optional for naming the saved files 
#' @returns List of data.frames
load_rel_csvs <- function(dist, str_mod) { # this should return a list of 4 dataframes (one for each metric)
  gini_df <- read.csv(paste("Simulations/", dist, "_Gini_", str_mod, ".csv", sep = ""), header = T, row.names = 1, check.names = F)
  vmr_df<- read.csv(paste("Simulations/", dist, "_VMR_", str_mod, ".csv", sep = ""), header = T, row.names = 1, check.names = F)
  var_df <- read.csv(paste("Simulations/", dist, "_Var_", str_mod, ".csv", sep = ""), header = T, row.names = 1, check.names = F)
  entropy_df <- read.csv(paste("Simulations/", dist, "_Entropy_", str_mod, ".csv", sep = ""), header = T, row.names = 1, check.names = F)
  all_dfs <- list(gini_df, vmr_df, var_df, entropy_df)
  rev_dfs <- lapply(all_dfs, function(df) df %>% arrange(desc(row_number())))
  return(rev_dfs)
}

#' Loads in CSVs from simulations for all 4 metrics for a given family of distributions as data.frames, without reversing the direction of the rows.
#' This is used to make sure that the dispersion is increasing up the rows (rather than down the rows).
#' @param dist String name of the distribution used to sample
#' @param str_mod String, optional for naming the saved files 
#' @returns List of data.frames
load_wo_flip_csvs <- function(dist, str_mod) { # this should return a list of 4 dataframes (one for each metric)
  gini_df <- read.csv(paste("Simulations/", dist, "_Gini_", str_mod, ".csv", sep = ""), header = T, row.names = 1, check.names = F)
  vmr_df<- read.csv(paste("Simulations/", dist, "_VMR_", str_mod, ".csv", sep = ""), header = T, row.names = 1, check.names = F)
  var_df <- read.csv(paste("Simulations/", dist, "_Var_", str_mod, ".csv", sep = ""), header = T, row.names = 1, check.names = F)
  entropy_df <- read.csv(paste("Simulations/", dist, "_Entropy_", str_mod, ".csv", sep = ""), header = T, row.names = 1, check.names = F)
  all_dfs <- list(gini_df, vmr_df, var_df, entropy_df)
  return(all_dfs)
}

#' Loads in CSVs from simulations for all 4 metrics for a given family of distributions as data.frames, where the dispersion was fixed 
#' @param dist String name of the distribution used to sample
#' @returns List of data.frames
load_fix_dispersion_csvs <- function(dist) {
  gini_df <- read.csv(paste("Simulations/", dist, "_Gini_", "fix_dispersion.csv", sep = ""), header = T, row.names = 1, check.names = F)
  vmr_df<- read.csv(paste("Simulations/", dist, "_VMR_","fix_dispersion.csv", sep = ""), header = T, row.names = 1, check.names = F)
  var_df <- read.csv(paste("Simulations/", dist, "_Var_", "fix_dispersion.csv", sep = ""), header = T, row.names = 1, check.names = F)
  entropy_df <- read.csv(paste("Simulations/", dist, "_Entropy_", "fix_dispersion.csv", sep = ""), header = T, row.names = 1, check.names = F)
  all_dfs <- list(gini_df, vmr_df, var_df, entropy_df)
  rev_dfs <- lapply(all_dfs, function(df) df %>% tibble::rownames_to_column() %>% arrange(desc(row_number())))
  rev_dfs <- lapply(rev_dfs, fix_rownames)
  return(rev_dfs)
}

#' Changes data.frame rownames to match after reversing row order 
#' @param df data.frame
#' @returns data.frame with the rownames reversed
fix_rownames <- function(df) {
  out <- df[, -1]
  rownames(out) <- df[, 1]
  return(out)
}

#' Truncates the rownames of a data.frame to a specified number of significant figures
#' @param df data.frame
#' @param dig Number of significant figures to keep 
#' @returns data.frame with truncated rownames
truncate_df_rownames <- function(df, dig) {
  temp_rows <- as.numeric(rownames(df))
  round_rows <- sapply(temp_rows, signif, digits = dig)
  rownames(df) <- round_rows
  return(df)
}

## Function calls to load in CSVs generated from Simulations Clean.R ##
# Poisson
pois_ngene_x_ncell <- load_fix_dispersion_csvs("Poisson_2")
pois_ngene_x_dispersion <- lapply(load_rel_csvs("Poisson_3", "vary_ngenes"), truncate_df_rownames, 3)
pois_ncell_x_dispersion <- lapply(load_rel_csvs("Poisson_3", "vary_ncells"), truncate_df_rownames, 3)

# Negative Binomial
nb_ngene_x_dispersion <- lapply(load_wo_flip_csvs("Negative_Binomial_2", "vary_ngenes"), truncate_df_rownames, 4)
nb_ncell_x_dispersion <- lapply(load_wo_flip_csvs("Negative_Binomial_fix_prob", "vary_ncells"), truncate_df_rownames, 4)
nb_ngene_x_ncell <- load_fix_dispersion_csvs("Negative_Binomial_fix_prob")

# Beta Poisson
bp_ngene_x_dispersion <- lapply(load_rel_csvs("Beta-Poisson_3", "vary_ngenes"), truncate_df_rownames, 3)
bp_ncell_x_dispersion <- lapply(load_rel_csvs("Beta-Poisson_3", "vary_ncells"), truncate_df_rownames, 3)
bp_ngene_x_ncell <- lapply(load_fix_dispersion_csvs("Beta-Poisson"), truncate_df_rownames, 3)

# Uniform
unif_ngene_x_dispersion <- lapply(load_rel_csvs("Uniform_2", "vary_ngenes"), truncate_df_rownames, 3)
unif_ncell_x_dispersion <- lapply(load_rel_csvs("Uniform", "vary_ncells"), truncate_df_rownames, 3)
unif_ngene_x_ncell <- load_fix_dispersion_csvs("Uniform")

# jitter the Entropy values for uniform (they're all the same) so that we can plot heatmaps 
unif_ngene_x_dispersion[[4]] <- apply(unif_ngene_x_dispersion[[4]], c(1, 2), function(x) x + runif(1, min = 0, max = 0.00001))
unif_ncell_x_dispersion[[4]] <- apply(unif_ncell_x_dispersion[[4]], c(1, 2), function(x) x + runif(1, min = 0, max = 0.00001))

#### make the heatmaps ####
## set up plotting functions ## 
#' Generates a pheatmap
#' @param x data.frame of values to generate the heatmap
#' @param title String title of the heatmap
#' @returns pheatmap
plot_heatmap <- function(x, title) {
  rnames <- rownames(as.matrix(x))
  rnames[seq_along(rnames) %% 5 != 1] <- ""
  
  cnames <- colnames(as.matrix(x))
  cnames[seq_along(cnames) %% 5 != 1] <- ""
  
  hm <- pheatmap(as.matrix(x), main = title, cluster_rows = F, cluster_cols = F, border_color = NA, 
                 angle_col = 90, fontsize = 7, labels_row = as.character(rnames), labels_col = as.character(cnames), 
                 silent = T)
  return(hm)
}

#' Generates a list of 4 heatmaps corresponding to each of the 4 metrics
#' @param df_list list of data.frames from which to generate the heatmaps
#' @param title String title of the heatmap
#' @param x_title String title for x-axis
#' @param y_title String title for y-axis
#' @returns list of pheatmap objects
plot_raw_heatmaps <- function(df_list, title, x_title, y_title) { # NOTE: the order of the df_list has to match the order here
  gini_hm <- plot_heatmap(df_list[[1]], "Gini")
  vmr_hm <- plot_heatmap(df_list[[2]], "VMR")
  var_hm <- plot_heatmap(df_list[[3]], "Variance")
  entropy_hm <- plot_heatmap(df_list[[4]], "Entropy")
  p <- grid.arrange(gini_hm[[4]], vmr_hm[[4]], var_hm[[4]], entropy_hm[[4]], ncol = 2, top = textGrob(title), right = y_title, bottom = x_title, padding = unit(1, "line"))
  return(list(gini_hm, vmr_hm, var_hm, entropy_hm))
}

## function calls to generate heatmaps of the raw, unnormalized dispersion metrics ##
pois_raw_gene <- plot_raw_heatmaps(pois_ngene_x_dispersion, "Poisson: Number of genes vs Lambda", "Number of genes", "Lambda")
pois_raw_cell <- plot_raw_heatmaps(pois_ncell_x_dispersion, "Poisson: Number of cells vs Lambda", "Number of cells", "Lambda")

nb_raw_gene <- plot_raw_heatmaps(nb_ngene_x_dispersion, "Negative Binomial: Number of genes vs r", "Number of genes", "r")
nb_raw_cell <- plot_raw_heatmaps(nb_ncell_x_dispersion, "Negative Binomial: Number of cells vs r", "Number of cells", "r")

unif_raw_gene <- plot_raw_heatmaps(unif_ngene_x_dispersion, "Uniform: Number of genes vs max", "Number of genes", "Max")
unif_raw_cell <- plot_raw_heatmaps(unif_ncell_x_dispersion, "Uniform: Number of cells vs max", "Number of cells", "Max")

bp_raw_gene <- plot_raw_heatmaps(bp_ngene_x_dispersion, "Beta-Poisson: Number of genes vs scale factor", "Number of genes", "Scale factor")
bp_raw_cell <- plot_raw_heatmaps(bp_ncell_x_dispersion, "Beta-Poisson: Number of cells vs scale factor", "Number of cells", "Scale factor")

## generate the relative change heatmaps ## 
#' Applies relative change to normalize down a matrix
#' @param mat Matrix of values to normalize, normalize a given value against the value above it in the same column
#' @returns Normalized matrix
metric_rel_change <- function(mat) {
  return(apply(mat, 2, function(x) -diff(x) / (tail(x, -1) + 1e-10)))
}

#' Renames a data.frame
#' @param df data.frame
#' @param new_rownames Vector of rownames to use as new rownames 
#' @returns data.frame with updated rownames
update_rownames <- function(df, new_rownames) {rownames(df) <- new_rownames; return(df)}

# calculate the relative change in the metrics
pois_metric_rel_change_ngene <- lapply(pois_ngene_x_dispersion, metric_rel_change)
pois_metric_rel_change_ncell <- lapply(pois_ncell_x_dispersion, metric_rel_change)

nb_metric_rel_change_ngene <- lapply(nb_ngene_x_dispersion, metric_rel_change)
nb_metric_rel_change_ncell <- lapply(nb_ncell_x_dispersion, metric_rel_change)

unif_metric_rel_change_ngene <- lapply(unif_ngene_x_dispersion, metric_rel_change)
unif_metric_rel_change_ncell <- lapply(unif_ncell_x_dispersion, metric_rel_change)

bp_metric_rel_change_ngene <- lapply(bp_ngene_x_dispersion, metric_rel_change)
bp_metric_rel_change_ncell <- lapply(bp_ncell_x_dispersion, metric_rel_change)

# calculate the relative change in the dispersion (used for the labels on the heatmap)
# Negative Binomial
prob <- 0.5
nb_disp_fix_prob <- as.numeric(rownames(nb_ngene_x_dispersion[[1]]))*(1 - prob)/prob^2
nb_metric_rel_change_ngene <- lapply(nb_metric_rel_change_ngene, update_rownames, head(nb_disp_fix_prob, -1))
nb_metric_rel_change_ncell <- lapply(nb_metric_rel_change_ncell, update_rownames, head(nb_disp_fix_prob, -1))

# Uniform 
unif_disp <- sapply(sapply(as.numeric(rownames(unif_ngene_x_dispersion[[1]])), function(x) 1/12 * x**2), signif, 3)
unif_metric_rel_change_ngene <- lapply(unif_metric_rel_change_ngene, update_rownames, head(unif_disp, -1))
unif_metric_rel_change_ncell <- lapply(unif_metric_rel_change_ncell, update_rownames, head(unif_disp, -1))

# Beta-Poisson 
lam1 <- 50  * 10
lam2 <- 50 * 0.01
alp <- 50 * 0.1
bp_var <- sapply(sapply(as.numeric(rownames(bp_ncell_x_dispersion[[1]])), function(x) varBP(alp, 50 * 0.1*x, lam1 = lam1, lam2 = lam2)), signif, 3)

bp_metric_rel_change_ngene <- lapply(bp_metric_rel_change_ngene, update_rownames, head(bp_var, -1))
bp_metric_rel_change_ncell <- lapply(bp_metric_rel_change_ncell, update_rownames, head(bp_var, -1))

# generate heatmaps
pois_rel_change_gene <- plot_raw_heatmaps(pois_metric_rel_change_ngene, "Poisson: Relative Change", "Number of genes", "Theoretical variance")
pois_rel_change_cell <- plot_raw_heatmaps(pois_metric_rel_change_ncell, "Poisson: Relative Change", "Number of cells", "Theoretical variance")

nb_rel_gene <- plot_raw_heatmaps(nb_metric_rel_change_ngene, "Negative Binomial: Relative Change", "Number of genes", "Theoretical variance")
nb_rel_cell <- plot_raw_heatmaps(nb_metric_rel_change_ncell, "Negative Binomial: Relative Change", "Number of cells", "Theoretical variance")

bp_rel_gene <- plot_raw_heatmaps(bp_metric_rel_change_ngene, "Beta-Poisson: Relative Change", "Number of genes", "Theoretical variance")
bp_rel_cell <- plot_raw_heatmaps(bp_metric_rel_change_ncell, "Beta-Poisson: Relative Change", "Number of cells", "Theoretical variance")

unif_rel_gene <- plot_raw_heatmaps(unif_metric_rel_change_ngene, "Uniform: Relative Change", "Number of genes", "Theoretical variance")
unif_rel_cell <- plot_raw_heatmaps(unif_metric_rel_change_ncell, "Uniform: Relative Change", "Number of cells", "Theoretical variance")

## generate ratio change heatmaps ## 
#' Applies ratio change to normalize down a matrix
#' @param mat Matrix of values to normalize, normalize a given value against the value above it in the same column
#' @returns Normalized matrix
metric_ratio_change <- function(mat) {return(apply(mat, 2, function(x) head(x, -1) / (tail(x, -1) + 1e-10)))}

# calculate the ratio change in the metrics
pois_metric_ratio_change_ngene <- lapply(pois_ngene_x_dispersion, metric_ratio_change)
pois_metric_ratio_change_ncell <- lapply(pois_ncell_x_dispersion, metric_ratio_change)

nb_metric_ratio_change_ngene <- lapply(nb_ngene_x_dispersion, metric_ratio_change)
nb_metric_ratio_change_ncell <- lapply(nb_ncell_x_dispersion, metric_ratio_change)

unif_metric_ratio_change_ngene <- lapply(unif_ngene_x_dispersion, metric_ratio_change)
unif_metric_ratio_change_ncell <- lapply(unif_ncell_x_dispersion, metric_ratio_change)

bp_metric_ratio_change_ngene <- lapply(bp_ngene_x_dispersion, metric_ratio_change)
bp_metric_ratio_change_ncell <- lapply(bp_ncell_x_dispersion, metric_ratio_change)

# calculate ratio change in the theoretical variance (for labeling the heatmaps)
nb_metric_ratio_change_ngene <- lapply(nb_metric_ratio_change_ngene, update_rownames, head(nb_disp_fix_prob, -1))
nb_metric_ratio_change_ncell <- lapply(nb_metric_ratio_change_ncell, update_rownames, head(nb_disp_fix_prob, -1))

unif_metric_ratio_change_ngene <- lapply(unif_metric_ratio_change_ngene, update_rownames, head(unif_disp, -1))
unif_metric_ratio_change_ncell <- lapply(unif_metric_ratio_change_ncell, update_rownames, head(unif_disp, -1))

bp_metric_ratio_change_ngene <- lapply(bp_metric_ratio_change_ngene, update_rownames, head(bp_var, -1))
bp_metric_ratio_change_ncell <- lapply(bp_metric_ratio_change_ncell, update_rownames, head(bp_var, -1))

# generate heatmaps
pois_ratio_gene <- plot_raw_heatmaps(pois_metric_ratio_change_ngene, "Poisson: Ratio Change", "Number of genes", "Theoretical variance")
pois_ratio_cell <- plot_raw_heatmaps(pois_metric_ratio_change_ncell, "Poisson: Ratio Change", "Number of cells", "Theoretical variance")

nb_ratio_gene <- plot_raw_heatmaps(nb_metric_ratio_change_ngene, "Negative Binomial: Ratio Change", "Number of genes", "Theoretical variance")
nb_ratio_cell <- plot_raw_heatmaps(nb_metric_ratio_change_ncell, "Negative Binomial: Ratio Change", "Number of cells", "Theoretical variance")

bp_ratio_gene <- plot_raw_heatmaps(bp_metric_ratio_change_ngene, "Beta-Poisson: Ratio Change", "Number of genes", "Theoretical variance")
bp_ratio_cell <- plot_raw_heatmaps(bp_metric_ratio_change_ncell, "Beta-Poisson: Ratio Change", "Number of cells", "Theoretical variance")

unif_ratio_gene <- plot_raw_heatmaps(unif_metric_ratio_change_ngene, "Uniform: Ratio Change", "Number of genes", "Theoretical variance")
unif_ratio_cell <- plot_raw_heatmaps(unif_metric_ratio_change_ncell, "Uniform: Ratio Change", "Number of cells", "Theoretical variance")

## generate heatmaps for fixed dispersion simulations ##
pois_fix_disp <- plot_raw_heatmaps(pois_ngene_x_ncell, "Poisson: Number of genes vs Number of cells", "Number of genes", "Number of cells")
nb_fix_disp <- plot_raw_heatmaps(nb_ngene_x_ncell, "Negative Binomial: Number of genes vs Number of cells", "Number of genes", "Number of cells")
unif_fix_disp <- plot_raw_heatmaps(unif_ngene_x_ncell, "Uniform: Number of genes vs Number of cells", "Number of genes", "Number of cells")
bp_fix_disp <- plot_raw_heatmaps(bp_ngene_x_ncell, "Beta-Poisson: Number of genes vs Number of cells", "Number of genes", "Number of cells")

## Additional plotting for figure generation ##
#' Plots heatmaps of raw, relative change normalized, and ratio change normalized values
#' @param raw List of pheatmaps for the raw metric values
#' @param rel List of pheatmaps for the relative change normalized metric values
#' @param ratio List of pheatmaps for the ratio change normalized metric values
plot_raw_and_normalized <- function(raw, rel, ratio) {
  (as.ggplot(raw[[1]]) | as.ggplot(raw[[2]]) | as.ggplot(rel[[1]]) | as.ggplot(rel[[2]]) | as.ggplot(ratio[[1]]) | as.ggplot(ratio[[2]])) /
    (as.ggplot(raw[[3]]) | as.ggplot(raw[[4]]) | as.ggplot(rel[[3]]) | as.ggplot(rel[[4]]) | as.ggplot(ratio[[3]]) | as.ggplot(ratio[[4]])) & 
    theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))
}

# function calls to plot the raw, relative change normalized, and ratio change normalized values together 
plot_raw_and_normalized(pois_raw_gene, pois_rel_change_gene, pois_ratio_gene)
plot_raw_and_normalized(pois_raw_cell, pois_rel_change_cell, pois_ratio_cell)

plot_raw_and_normalized(nb_raw_gene, nb_rel_gene, nb_ratio_gene)
plot_raw_and_normalized(nb_raw_cell, nb_rel_cell, nb_ratio_cell)

plot_raw_and_normalized(unif_raw_gene, unif_rel_gene, unif_ratio_gene)
plot_raw_and_normalized(unif_raw_cell, unif_rel_cell, unif_ratio_cell)

plot_raw_and_normalized(bp_raw_gene, bp_rel_gene, bp_ratio_gene)
plot_raw_and_normalized(bp_raw_cell, bp_rel_cell, bp_ratio_cell)

# plot the fixed dispersion heatmaps together (broken up into 2 figures)
(as.ggplot(pois_fix_disp[[1]]) | as.ggplot(pois_fix_disp[[2]]) | as.ggplot(nb_fix_disp[[1]]) | as.ggplot(nb_fix_disp[[2]])) /
  (as.ggplot(pois_fix_disp[[3]]) | as.ggplot(pois_fix_disp[[4]]) | as.ggplot(nb_fix_disp[[3]]) | as.ggplot(nb_fix_disp[[4]])) &
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))
(as.ggplot(bp_fix_disp[[1]]) | as.ggplot(bp_fix_disp[[2]]) | as.ggplot(unif_fix_disp[[1]]) | as.ggplot(unif_fix_disp[[2]])) /
  (as.ggplot(bp_fix_disp[[3]]) | as.ggplot(bp_fix_disp[[4]]) | as.ggplot(unif_fix_disp[[3]]) | as.ggplot(unif_fix_disp[[4]])) & 
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))

# Scatter plots 
#' Plots scatter plots of raw metric values vs theoretical variance 
#' @param raw_mean List of pheatmaps for the raw metric values
#' @param theoretical_var List of pheatmaps for the relative change normalized metric values
raw_scatter_plot <- function(raw_mean, theoretical_var) {
  raw_mean_dfs <- lapply(raw_mean, function(x) data.frame(theoretical_var = as.numeric(theoretical_var), raw_mean = x))
  raw_scatter_plots <- lapply(raw_mean_dfs, function(x) ggplot(x) + 
                                geom_point(aes(x = theoretical_var, y = raw_mean)) + 
                                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                      panel.background = element_blank(), axis.line = element_line(colour = "black")) +  
                                labs(x = "Theoretical Variance", y = "Dispersion Metric"))
  (raw_scatter_plots[[1]] + ggtitle("Gini") + raw_scatter_plots[[2]] + ggtitle("VMR"))/
    (raw_scatter_plots[[3]] + ggtitle("Variance") + raw_scatter_plots[[4]] + ggtitle("Entropy"))
}

nb_raw_mean <- lapply(nb_ncell_x_dispersion, rowMeans)
raw_scatter_plot(nb_raw_mean, nb_disp_fix_prob)

#### plot PDF of distributions for visualization in Fig. 1 #### 
poisson_df <- expand.grid(lambda = as.numeric(rownames(pois_ngene_x_dispersion[[1]]))[c(T, F, F, F, F, F, F, F, F)], x = 0:100) %>%
  mutate(prob = dpois(x, lambda))
pois_p <- ggplot(poisson_df, aes(x = x, y = prob, color = as.factor(lambda))) +
  geom_line() +
  labs(
    x = "x", y = "Probability",
    color = "Lambda") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

unif_df <- expand.grid(max = as.numeric(rownames(unif_ngene_x_dispersion[[1]]))[c(T, F, F, F, F, F, F, F, F)], x = 0:100) %>%
  mutate(prob = dunif(x, 0, max))
unif_p <- ggplot(unif_df, aes(x = x, y = prob, color = as.factor(max))) +
  geom_line() +
  labs(
    x = "x", y = "Probability",
    color = "Max") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

nb_df <- expand.grid(size = as.numeric(rownames(nb_ngene_x_dispersion[[1]]))[c(T, F, F, F, F, F, F, F, F)], x = seq(0, 25, 1)) %>%
  mutate(p = dnbinom(x, size, prob = 0.5))
nb_p <- ggplot(nb_df, aes(x = x, y = p, color = as.factor(size))) +
  geom_line() +
  labs(
    x = "x", y = "Probability",
    color = "Size") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

bp_df <-  expand.grid(scale = as.numeric(rownames(bp_ngene_x_dispersion[[1]]))[c(T, F, F, F, F, F, F, F, F)], x = 0:250)
bp_df$prob <- apply(bp_df, 1, function(x) dpb(x[2], 5, x[1]*5, 250))
bp_p <- ggplot(bp_df, aes(x = x, y = prob, color = as.factor(scale))) +
  geom_line() +
  labs(x = "x", y = "Probability", color = "Scale factor") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# put plots together
pois_p | nb_p | bp_p | unif_p

# generate example distribution in Fig. 1b
example_dist_df <- expand.grid(size = c(8), x = seq(0, 50, 1)) %>%
  mutate(prob = dnbinom(x, size, prob = 0.5))
ggplot(example_dist_df, aes(x = x, y = prob)) +
  geom_line() + labs(x = "", y = "") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

## Generate figures for paradoxical behavior of the Gini (Fig. 3) ##
# PDFs 
# Case 1: many 0s, size < 1
zeros <- data.frame(x = 0:10, prob = dnbinom(0:10, 1, 0.5))
pdf_1 <- ggplot(zeros) + geom_line(aes(x = x, y = prob)) + 
  labs(y = "Probability Density", x = "x") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Case 2: size > 1
not_zeros <- data.frame(x = 0:30, prob = dnbinom(0:30, 10, 0.5))
pdf_2 <-ggplot(not_zeros) + geom_line(aes(x = x, y = prob)) + 
  labs(y = "Probability Density", x = "x") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

pdf_1 + pdf_2

# Lorenz Curves
# Case 1
zeros_sample <- rnbinom(1000, 1, 0.5)
plot(Lc(zeros_sample), col = "blue", lwd = 2, main = "", xlab = "Cumulative share of people, ranked by increasing income", ylab = "Cumulative share of income")
abline(0, 1, col = "red", lwd = 2, lty = 2)  # Dashed red line

# Case 2
not_zeros_sample <- rnbinom(1000, 10, 0.5)
plot(Lc(not_zeros_sample), col = "blue", lwd = 2, main = "", xlab = "Cumulative share of people, ranked by increasing income", ylab = "Cumulative share of income")
abline(0, 1, col = "red", lwd = 2, lty = 2)  # Dashed red line



