library(dplyr)
library(tidyr)
library(devtools)
library(BPSC)
library(DescTools)

#### set up ####
path_to_data <- "~/Simulations/"

n_genes <- seq(1000, 10800, 200)
n_cells <- seq(500, 5400, 100)

vmr <- function(x) var(x)/mean(x)

#' Converts matrix to dataframe using given row and column names
#' @param mat The matrix to be converted
#' @param rnames Vector of row names 
#' @param cnames Vector of column names 
#' @return data.frame 
mat_to_df <- function(mat, rnames, cnames) {
  as_df <- as.data.frame(mat)
  rownames(as_df) <- rnames
  colnames(as_df) <- cnames
  return(as_df)
}

## set up Poisson parameters ##
lambdas <- numeric(100)
lambdas[1] <- 0.2
for (i in 2:length(lambdas)) {
  lambdas[i] <- lambdas[i - 1] * 1.075
}

#' Random generation from a Poisson distribution
#' @param n The number of samples to generate
#' @param lambda Lambda parameter of the Poisson
#' @param nr Integer, number of rows in the matrix to return 
#' @return Matrix with nr rows 
sample_pois <- function(n, lambda, nr) {
  data <- as.numeric(rpois(n, lambda))
  return(matrix(data, nrow = nr))
}

## set up negative binomial parameters ##
#' Random generation from a negative binomial with probability 0.5
#' @param n The number of samples to generate
#' @param size Target for number of successful trials
#' @param nr Integer, number of rows in the matrix to return 
#' @return Matrix with nr rows
sample_rnbinom <- function(n, size, nr) {
  data <- rnbinom(n, size, prob = 0.5)
  return(matrix(data, nrow = nr))
}

sizes <- numeric(100) 
sizes[1] <- 10
for (i in 2:length(sizes)) {
  sizes[i] <- 0.95*sizes[i - 1]
}

## set up beta-Poisson parameters ##
#' Random generation of a matrix from a beta-Poisson distribution using BPSC package (Vu et al., 2016)
#' @param n Number of samples to generate
#' @param scale_factor Scale factor on the beta parameter, changes the dispersion in the distribution 
#' @param n_genes Number of genes to simulate 
#' @return Matrix with n_genes rows
generate_betapois <- function(n, scale_factor, n_genes) { 
  alps <- sample(100, n_genes, replace = TRUE) * 0.1;
  bets <- sample(100, n_genes, replace = TRUE) * 0.1 * scale_factor;
  lam1s <- sample(100, n_genes, replace = TRUE) * 10;
  lam2s <-  sample(100, n_genes, replace = TRUE) * 0.01
  mat <- NULL
  n_cells <- n / n_genes
  for (i in 1:n_genes) mat <- cbind(mat, rBP(n_cells, alp = alps[i], bet = bets[i],
                                             lam1 = lam1s[i],lam2 = lam2s[i]))
  return(mat)
}

lam1 <- 50  * 10
lam2 <- 50 * 0.01
alp <- 50 * 0.1

scale_factors <- numeric(100)
scale_factors[1] <- 100
k <- 0.0684
for (i in 2:100) { # since relative increase in theoretical dispersion is fixed, need to calculate the needed scale factors
  RHS <- (k + 1)/(lam1 * (lam2**2)) * varBP(5, 5*scale_factors[i - 1], 500, 0.5)
  a <- RHS 
  b <- RHS*(3*alp + 1) - alp
  c <- RHS*(3*alp**2 + 2*alp) - 2*alp**2 - alp - lam1*alp
  d <- RHS*(alp**3 + alp**2) - alp**3 - alp**2
  
  coefficients <- c(d, c, b, a)
  roots <- sort(polyroot(coefficients), decreasing = T)
  for (r in Re(roots)) {
    if (r >= 0) {
      rel_change <- (varBP(5, 5*r/5, 500, 0.5) - varBP(5, 5*scale_factors[i - 1], 500, 0.5))/varBP(5, 5*scale_factors[i - 1], 500, 0.5)
      if (abs(rel_change - k) < 1e-6) {
        scale_factors[i] <- Re(r/5)
        break
      }
    }
  }
}

# set up uniform parameters 
bs <- numeric(99)
bs[1] <- 1
for (i in 2:length(bs)) {
  bs[i] <- bs[i - 1]*(sqrt(1.1))
}

#' Random generation from a uniform distribution
#' @param n The number of samples to generate
#' @param b Max parameter of the uniform
#' @param nr Integer, number of rows in the matrix to return 
#' @return Matrix with nr rows 
sample_unif <- function(n, b, nr) {
  data <- runif(n, max = b)
  return(matrix(data, nrow = nr))
}

#' Generates & saves noise matrices for all 4 metrics
#' @param dispersion_params Vector of dispersion parameters
#' @param dist_func Function from which to sample 
#' @param dist_str String name of the family of distributions used to sample, for naming the saved files 
#' @param str_mod String, optional for naming the saved files 
#' @param ncells Vector or integer of number of cells
#' @param ngenes Vector or integer of number of genes
#' @param is_ncell_change Logical; set to T if changing the number of cells, F is changing the number of genes 
generate_mats <- function(dispersion_params, dist_func, dist_str, str_mod = "", ncells, ngenes, is_ncell_change) {
  # set up matrices to save results
  if (is_ncell_change) {
    changing_n <- ncells
    fixed_n <- ngenes
  } else {
    changing_n <- ngenes
    fixed_n <- ncells
  }
  n_d <- length(dispersion_params)
  n_n <- length(changing_n)

  gini_mat <- matrix(NA, nrow = n_d, ncol = n_n)
  var_mat <- matrix(NA, nrow = n_d, ncol = n_n)
  vmr_mat <- matrix(NA, nrow = n_d, ncol = n_n)
  entropy_mat <- matrix(NA, nrow = n_d, ncol = n_n)
  
  for (i in 1:length(dispersion_params)) {
    d <- dispersion_params[i]
    
    for (j in 1:length(changing_n)) { # length(n_genes)
      n <- changing_n[j]
      
      # generate the counts 
      if (is_ncell_change) {nr <- ngenes} else {nr <- n}
      counts_mat <- dist_func(fixed_n * n, d, nr)

      # calculate the noise, take the median, and save into the matrices
      noise <- apply(counts_mat, 1, function(x) c(gini = Gini(x), variance = var(x), entropy = Entropy(table(x)), vmr = vmr(x)))
      noise_meds <- apply(noise, 1, median, na.rm = T)
      gini_mat[i, j] <- noise_meds[1]
      var_mat[i, j] <- noise_meds[2]
      entropy_mat[i, j] <- noise_meds[3]
      vmr_mat[i, j] <- noise_meds[4]
    }
    print(paste("round", i, "done"))
  }
  
  gini_df <- mat_to_df(gini_mat, dispersion_params, changing_n)
  var_df <- mat_to_df(var_mat, dispersion_params, changing_n)
  vmr_df <- mat_to_df(vmr_mat, dispersion_params, changing_n)
  entropy_df <- mat_to_df(entropy_mat, dispersion_params, changing_n)
  
  write.csv(gini_df, paste(path_to_data, dist_str, "_Gini_", str_mod, ".csv", sep = ""))
  write.csv(var_df, paste(path_to_data, dist_str, "_Var_", str_mod, ".csv", sep = ""))
  write.csv(vmr_df, paste(path_to_data, dist_str, "_VMR_", str_mod, ".csv", sep = ""))
  write.csv(entropy_df, paste(path_to_data, dist_str, "_Entropy_", str_mod, ".csv", sep = ""))
}

#' Generates & saves noise matrices for all 4 metrics where the dispersion is fixed, and the number of cells and number of genes varies
#' @param dist_func Function from which to sample 
#' @param dist_str String name of the family of distributions used to sample, for naming the saved files 
#' @param str_mod String, optional for naming the saved files 
#' @param ncells Vector or integer of number of cells
#' @param ngenes Vector or integer of number of genes
#' @param d Dispersion parameter to be used for all matrices 
generate_mats_fix_dispersion <- function(dist_func, dist_str, str_mod, n_cells, n_genes, d) {
  n_g <- length(n_genes)
  n_c <- length(n_cells)
  gini_mat <- matrix(NA, nrow = n_g, ncol = n_c)
  var_mat <- matrix(NA, nrow = n_g, ncol = n_c)
  vmr_mat <- matrix(NA, nrow = n_g, ncol = n_c)
  entropy_mat <- matrix(NA, nrow = n_g, ncol = n_c)
  
  for (i in 1:length(n_cells)) {
    n_cell <- n_cells[i]
    for (j in 1:length(n_genes)) { # length(n_genes)
      n_gene <- n_genes[j]
      
      # generate the counts 
      counts_mat <- dist_func(n_cell * n_gene, d, n_gene)
      
      # calculate the noise, take the median, and save into the matrices
      noise <- apply(counts_mat, 1, function(x) c(gini = Gini(x), variance = var(x), entropy = Entropy(table(x)), vmr = vmr(x)))
      noise_meds <- apply(noise, 1, median, na.rm = T)
      gini_mat[i, j] <- noise_meds[1]
      var_mat[i, j] <- noise_meds[2]
      entropy_mat[i, j] <- noise_meds[3]
      vmr_mat[i, j] <- noise_meds[4]
    }
    print(paste("round", i, "done"))
  }
  
  gini_df <- mat_to_df(gini_mat, n_cells, n_genes)
  var_df <- mat_to_df(var_mat, n_cells, n_genes)
  vmr_df <- mat_to_df(vmr_mat, n_cells, n_genes)
  entropy_df <- mat_to_df(entropy_mat, n_cells, n_genes)
  
  write.csv(gini_df, paste(path_to_data, dist_str, "_Gini_", str_mod, ".csv", sep = ""))
  write.csv(var_df, paste(path_to_data, dist_str, "_Var_", str_mod, ".csv", sep = ""))
  write.csv(vmr_df, paste(path_to_data, dist_str, "_VMR_", str_mod, ".csv", sep = ""))
  write.csv(entropy_df, paste(path_to_data, dist_str, "_Entropy_", str_mod, ".csv", sep = ""))
}

## Function calls to generate the simulated metrics ## 
# for Poisson
generate_mats(lambdas, sample_pois, "Poisson", "vary_ngenes", 500, n_genes, F)
generate_mats(lambdas, sample_pois, "Poisson", "vary_ncells", n_cells, 5000, T)
generate_mats_fix_dispersion(sample_pois, "Poisson", "fix_dispersion", n_cells, n_genes, 1)

# for negative binomial
generate_mats(sizes, sample_rnbinom, "Negative_Binomial", "vary_ngenes", 500, n_genes, F)
generate_mats(sizes, sample_rnbinom, "Negative_Binomial", "vary_ncells", n_cells, 5000, T)
generate_mats_fix_dispersion(sample_rnbinom, "Negative_Binomial", "fix_dispersion", n_cells, n_genes, 1)

# for beta-Poisson
generate_mats(scale_factors, generate_betapois, "Beta-Poisson", "vary_ngenes", 500, n_genes, F)
generate_mats(scale_factors, generate_betapois, "Beta-Poisson", "vary_ncells", n_cells, 5000, T)
generate_mats_fix_dispersion(generate_betapois, "Beta-Poisson", "fix_dispersion", n_cells, n_genes, 1)

# for Uniform 
generate_mats(bs, sample_unif, "Uniform", "vary_ngenes", 500, n_genes, F)
generate_mats(bs, sample_unif, "Uniform", "vary_ncells", n_cells, 5000, T)
generate_mats_fix_dispersion(sample_unif, "Uniform", "fix_dispersion", n_cells, n_genes, 2)

