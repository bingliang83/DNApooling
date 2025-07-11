#' Run Family Contribution Analysis
#'
#' This function reads genomic and phenotypic input files, runs DEoptim optimization,
#' and saves results including allele frequencies and parent contributions.
#'
#' @param geno_parents_file Path to geno_parents.txt
#' @param pheno_parents_file Path to pheno_parents.txt
#' @param geno_off_file Path to geno_off.txt
#' @param pheno_off_file Path to pheno_off.txt
#' @param af_pool_file Path to af_pool.txt
#' @param out_dir Directory to save result files (default: current directory)
#' @param maxgen Maximum number of generations (default: 100000)
#' @param nrep Number of replicates (default: 5)
#' @param popsize_factor Multiplier to determine population size (default: 10)
#' @param use_example If TRUE, use example data from inst/extdata/ (default: FALSE)
#'
#' @return Saves result files to out_dir
#' @export
run_analysis <- function(geno_parents_file = NULL,
                         pheno_parents_file = NULL,
                         geno_off_file = NULL,
                         pheno_off_file = NULL,
                         af_pool_file = NULL,
                         out_dir = ".",
                         maxgen = 100000,
                         nrep = 5,
                         popsize_factor = 10,
                         use_example = FALSE) {

  library(DEoptim)
  library(ggplot2)
  library(gaston)
  library(dplyr)

  # Use packaged example data if use_example = TRUE
  if (use_example) {
    geno_parents_file <- system.file("extdata", "geno_parents.txt", package = "DNApooling")
    pheno_parents_file <- system.file("extdata", "pheno_parents.txt", package = "DNApooling")
    geno_off_file      <- system.file("extdata", "geno_off.txt", package = "DNApooling")
    pheno_off_file     <- system.file("extdata", "pheno_off.txt", package = "DNApooling")
    af_pool_file       <- system.file("extdata", "af_pool.txt", package = "DNApooling")
  }

  # Validate paths
  if (any(sapply(list(geno_parents_file, pheno_parents_file, geno_off_file, pheno_off_file, af_pool_file), is.null))) {
    stop("Missing one or more input file paths. Set use_example = TRUE to use built-in example files.")
  }

  # Read input files
  geno_parents <- read.table(geno_parents_file, stringsAsFactors = FALSE)
  pheno_parents <- read.table(pheno_parents_file, header = TRUE)
  geno_off <- read.table(geno_off_file, stringsAsFactors = FALSE)
  pheno_off <- read.table(pheno_off_file, header = TRUE, stringsAsFactors = FALSE)
  Mpool <- as.matrix(read.table(af_pool_file, header = FALSE))

  nsnps <- ncol(geno_parents)
  np <- 1
  sample_perpool <- sum(pheno_off$pool == 1)

  #create all possible combinations among the parents and use all of them as possible families
  ids <- pheno_parents$id
  combinations <- combn(ids, 2)
  combinations_df <- as.data.frame(t(combinations))
  families_all <- combinations_df
  families_all <- families_all %>%
    rename(siresim=V1,damsim=V2)
  nfam <- length(families_all$siresim)
  M <- matrix(c(0),ncol=nfam,nrow=nsnps)

  for(i in 1:nfam)
  {
    gen1 <- geno_parents[which(rownames(geno_parents)==families_all$siresim[i]),]
    gen2 <- geno_parents[which(rownames(geno_parents)==families_all$damsim[i]),]

    d <- 1

    while(d <= nsnps)
    {

      M[d,i] <- sum(gen1[d],gen2[d])/4 #e.g. for family 1 snp 1, the allele frequency is
      #(V1[1] + V2[1])/ 4 = 0.25

      d <- d + 1
    }
    d <- 1

  }

  F <- 0.8
  CR <- 0.5
  popsize <- popsize_factor * nfam

  of <- function(x) {
    v <- matrix(x, ncol = np, nrow = nfam, byrow = TRUE)
    z <- M %*% v
    sum(abs(z - Mpool)) + (sum(v) - 1)^2 * 1e3
  }

  est_parent_contrib_final <- list()
  siresim_list <- list()
  damsim_list <- list()
  res_list <- list()

  for(r in 1:nrep) {
    outDEoptim <- DEoptim(of, lower = rep(0, nfam*np), upper = rep(1, nfam*np),
                          DEoptim.control(strategy = 1, VTR = 0.5, itermax = maxgen, NP = popsize,
                                          trace = TRUE, F = F, CR = CR))

    best_v <- matrix(outDEoptim$optim$bestmem, ncol = np, nrow = nfam, byrow = TRUE)
    colnames(best_v) <- "contrib"  # assign column name

    Yest <- M %*% best_v

    write.table(best_v, file.path(out_dir, paste0("ContribSolutionRepli", r, ".txt")),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(Yest, file.path(out_dir, paste0("AlleleFreqSolutionRepli", r, ".txt")),
                quote = FALSE, row.names = FALSE, col.names = FALSE)

    est_contrib <- cbind(families_all, best_v)

    siresim_list[[r]] <- est_contrib %>%
      group_by(siresim) %>%
      summarise(contrib_percent = sum(contrib), .groups = "drop") %>%
      rename(parent = siresim) %>%
      mutate(parent_type = 'sire', replicate = r)

    damsim_list[[r]] <- est_contrib %>%
      group_by(damsim) %>%
      summarise(contrib_percent = sum(contrib), .groups = "drop") %>%
      rename(parent = damsim) %>%
      mutate(parent_type = 'dam', replicate = r)

    est_parent_contrib_final[[r]] <- rbind(siresim_list[[r]], damsim_list[[r]]) %>% group_by(parent) %>%
      summarise(contrib_percent = sum(contrib_percent), .groups = "drop") %>%
      arrange(parent)

    res_list[[r]] <- data.frame(replicate = r, pop = popsize, gen = maxgen,
                                lastiter = min(which(outDEoptim$member$bestvalit == outDEoptim$optim$bestval)))
  }

  combined_contrib <- do.call(rbind, est_parent_contrib_final)
  write.csv(combined_contrib, file.path(out_dir, "est_parent_contrib_final_all.csv"), row.names = FALSE)

  combined_res <- do.call(rbind, res_list)
  write.table(combined_res, file.path(out_dir, "result.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE)
}
