#' =============================================================================
#' Ex2 Main Runner - IJDS Paper
#' =============================================================================
#' 
#' Dispatcher script to run different IS methods for Example 2.
#' 
#' Usage:
#'   1. Interactive: Set `method` below and source this file
#'   2. Command line: Rscript ex2_main.R --method=Lasso
#' =============================================================================

# Configuration
ROOT <- "E:/IJDS_paper_code"  # Assumes running from IJDS_paper_code directory
RESULTS_DIR <- file.path(ROOT, "Ex2", "results")

# Create results directory if needed
if (!dir.exists(RESULTS_DIR)) {
  dir.create(RESULTS_DIR, recursive = TRUE)
}

# Available methods and their source files
METHOD_FILES <- list(
  "Lasso"      = "ex2_Lasso.R",
  "IS-VS"      = "ex2_IS_VS.R",
  "IS-CE"      = "ex2_IS_CE.R",
  "IS-Pareto"  = "ex2_IS_Pareto.R",
  "WAMK-SIS"   = "ex2_wamk.R",
  "SpAM"       = "ex2_SpAM.R",
  "RF-RFE"     = "ex2_randomForest.R",
  "IS-VS-GP"   = "ex2_GP.R",
  "IS-VS-RF"   = "ex2_randomForest_VS.R"
)

# Function name convention: run_<lowercase method name with underscores>
METHOD_FUNCTIONS <- list(
  "Lasso"      = "run_lasso",
  "IS-VS"      = "run_is_vs",
  "IS-CE"      = "run_is_ce",
  "IS-Pareto"  = "run_is_pareto",
  "WAMK-SIS"   = "run_wamk_sis",
  "SpAM"       = "run_spam",
  "RF-RFE"     = "run_rf_rfe",
  "IS-VS-GP"   = "run_is_vs_gp",
  "IS-VS-RF"   = "run_is_vs_rf"
)

#' Run a single method
#' 
#' param method Character string specifying which method to run
#' param n_reps Number of replications
#' param ncores Number of parallel cores
#' param save_results If TRUE, saves results to file
#' param ... Additional arguments passed to the method function
run_ex2 <- function(method, 
                    n_reps = 25, 
                    ncores = 25, 
                    save_results = TRUE,
                    ...) {
  
  if (!method %in% names(METHOD_FILES)) {
    stop("Unknown method: ", method, "\n",
         "Available methods: ", paste(names(METHOD_FILES), collapse = ", "))
  }
  
  message("=== Running method: ", method, " ===")
  message("Replications: ", n_reps, ", Cores: ", ncores)
  
  # Source the method file
  method_file <- file.path(ROOT, "Ex2", "methods", METHOD_FILES[[method]])
  if (!file.exists(method_file)) {
    stop("Method file not found: ", method_file)
  }
  source(method_file)
  
  # Get the run function
  run_fn <- get(METHOD_FUNCTIONS[[method]])
  
  # Time the run
  start_time <- Sys.time()
  
  # Execute
  results <- run_fn(n_reps = n_reps, 
                    ncores = ncores, 
                    ROOT = ROOT,
                    ...)
  
  end_time <- Sys.time()
  elapsed <- end_time - start_time
  message("Elapsed time: ", format(elapsed))
  
  # Package results
  output <- list(
    method = method,
    POE_est = results,
    n_reps = n_reps,
    elapsed_time = elapsed,
    timestamp = Sys.time()
  )
  
  # Save if requested
  if (save_results) {
    output_file <- file.path(RESULTS_DIR, sprintf("ex2_%s.Rdata", 
                                                   gsub("-", "_", method)))
    save(output, file = output_file)
    message("Results saved to: ", output_file)
  }
  
  return(output)
}

#' Run all methods
#' 
#' param methods Vector of method names (default: all methods)
#' param ... Arguments passed to run_ex2
run_all_methods <- function(methods = names(METHOD_FILES), ...) {
  
  all_results <- list()
  
  for (method in methods) {
    tryCatch({
      all_results[[method]] <- run_ex2(method, ...)
    }, error = function(e) {
      warning("Method ", method, " failed: ", e$message)
      all_results[[method]] <<- list(method = method, error = e$message)
    })
  }
  
  return(all_results)
}

# =============================================================================
# Command-line interface
# =============================================================================
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Parse arguments
  method <- "Lasso"  # defaultresul
  n_reps <- 25
  ncores <- 25
  
  for (arg in args) {
    if (grepl("^--method=", arg)) {
      method <- sub("^--method=", "", arg)
    } else if (grepl("^--reps=", arg)) {
      n_reps <- as.integer(sub("^--reps=", "", arg))
    } else if (grepl("^--cores=", arg)) {
      ncores <- as.integer(sub("^--cores=", "", arg))
    }
  }
  
  run_ex2(method, n_reps = n_reps, ncores = ncores)
}

# =============================================================================
# Interactive usage example (uncomment to run)
# =============================================================================
# result <- run_ex2("Lasso", n_reps = 25, ncores = 25) # Done
# result <- run_ex2("IS-CE", n_reps = 25, ncores = 25) # Done
# result <- run_ex2("IS-Pareto", n_reps = 25, ncores = 25) # Done
# result <- run_ex2("RF-RFE", n_reps = 25, ncores = 25) # Done
# result <- run_ex2("WAMK-SIS", n_reps = 25, ncores = 25) # Done
# result <- run_ex2("SpAM", n_reps = 25, ncores = 25) # Done
result <- run_ex2("IS-VS", n_reps = 25, ncores = 25) #Done
print(result$POE_est)
