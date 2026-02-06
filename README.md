# IJDS Paper - Importance Sampling with Variable Selection

## Overview

This repository contains the code for numerical experiments in the IJDS paper on importance sampling with variable selection (IS-VS). The experiments compare our proposed IS-VS method against several baseline methods across four examples of increasing dimensionality.

## Directory Structure

```
IJDS_paper_code/
├── README.md
├── AMISE_bivariate.R          # Bandwidth selection via AMISE
├── Additive_S.R               # Additive kernel estimator
├── cv.R                       # Cross-validation utilities
├── shapley_value.R            # Shapley value computation
├── functions_ex1.R            # Example 1 specific functions (D=4)
├── functions_ex2.R            # Example 2 specific functions (D=5, correlated)
├── functions_ex3.R            # Example 3 specific functions (D=10, independent)
├── functions_ex4.R            # Example 4 specific functions (D=50)
│
├── Ex1/
│   ├── ex1_main.R             # Dispatcher script
│   └── methods/
│       ├── ex1_IS_VS.R        # Proposed method
│       ├── ex1_IS_CE.R
│       ├── ex1_IS_Pareto.R
│       ├── ex1_wamk.R
│       ├── ex1_Lasso.R
│       ├── ex1_SpAM.R
│       ├── ex1_randomForest.R
│       ├── ex1_randomForest_VS.R
│       └── ex1_GP_VS.R
│
├── Ex2/
│   ├── ex2_main.R
│   └── methods/
│       ├── ex2_IS_VS.R        # Proposed method
│       ├── ex2_IS_CE.R
│       ├── ex2_IS_Pareto.R
│       ├── ex2_wamk.R
│       ├── ex2_Lasso.R
│       ├── ex2_SpAM.R
│       ├── ex2_randomForest.R
│       ├── ex2_randomForest_VS.R
|       └── ex2_GP_VS.R
│
├── Ex3/
│   ├── ex3_main.R
│   └── methods/
│       ├── ex3_IS_VS.R        # Proposed method
│       ├── ex3_IS_CE.R
│       ├── ex3_IS_Pareto.R
│       ├── ex3_wamk.R
│       ├── ex3_SpAM.R
│       |── ex3_randomForest.R
│       ├── ex3_randomForest_VS.R
|       └── ex3_GP_VS.R
│
└── Ex4/
    ├── ex4_main.R
    └── methods/
        ├── ex4_IS_VS.R        # Proposed method
        ├── ex4_IS_CE.R
        ├── ex4_IS_Pareto.R
        ├── ex4_wamk.R
        ├── ex4_Lasso.R
        ├── ex4_SpAM.R
        ├── ex4_randomForest.R
        ├── ex4_randomForest_VS.R
        └── ex4_GP_VS.R
```

## Example Parameters

| Parameter | Example 1 | Example 2 | Example 3 | Example 4 |
|-----------|-----------|-----------|-----------|-----------|
| **Dimension (D)** | 4 | 5 | 10 | 50 |
| **Threshold (l_target)** | 18.99 | 24.2 | 18.74 | 18.42 |
| **Kernel Pairs C(D,2)** | 6 | 10 | 45 | 1225 |
| **Seed Base** | 123 | 1234 | 123 | 123 |
| **Sampling Structure** | Independent | Correlated | Independent | Independent |
| **Functions File** | functions_ex1.R | functions_ex2.R | functions_ex3.R | functions_ex4.R |

### Sampling Structures

**Independent (Ex1, Ex3, Ex4):**
```r
X ~ N(0, I_D)  # iid standard normal
```

**Correlated (Ex2):**
```r
X1, X4, X5 ~ iid N(0, 1)
X2 | X1 ~ N(X1, 1)
X3 | X1 ~ N(X1, 1)
```

## Methods Overview

| Method | Function | Description |
|--------|----------|-------------|
| **IS-VS** | `run_is_vs()` | **Proposed method.** Greedy CE ordering + Pareto-k filtering + validation-based subset selection |
| **IS-CE** | `run_is_ce()` | Greedy CE ordering only (no Pareto-k filtering or validation) |
| **IS-Pareto** | `run_is_pareto()` | Greedy CE ordering + Pareto-k filtering (takes first passing subset) |
| **WAMK-SIS** | `run_wamk_sis()` | Weighted Additive Multi-Kernel (uses all kernel pairs) |
| **Lasso** | `run_lasso()` | L1-regularized logistic regression |
| **SpAM** | `run_spam()` | Sparse Additive Models with CV-based λ selection |
| **RF-RFE** | `run_rf_rfe()` | Random Forest with Recursive Feature Elimination |
| **IS-VS-RF** | `run_is_vs_rf()` | Random Forest with greedy CE ordering + validation selection |
| **GP-VS** | `run_gp_vs()` | Gaussian Process with Pareto-k filtering + validation selection |

## Usage

### Running a Single Method

```r
# Set working directory to IJDS_paper_code
setwd("path/to/IJDS_paper_code")

# Source the main dispatcher
source("Ex1/ex1_main.R")

# Run a specific method
result <- run_ex1("IS-VS", n_reps = 25, ncores = 25)

# Access results
print(result$POE_est)
print(result$elapsed_time)
```

### Running All Methods

```r
source("Ex1/ex1_main.R")
all_results <- run_all_methods(n_reps = 25, ncores = 25)
```

### Command Line Interface

```bash
Rscript Ex1/ex1_main.R --method=IS-VS --reps=25 --cores=25
```

### Function Parameters

All method functions share the same signature:

```r
run_method_name(
  n_reps = 25,      # Number of independent replications

  ncores = 25,      # Number of parallel cores
  seed_base = 123,  # Base seed for reproducibility
  ROOT = getwd()    # Root directory for sourcing dependencies
)
```

## Key Technical Details

### Greedy Subset Search Functions

| Function | Used In | Description |
|----------|---------|-------------|
| `S_find_h_s_ce_order()` | Ex1, Ex2 | Exhaustive CE ordering for small D |
| `sampling_ce_order_greedy()` | Ex3, Ex4 | Greedy search with Pareto-k filter |
| `find_ce_order_greedy()` | IS-CE methods | CE ordering without Pareto-k |

### Validation Functions

| Function | Description |
|----------|-------------|
| `S_find_index_validation()` | Kernel-based validation for IS-VS |
| `S_find_index_validation_GP()` | GP-based validation |
| `S_find_index_validation_randomForest()` | RF-based validation |

### Key Differences by Dimension

- **Low-D (Ex1, D=4):** Enumerate all 2^D - 1 = 15 subsets
- **Medium-D (Ex2, D=5):** Enumerate all 2^D - 1 = 31 subsets with correlated structure
- **High-D (Ex3, D=10):** Greedy search with `sampling_ce_order_greedy()` for efficient exploration
- **Very High-D (Ex4, D=50):** No bandwidth optimization; use `h_new_mat = h_init_mat`

## Dependencies

### R Packages

```r
# Core
install.packages(c("doParallel", "foreach"))

# Numerical
install.packages(c("pracma", "mvtnorm", "cubature", "lhs"))

# Diagnostics
install.packages("loo")

# Machine Learning
install.packages(c("glmnet", "SAM", "randomForest", "caret"))

# Gaussian Process
install.packages(c("GauPro", "DiceKriging", "kernlab"))
```

### Shared Source Files

All methods require these files in the ROOT directory:

| File | Description |
|------|-------------|
| `AMISE_bivariate.R` | Bandwidth selection via Asymptotic MISE |
| `Additive_S.R` | Additive kernel density estimator |
| `cv.R` | Cross-validation utilities (includes `S_find_index_validation`) |
| `functions_exN.R` | Example-specific functions (sampling, `get_Y`, `mu`, etc.) |

## Output Format

Each method returns a list containing:

```r
list(
  method = "IS-VS",
  POE_est = <list of replication results>,
  n_reps = 25,
  elapsed_time = <difftime>,
  timestamp = <POSIXct>
)
```

Results are automatically saved to `ExN/results/exN_method_name.Rdata`.

## Notes

1. **Parallel Execution:** All methods use `doParallel` with `foreach`. The cluster is properly cleaned up with `on.exit(stopCluster(cl))`.

2. **Reproducibility:** Each replication uses `set.seed(seed_base + j)` for reproducibility.

3. **Error Handling:** Methods use `.errorhandling = 'remove'` in foreach to skip failed replications.

4. **Source Files in Workers:** All `source()` calls are inside the foreach loop to ensure worker processes have access to required functions.

## Citation

If you use this code, please cite:

```
[Paper citation to be added]
```

## Contact

For questions or issues, please contact the authors.
