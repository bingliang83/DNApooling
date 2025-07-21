# DNApooling

**DNApooling** is an R package for estimating parental contributions to DNA pools using SNP genotype data. It applies differential evolution optimization to infer family proportions from allele frequencies observed in pooled offspring samples. This is particularly useful for aquaculture and livestock breeding programs where parentage is unknown but genotype information is available.

---

## 📦 Installation

```r
# Install from GitHub
remotes::install_github("bingliang83/DNApooling")
```

---

## 🧪 Usage

### 🔹 Option 1: Run Example Analysis (with or without known parent sex)

```r
library(DNApooling)
run_analysis(use_example = TRUE)
```

```r
library(DNApooling)
run_analysis_unsexed(use_example = TRUE)
```


This loads example input files bundled in the package and outputs results to the current working directory.

---

### 🔹 Option 2: Run Analysis with Known Parent Sex

Use `run_analysis()` if your `pheno_parents.txt` file includes a `sex` column.

```r
run_analysis(
  geno_parents_file = "geno_parents.txt",
  pheno_parents_file = "pheno_parents.txt",
  pheno_off_file = "pheno_off.txt",
  af_pool_file = "af_pool.txt",
  out_dir = "output"
)
```

---

### 🔹 Option 3: Run Analysis Without Parent Sex Info

Use `run_analysis_unsexed()` if your `pheno_parents.txt` file does **not** include a sex column or the sex information is unknown.

```r
run_analysis_unsexed(
  geno_parents_file = "geno_parents.txt",
  pheno_parents_file = "pheno_parents.txt",
  pheno_off_file = "pheno_off.txt",
  af_pool_file = "af_pool.txt",
  out_dir = "output"
)
```

---

## 📁 Input File Descriptions

Before running the analysis, make sure your input files are accessible by R. You can:

- Place the files in your **current working directory** (check with `getwd()` in R), or  
- Provide **relative or absolute paths** to the input files when calling the function.

For example:

```r
geno_parents_file = "data/geno_parents.txt"        # relative path
geno_parents_file = "/home/user/files/geno_parents.txt"  # absolute path
```

### 🧬 `geno_parents.txt`
- SNP genotype matrix for parents
- **Rows**: Parent IDs (e.g., `id0023`, `id0024`)
- **Columns**: SNP IDs
- Genotypes must be coded as `0`, `1`, or `2` (allele dosage)

---

### 📋 `pheno_parents.txt`
- Metadata for parental individuals
- **Required column**:
  - `ID`: Unique parent identifier (must match `geno_parents.txt`)
- **Optional column**:
  - `sex`: Biological sex of parent (`1` = sire/male, `2` = dam/female)
    - Required only when using `run_analysis()`
    - Not needed for `run_analysis_unsexed()`
- *Note: The column must be explicitly named `sex` if included.*

---

### 🧒 `pheno_off.txt`
- Metadata for offspring used in DNA pooling
- **Required columns**:
  - `ID`: Unique offspring identifier
  - `pool`: Pooling status (`1` = included in pool, `0` = excluded)
- *Note: This file must contain a column named `pool` for correct processing.*

---

### 📈 `af_pool.txt`
- Observed allele frequencies from pooled offspring
- Derived from quantitative genotyping of the DNA pool
- **Format**:
  - Single column
  - **No header**
  - Each row represents a SNP
  - SNPs must be in the same order as in `geno_parents.txt`



---

## 📤 Output Files

All output files are written to the directory specified by the `out_dir` parameter. By default, this is the current working directory (`out_dir = "."`). You can set a custom output folder:

```r
run_analysis(..., out_dir = "output")
```

If the directory does not exist, R will return an error. You can create it before running the analysis with:

```r
dir.create("output")
```

### Output Includes:

- `ContribSolutionRepliX.txt`: Estimated family contributions (replicate X)
- `AlleleFreqSolutionRepliX.txt`: Estimated allele frequencies (replicate X)
- `est_parent_contrib_final_all.csv`: Combined parent-level contributions across replicates
- `result.txt`: Summary of DEoptim parameters and convergence information

---

## 📚 Function Documentation

### `run_analysis()`

Use when `pheno_parents.txt` includes a `sex` column.

```r
run_analysis(
  geno_parents_file,
  pheno_parents_file,
  pheno_off_file,
  af_pool_file,
  out_dir = ".",
  maxgen = 100000,
  nrep = 5,
  popsize_factor = 10,
  use_example = FALSE
)
```

---

### `run_analysis_unsexed()`

Use when `pheno_parents.txt` does **not** include sex info. Assumes any × any parent combinations.

```r
run_analysis_unsexed(
  geno_parents_file,
  pheno_parents_file,
  pheno_off_file,
  af_pool_file,
  out_dir = ".",
  maxgen = 100000,
  nrep = 5,
  popsize_factor = 10
)
```

---

## 🔍 Example Input Files

To inspect example input formats, check:

```
inst/extdata/
```

---

## 📬 Support

If you encounter any issues or have suggestions, please [open a GitHub Issue](https://github.com/bingliang83/DNApooling/issues).
