# DNApooling

**DNApooling** is an R package for estimating parental contributions to DNA pools using SNP genotype data. It applies differential evolution optimization to infer family proportions from allele frequencies observed in pooled offspring samples. This is particularly useful for aquaculture and livestock breeding programs.

---

## 📦 Installation

```r
# From GitHub
remotes::install_github("bingliang83/DNApooling")
```

---

## 🧪 Usage

### 🔹 Run Example Analysis
The package includes built-in example data.

```r
library(DNApooling)
run_analysis(use_example = TRUE)
```

This will read example files from the package and output results to the current directory.

---

### 🔹 Run Your Own Analysis
Provide paths to your own input files:

```r
run_analysis(
  geno_parents_file = "geno_parents.txt",
  pheno_parents_file = "pheno_parents.txt",
  geno_off_file = "geno_off.txt",
  pheno_off_file = "pheno_off.txt",
  af_pool_file = "af_pool.txt",
  out_dir = "output"
)
```

Ensure your input files follow the format shown in the example files in [`inst/extdata/`](https://github.com/bingliang83/DNApooling/tree/main/inst/extdata).

---

## 📁 Input File Formats

### 🧬 `geno_parents.txt`
- **Matrix of parental SNP genotypes**
- **Row names:** Parent IDs (e.g., `id0023`, `id0024`)
- **Column names:** SNP IDs
- Values: `0`, `1`, or `2` (allele dosage)

### 📋 `pheno_parents.txt`
- **Metadata for parents**
- Required columns:
  - `ID` (e.g., `id0023`, `id0024`)
  - `sex` (1 = male/sire, 2 = female/dam)

### 🧬 `geno_off.txt`
- **Matrix of offspring SNP genotypes**
- **Row names:** Offspring IDs
- **Column names:** SNP IDs (must match `geno_parents.txt`)

### 📋 `pheno_off.txt`
- **Metadata for offspring**
- Required columns:
  - `ID`
  - `pool` (1 = included in pool, 0 = not included)

### 📈 `af_pool.txt`
- **Numeric matrix of allele frequencies in the pool**
- One column (for one pool), **no header**
- Each row corresponds to a SNP (same order as in genotype files)

---

## 📤 Output Files
- `ContribSolutionRepliX.txt`: Estimated family contributions
- `AlleleFreqSolutionRepliX.txt`: Estimated allele frequencies
- `est_parent_contrib_final_all.csv`: Combined parent-level contributions
- `result.txt`: Summary of replicates

---

## 📚 Documentation

### Function: `run_analysis()`

```r
run_analysis(
  geno_parents_file,
  pheno_parents_file,
  geno_off_file,
  pheno_off_file,
  af_pool_file,
  out_dir = ".",
  maxgen = 100000,
  nrep = 5,
  popsize_factor = 10,
  use_example = FALSE
)
```

#### @details
Input files must follow the formats:
- **geno_parents.txt**: SNP genotype matrix with parent IDs as rownames, SNP IDs as column names.
- **pheno_parents.txt**: Data frame with columns `ID` (e.g. `id0023`) and `sex` (1 = sire, 2 = dam).
- **geno_off.txt**: Genotype matrix of offspring, with same SNPs and order as `geno_parents.txt`, rownames = offspring IDs.
- **pheno_off.txt**: Metadata with columns `ID` and `pool` (1 = in pool, 0 = not).
- **af_pool.txt**: Matrix of observed pooled allele frequencies (no header, one column).

#### @examples
```r
# Example with package data
run_analysis(use_example = TRUE)

# Example with user-provided input
run_analysis(
  geno_parents_file = "geno_parents.txt",
  pheno_parents_file = "pheno_parents.txt",
  geno_off_file = "geno_off.txt",
  pheno_off_file = "pheno_off.txt",
  af_pool_file = "af_pool.txt",
  out_dir = "output"
)
```

---

Let us know via GitHub Issues if you encounter bugs or have suggestions!
