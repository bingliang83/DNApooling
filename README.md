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

### 🔹 Option 1: Run Example Analysis (with known parent sex)

```r
library(DNApooling)
run_analysis(use_example = TRUE)
```

This loads example input files bundled in the package and outputs results to the current working directory.

---

### 🔹 Option 2: Run Analysis with Known Parent Sex

Use `run_analysis()` if your `pheno_parents.txt` file includes a `sex` column.

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

---

### 🔹 Option 3: Run Analysis Without Parent Sex Info

Use `run_analysis_unsexed()` if your `pheno_parents.txt` file does **not** include a sex column or the sex information is unknown.

```r
run_analysis_unsexed(
  geno_parents_file = "geno_parents.txt",
  pheno_parents_file = "pheno_parents.txt",
  geno_off_file = "geno_off.txt",
  pheno_off_file = "pheno_off.txt",
  af_pool_file = "af_pool.txt",
  out_dir = "output"
)
```

---

## 📁 Input File Formats

### 🧬 `geno_parents.txt`
- SNP genotype matrix for parents
- **Row names**: Parent IDs (e.g., `id0023`, `id0024`)
- **Column names**: SNP IDs
- Genotypes coded as `0`, `1`, or `2` (allele dosage)

### 📋 `pheno_parents.txt`
- Metadata for parents
- Required column:
  - `ID`
  - `sex` (1 = sire/male, 2 = dam/female) — *optional if using `run_analysis_unsexed()`*

### 🧬 `geno_off.txt`
- SNP genotype matrix for offspring
- **Row names**: Offspring IDs
- **Columns**: SNPs (same order as in `geno_parents.txt`)

### 📋 `pheno_off.txt`
- Metadata for offspring
- Required columns:
  - `ID`
  - `pool` (1 = included in DNA pool, 0 = excluded)

### 📈 `af_pool.txt`
- A matrix of observed allele frequencies from the pool
- **One column**, **no header**
- Each row corresponds to a SNP (same order as genotype files)

---

## 📤 Output Files

Written to the specified `out_dir`:

- `ContribSolutionRepliX.txt`: Estimated family contributions (replicate X)
- `AlleleFreqSolutionRepliX.txt`: Estimated allele frequencies (replicate X)
- `est_parent_contrib_final_all.csv`: Combined parent-level contributions
- `result.txt`: DEoptim summary including parameters and convergence info

---

## 📚 Function Documentation

### `run_analysis()`

Use when `pheno_parents.txt` includes a `sex` column.

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

---

### `run_analysis_unsexed()`

Use when `pheno_parents.txt` does **not** include sex info. Assumes any × any parent combinations.

```r
run_analysis_unsexed(
  geno_parents_file,
  pheno_parents_file,
  geno_off_file,
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
