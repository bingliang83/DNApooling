# DNApooling

**DNApooling** is an R package for estimating parental contributions to DNA pools using SNP genotype data. It applies **Differential Evolution (DE)** optimization to infer family proportions based on allele frequencies observed in pooled offspring samples.

This method is especially useful for **mass-spawning species**, where tracking individual mating pairs is often impractical. DNApooling enables researchers and breeders to analyze pooled DNA data efficiently, even without prior knowledge of parental sex. However, including sex information can improve speed and accuracy, particularly in large spawning populations.

## Key Features

* Estimates relative parental contributions in pooled samples.
* Compatible with datasets with or without known parental sex.
* Uses observed allele frequencies and parent genotypes as input.
* Allows optional tuning of **Differential Evolution parameters (`F` and `CR`)** for optimization.

## Applications

DNApooling is particularly beneficial in breeding programs and genetic studies involving mass-spawning aquatic species. Potential use cases include:

* **Progeny testing** without individual tagging.
* **Investigating spawning patterns** and contributions in mass spawning events.
* **Estimating genetic diversity** within pooled offspring samples.
* Supporting **broodstock management** and **selective breeding** decisions.

## Requirements

* Genotype data for all candidate parents.
* Allele frequency data from pooled offspring.
* Optional: sex information for each parent.

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

## ⚙️ Optional: Adjust Differential Evolution Parameters

DNApooling uses **Differential Evolution optimization** via the `DEoptim` package.
Advanced users may optionally modify the algorithm parameters:

* **F** – mutation factor
* **CR** – crossover probability

Default values used in the package:

```
F = 0.8
CR = 0.5
```

Users can override these values when calling the functions.

### Example

```r
run_analysis(
  geno_parents_file = "geno_parents.txt",
  pheno_parents_file = "pheno_parents.txt",
  pheno_off_file = "pheno_off.txt",
  af_pool_file = "af_pool.txt",
  out_dir = "output",
  F = 0.1,
  CR = 0.9
)
```

Similarly for unsexed analysis:

```r
run_analysis_unsexed(
  geno_parents_file = "geno_parents.txt",
  pheno_parents_file = "pheno_parents.txt",
  pheno_off_file = "pheno_off.txt",
  af_pool_file = "af_pool.txt",
  out_dir = "output",
  F = 0.1,
  CR = 0.9
)
```

If these parameters are not specified, the default values **F = 0.8** and **CR = 0.5** will be used.

---

## 📁 Input File Descriptions

Before running the analysis, make sure your input files are accessible by R. You can:

* Place the files in your **current working directory** (check with `getwd()` in R), or
* Provide **relative or absolute paths** to the input files when calling the function.

Example:

```r
geno_parents_file = "data/geno_parents.txt"
geno_parents_file = "/home/user/files/geno_parents.txt"
```

---

### 🧬 `geno_parents.txt`

* SNP genotype matrix for parents
* **Rows**: Parent IDs (e.g., `id0023`, `id0024`)
* **Columns**: SNP IDs
* Genotypes must be coded as `0`, `1`, or `2` (allele dosage)

---

### 📋 `pheno_parents.txt`

Metadata for parental individuals.

Required column:

```
ID
```

Optional column:

```
sex
```

Coding:

```
1 = sire (male)
2 = dam (female)
```

Notes:

* Required only when using `run_analysis()`
* Not needed for `run_analysis_unsexed()`
* Column must be explicitly named `sex` if included.

---

### 🧒 `pheno_off.txt`

Metadata for offspring used in DNA pooling.

Required columns:

```
ID
pool
```

Coding:

```
1 = included in pool
0 = excluded
```

---

### 📈 `af_pool.txt`

Observed allele frequencies from pooled offspring.

Format:

* Single column
* **No header**
* One row per SNP
* SNP order must match `geno_parents.txt`

---

## 📤 Output Files

All output files are written to the directory specified by the `out_dir` parameter. By default, this is the current working directory (`out_dir = "."`).

Example:

```r
run_analysis(..., out_dir = "output")
```

If the directory does not exist:

```r
dir.create("output")
```

### Output Files

```
ContribSolutionRepliX.txt
AlleleFreqSolutionRepliX.txt
est_parent_contrib_final_all.csv
result.txt
```

Descriptions:

* **ContribSolutionRepliX.txt** – estimated family contributions for replicate X
* **AlleleFreqSolutionRepliX.txt** – estimated allele frequencies for replicate X
* **est_parent_contrib_final_all.csv** – combined parent-level contributions across replicates
* **result.txt** – summary of optimization parameters and convergence information

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
  F = 0.8,
  CR = 0.5,
  use_example = FALSE
)
```

---

### `run_analysis_unsexed()`

Use when `pheno_parents.txt` does **not** include sex information.

```r
run_analysis_unsexed(
  geno_parents_file,
  pheno_parents_file,
  pheno_off_file,
  af_pool_file,
  out_dir = ".",
  maxgen = 100000,
  nrep = 5,
  popsize_factor = 10,
  F = 0.8,
  CR = 0.5,
  use_example = FALSE
)
```

---

## 🔍 Example Input Files

Example input files are included in the package:

```
inst/extdata/
```

---

## 📬 Support

If you encounter any issues or have suggestions, please open a GitHub Issue:

https://github.com/bingliang83/DNApooling/issues
