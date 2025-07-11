# DNApooling

**DNApooling** is an R package for estimating parental contributions to DNA pools using SNP genotype data. It applies differential evolution optimization to infer family proportions from allele frequencies observed in pooled offspring samples. This is particularly useful for aquaculture and livestock breeding programs where parentage is unknown but genotype information is available.

---

## 📦 Installation

```r
# Install from GitHub
remotes::install_github("bingliang83/DNApooling")
🧪 Usage
🔹 Option 1: Run Example Analysis (with known parent sex)
r
Copy
Edit
library(DNApooling)
run_analysis(use_example = TRUE)
This loads example input files bundled in the package and outputs results to the current working directory.

🔹 Option 2: Run Analysis with Known Parent Sex
Use run_analysis() if you have sex information for each parent (i.e., pheno_parents.txt includes a sex column).

r
Copy
Edit
run_analysis(
  geno_parents_file = "geno_parents.txt",
  pheno_parents_file = "pheno_parents.txt",
  geno_off_file = "geno_off.txt",
  pheno_off_file = "pheno_off.txt",
  af_pool_file = "af_pool.txt",
  out_dir = "output"
)
🔹 Option 3: Run Analysis without Parent Sex Info
Use run_analysis_unsexed() if you do not have sex information for each parent (i.e., pheno_parents.txt does not include or uses a placeholder for sex).

r
Copy
Edit
run_analysis_unsexed(
  geno_parents_file = "geno_parents.txt",
  pheno_parents_file = "pheno_parents.txt",
  geno_off_file = "geno_off.txt",
  pheno_off_file = "pheno_off.txt",
  af_pool_file = "af_pool.txt",
  out_dir = "output"
)
📁 Input File Formats
🧬 geno_parents.txt
Matrix of SNP genotypes for parents

Row names: Parent IDs (e.g., id0023, id0024)

Column names: SNP IDs

Genotypes coded as 0, 1, or 2 (allele dosage)

📋 pheno_parents.txt
Parent metadata file

Required columns:

ID

sex (1 = sire/male, 2 = dam/female) — optional if using run_analysis_unsexed()

🧬 geno_off.txt
Matrix of SNP genotypes for offspring

Row names: Offspring IDs

Same SNPs (columns) and order as in geno_parents.txt

📋 pheno_off.txt
Offspring metadata file

Required columns:

ID

pool (1 = included in DNA pool, 0 = excluded)

📈 af_pool.txt
A numeric matrix with observed allele frequencies in the pooled sample

One column, no header

Each row corresponds to a SNP (same order as genotype files)

📤 Output Files
The following files are written to the specified out_dir:

ContribSolutionRepliX.txt: Estimated family contributions for replicate X

AlleleFreqSolutionRepliX.txt: Estimated allele frequencies for replicate X

est_parent_contrib_final_all.csv: Combined parental contributions across replicates

result.txt: Summary table with DEoptim parameters and convergence info

📚 Function Documentation
run_analysis()
r
Copy
Edit
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
Use this function if you have sex column in pheno_parents.txt. Parent combinations will be formed as sire × dam.

run_analysis_unsexed()
r
Copy
Edit
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
Use this function if you do not have sex info. It assumes any × any parent combination.

🔍 Example Input Files
To inspect the example input formats, visit:

📁 inst/extdata/

📬 Support
If you encounter issues or have suggestions, please open a GitHub Issue.