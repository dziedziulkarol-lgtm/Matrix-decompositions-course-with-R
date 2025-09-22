# R_codes — Image Decompositions in R

This folder contains R scripts for experimenting with **matrix decompositions on images**.

1. **`loadjpg.R`**
   - Loads a `.jpg` image (e.g., your portrait or *Stańczyk* by Jan Matejko).
   - Converts it to grayscale and rescales it to an `N×N` matrix (default: `128×128`).
   - Saves the result as a CSV file (`img_small.csv`).

2. **`decomposition_face.R`**
   - Reads the CSV matrix from step 1.
   - Performs matrix decompositions (SVD, CUR, optionally NMF).
   - Compares reconstruction error, PSNR, runtime, and compression gain.
   - Saves a report (`decomposition_report.csv`) and reconstructed images.

---

## Requirements

- R (≥ 4.0)  
- Packages:

```r
install.packages(c("jpeg","MASS"))   # required
install.packages("NMF")              # optional, Nonnegative Matrix Factorization
install.packages("NNLM")             # alternative package for NMF


