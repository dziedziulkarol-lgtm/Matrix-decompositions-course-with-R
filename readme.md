# Matrix Decompositions with R

Lecture notes and R codes for a 15-hour course on **matrix decompositions**  
taught at Gdańsk University of Technology (PG).  
The materials combine **linear algebra**, **matrix analysis**, and **hands-on R programming**.

---

## 📚 Structure

The repository is organized into thematic blocks:

### 1. Introduction and Linear Algebra Recap (`01_intro/`)
- Row/column picture of linear equations  
- Least squares and pseudoinverse (Moore–Penrose)  
- Motivation for matrix decompositions

### 2. Classical Decompositions (`02_classical/`)
- LU decomposition (Banachiewicz, Crout, Doolittle)  
- Cholesky factorization  
- Vandermonde matrices and polynomial interpolation

### 3. Low-Rank Approximations (`03_lowrank/`)
- Singular Value Decomposition (SVD)  
- CUR decomposition (column-row sampling)  
- Nonnegative Matrix Factorization (NMF)  
- Positive semidefinite approximations, Frobenius vs operator norm

### 4. Symplectic and Advanced Topics (`04_symplectic/`)
- Symplectic matrices and canonical forms  
- Williamson’s theorem  
- Infinite-dimensional symplectic geometry (brief introduction)  
- Gromov’s non-squeezing phenomenon

### 5. Open Problems (`05_open_problems/`)
- Hadamard matrices  
- Nonnegative Inverse Eigenvalue Problem (NIEP)  
- Column Subset Selection Problem (CSSP)  
- Connections to functional analysis and combinatorics

### 6. R Codes (`R_codes/`)
- R scripts for computing decompositions  
- Image approximation with SVD, NMF, CUR  
- Numerical experiments and visualization

---

## 🚀 Usage

Clone the repository:
```bash
git clone https://github.com/dziedziulkarol-lgtm/Matrix-decompositions-with-R.git
