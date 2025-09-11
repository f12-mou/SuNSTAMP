# SuNSTAMP Analysis Repository

## Repository Structure

This repository contains replicate analyses organized in a folder-wise structure. Each replicate follows the pattern: `SuNSTAMP/SuNSTAMP-AD/replicate/[month]_[replicate_number]/`

### Folder Organization

Each replicate directory contains four main folders:
- `data/` - Raw and processed data files
- `src/` - Source code and analysis scripts  
- `result/` - Analysis outputs and results
- `figure/` - Generated plots and visualizations

### Data Access

**Important**: Due to file size limitations, the repository structure is split across two locations:

1. **This GitHub Repository**: Contains complete folder structure for `src/`, `result/`, and `figure/` directories for all replicates
2. **Google Drive**: Contains only the `data/` folders for each replicate

**To use this repository**: Download both the GitHub repository and the corresponding `data/` folders from the Google Drive link below, then merge them to recreate the complete directory structure.

📁 **Data folders available at**: https://drive.google.com/drive/folders/1jXY8URujCCVqC4V-0jcnMtVKxvr5HCgh?usp=sharing

### Example Structure
```
SunSTAMP/SuNSTAMP-AD/replicate/8_1/
├── data/          # ← Download from Google Drive
├── src/           # ← Available in this repo
├── result/        # ← Available in this repo
└── figure/        # ← Available in this repo
```

### Setup Instructions

1. Clone this repository:
   ```bash
   git clone https://github.com/f12-mou/SuNSTAMP.git
   ```

2. Download the data folders from the Google Drive link above

3. Extract and place the `data/` folders into their corresponding replicate directories

4. Verify the structure matches the example above

## Pipeline Steps  

The analysis pipeline is organized as a sequence of scripts inside each replicate’s `src/` folder. Below, we describe each script, its inputs/outputs, and its role in the workflow.  

---

<details>
  <summary> Step-1:  extraction.py</summary>

**Purpose:**  
Extracts and preprocesses a subset of the gene expression matrix corresponding to a given replicate. It selects spots (cells) by index range, transposes the matrix, and saves the result as a new CSV file for downstream analysis. The gene expression csv file should be placed in the `data` folder of the corresponding replicate.

**Inputs:**  
- `expression_matrix_normalized.csv` – Full normalized gene expression matrix (rows = genes, columns = spot IDs).  
- Parameters inside the script: `start_index_spot` and `end_index_spot` define the range of spot IDs to extract.  

**Outputs:**  
- `gene_exp_13_1.csv` (example for replicate `13_1`) – Transposed CSV file where rows = spots, columns = genes.  

**Usage:**  

You have two options:

1. **Run extraction yourself**  
   - Download the main dataset (`expression_matrix_normalized.csv`) from:  
     [Main dataset link](https://zenodo.org/records/7332091) [1]
   - Place it in the same directory as extraction.py.  
   - Run the script:  
     ```bash
     python extraction.py
     ```  

2. **Skip this step**  
   - Directly use the processed gene expression files provided in the Google Drive folder:  
     [Processed data (Google Drive)](https://drive.google.com/drive/folders/1jXY8URujCCVqC4V-0jcnMtVKxvr5HCgh?usp=sharing)  

For MI, there is a similar file named process_adata.py. If you want to run this file, you need to download the original dataset from [2]. Otherwise, you may use processed datasets from the Google Drive. 

</details>

---

<details>
  <summary>Step-2: spatial_decay_model.py (for Alzheimer's analysis)</summary>

**Purpose:**  
Calculates a spatially smoothened supervision value for each cell based on nearby plaques, normalizes the values and initial intensity, and generates visualizations. This step is part of the supervision pipeline for SuNSTAMP in AD analysis. For MI, we directly used Pathological scores from SPaSE [3] pipeline. The python file also creates two key plots:  

1. Colormap of normalized attraction values for all cells.  
2. Overlay of cells and plaques with plaque radii.  

**Inputs:**  
- `spatial_13_1.csv` – Coordinates and metadata for cells (from the `data/` folder).  
- `plaque_13_1.csv` – Coordinates and radii of plaques (from the `data/` folder).  

**Outputs:**  
- `supervision_13_1.csv` – CSV file containing:  
  - Cell IDs  
  - X/Y coordinates  
  - Supervision values  
  - Normalized supervision values  
- Figures saved in `figure/`:  
  - `normalized_attraction_values_13_1.png`  
  - `cells_and_plaques_with_radii_13_1.png`  

**Usage:**  

1. Ensure the required `data/` files for the replicate are present (`spatial_13_1.csv` and `plaque_13_1.csv`).  
2. Run the script from the replicate `src/` folder:  

```bash
python spatial_decay_model.py
```

</details>

---

<details>
  <summary>Step-3: run_nmf.sh</summary>

**Purpose:**  
Performs **k-fold nested cross-validation** for SuNSTAMP using a wide range of hyperparameters (ranks and lambda values). This step is computationally intensive and was run on an HPC cluster.  

**Inputs:**  
- `nmf_final_2.py` – Python script performing the NMF decomposition with cross-validation.  
- Preprocessed data files from the replicate `data/` folder (e.g., `gene_exp_13_1.csv`).  
- Hyperparameter configuration inside `nmf_final_2.py`.  

**Outputs:**  
- Cross-validation results for each rank and lambda combination (saved in the `result/` folder).  
- Due to the heavy computation (~23 hours runtime on 16 CPUs, 128GB RAM), **these results are already provided in the `result/` folder** for readers/users.  

**Usage:**  

- If you have access to an HPC cluster with SLURM, you can submit the job using:  

```bash
sbatch run_nmf.sh
```

**Notes:**

- All scripts assume the folder structure is maintained (data, src, result, figure).
- Users can skip heavy computations if precomputed results are available in the `result` folder.

</details>

---

<details>
  <summary>Step-4. errors_grouping.py</summary>

**Purpose:**  
Aggregates **loss values** from multiple cross-validation result files (different ranks and lambdas) into a single summary CSV.  

**Inputs:**  
- Folder containing result files: `../result/mou_errors_{months}_{replicate}/`  
  - Files are named as:  
    ```
    mou_errors_0_{months}_{replicate}.csv  
    mou_errors_1_{months}_{replicate}.csv  
    ...  
    mou_errors_20_{months}_{replicate}.csv
    ```
- Each file contains columns: `file_rank`, `lambda`, `loss`.  

**Outputs:**  
- A combined CSV file with summed losses across folds:  `../result/summed_loss_by_rank_lambda_{months}_{replicate}.csv`
  
After aggregating the cross-validation losses across all ranks and lambda values (Step 4),  
the next step is to **statistically evaluate whether consecutive ranks differ significantly**.  
We apply the **Mann-Whitney U Test** to compare the distribution of losses between rank *r* and rank *r+1*.  

  
</details>

---

<details>
  <summary>Step-5: rank_test.py</summary>

**Purpose:**  
Performs **pairwise statistical testing** across consecutive NMF ranks to check whether increasing the rank leads to a significant improvement in reconstruction loss.  

**Inputs:**  
- Aggregated loss file from Step 4: `../result/summed_loss_by_rank_lambda_{months}_{replicate}.csv`
  

**Outputs:**  
- CSV file with Mann-Whitney U test results: `../result/mann_whitney_results_summed_{months}_{replicate}.csv`
- Plot of `-log10(p-value)` vs rank, with a threshold line at **p=0.05**: `../figure/mann_whitney_log_10.png`
   

With the rank fixed, the next step is to explore how different lambda values (the supervision regularization parameter) affect the model. This allows us to tune the strength of supervision applied in SuNSTAMP and evaluate its stability across a wide range of configurations.


</details>

---

<details>
  <summary>Step-6: automate_W_matrix.sh</summary>

**Purpose:**  
Performs **supervised NMF runs** for the chosen rank while scanning across a list of lambda values. Each run produces the corresponding **W and H matrices**, which form the core of the SuNSTAMP interpretation.  

**Inputs:**  
- `new_supervised_final_run.py` – Python script implementing supervised NMF for a given `(rank, lambda)` pair.  
- data folder with Expression matrix (`gene_exp_[month]_[replicate].csv`) and supervision values (`supervision_[month]_[replicate].csv`).  

**Outputs:**  
- `result/W_Matrices/` → W matrices for each `(rank, lambda)` configuration.  
- `result/H_Matrices/` → First row of H matrix for each `(rank, lambda)` configuration.  
  
</details>

---


<details>
  <summary>Step-7: Coefficients.py </summary>

**Purpose:**  
In the previous step, we stored the **W matrices** and the first row of the **H matrices** across different configurations.  
Now, we want to assess how well the first factor (`Supervised Factor`) from the W matrix reflects the supervision signal.  
  
**Inputs:**  
- Stored W matrices from the previous step (`/result/W_matrices_{months}_{replicate}`).
- Supervision values' file (`data/supervision_{months}_{replicate}.csv`).

**Outputs:**  
- `result/Supervised_Factor_correlation_results_{months}_{replicate}.csv`
  
</details>

---

<details>
  <summary>Step-8: Predictive_Scoring.py </summary>


**Purpose:**  
Apply gene weights learned from a **source replicate** to a **target replicate** and evaluate their predictive power for **plaque proximity**. This quantifies the cross-replicate generalizability of SuNSTAMP factors.

### Inputs
- **Command-line arguments**
  - `top_n_genes`: Number of top-ranked genes from H-matrix to use.  
  - `output_file`: CSV file path to append results (Top-N genes vs. AUC).

**Inputs**
  - `gene_exp_{months}_{replicate}.csv` – Target replicate gene expression matrix (cells × genes).  
- `supervision_{months}_{replicate}.csv` – Target replicate spatial coordinates (`X-scaled`, `Y-scaled`).
-  `output_H_matrix_first_row_rank_{rank}_{lambda}_{months}_{genes_from_replicate}.csv` – Source replicate H-matrix scores (gene weights).  
  - `plaque_{months}_{replicate}.csv` – Plaque positions in the tissue (`m.cx`, `m.cy`).

### Outputs
- **Result CSV**: Rows with AUC-ROC scores  


This step demonstrates **SuNSTAMP’s cross-replicate generalization**:  
- Gene weights learned from a **source replicate** predict disease-relevant patterns (plaque proximity) in a **target replicate**.  
- High **AUC-ROC** indicates that SuNSTAMP captures **robust, biologically meaningful gene modules**.

</details>

---

<details>
  <summary>Step-9: predictive_score_automation.sh </summary>

  **Purpose:**  
Evaluate how the predictive performance of SuNSTAMP-derived gene weights varies as the number of top-ranked genes increases (from 1 to 100). This identifies the required number of genes needed for maximal plaque proximity prediction in the target replicate.

</details>

---



[1]: Zeng, H., Huang, J., Zhou, H., Meilandt, W.J., Dejanovic, B., Zhou, Y., Bohlen, C.J., Lee, S.-H., Ren, J., Liu, A., et al. *Integrative in situ mapping of single-cell transcriptional states and tissue histopathology in a mouse model of Alzheimer’s disease.* Nature Neuroscience, 26(3):430–446, 2023.


[2]: Calcagno, M., Taghdiri, N., Ninh, V.K., Mesfin, J.M., Toomu, A., Sehgal, R., Lee, J., Liang, Y., Duran, J.M., Adler, E., et al. *Single-cell and spatial transcriptomics of the infarcted heart define the dynamic onset of the border zone in response to mechanical destabilization.* Nature Cardiovascular Research, 1(11):1039–1055, 2022.


[3] Rahman, M.N., Abrar, M.A., Shaw, V.R., Martin, J.F., Rahman, M.S., Samee, M.A.H. *SPASE: Spatially resolved pathology scores using optimal transport on spatial transcriptomics data.* Cell Systems, 2025.

















