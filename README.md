# SuNSTAMP

SuNSTAMP is a supervised non-negative matrix factorization in which one factor is aligned with a supplied spatial supervision signal (e.g., a pathology score). This repository contains the code, saved results and figure scripts accompanying the manuscript *Supervised Factorization to Associate Spatial Transcriptomics with Complementary Molecular Readouts*.

SuNSTAMP has been applied to two datasets:

1. **Alzheimer's Disease (AD) dataset**: spatial transcriptomics data from mouse brain tissue (8- and 13-month-old mice, two replicates each).
2. **Myocardial Infarction (MI) dataset**: spatial transcriptomics data from infarcted mouse heart tissue (3 and 7 days post-MI, two replicates each).

**Paper release:** the version of this repository that accompanies the manuscript is the GitHub release tagged `v1.0.0-paper` (see the *Releases* page).

---

## Quick start

1. **Clone the repository**
   ```bash
   git clone https://github.com/f12-mou/SuNSTAMP.git
   cd SuNSTAMP
   ```
2. **Create a Python environment.** Use Python **3.10, 3.11 or 3.12**; the pinned `scipy==1.14.1` requires Python ≥ 3.10.
   ```bash
   python -m venv .venv
   source .venv/bin/activate          # Windows: .venv\Scripts\activate
   pip install -r requirements.txt
   ```
   See [Software environment](#software-environment) for details.
3. **Download the `data/` folders** from Google Drive (see [Data access](#data-access)) and place each one inside its replicate folder, e.g. `SuNSTAMP-AD/replicate_13_1/data/`.
4. **Run scripts from inside the replicate's `src/` folder.** All pipeline scripts read and write through relative paths (`../data`, `../result`, `../figure`):
   ```bash
   cd SuNSTAMP-AD/replicate_13_1/src
   python spatial_decay_model.py
   ```
5. **To reproduce figures without retraining**, see [`FIGURE_REPRODUCTION.md`](FIGURE_REPRODUCTION.md). The final W and H matrices for every analysed sample are already saved in the corresponding `result/` folder.

---

## Repository Structure

This repository contains replicate analyses organized in a folder-wise structure. Each replicate follows the pattern: `SuNSTAMP-AD[or SuNSTAMP-MI]/replicate_[time-point]_[replicate_number]/`

| Dataset | Replicate folders |
|---|---|
| AD (time point in months) | `SuNSTAMP-AD/replicate_8_1`, `replicate_8_2`, `replicate_13_1`, `replicate_13_2` |
| MI (time point in days post-MI) | `SuNSTAMP-MI/replicate_3_1`, `replicate_3_3`, `replicate_7_2`, `replicate_7_3` |

Other top-level folders:
- `Figure_Generation/`: one folder per figure panel (Figures 1–4 and Supplementary Figures S1–S4). See [`FIGURE_REPRODUCTION.md`](FIGURE_REPRODUCTION.md).
- `SuNSTAMP-AD/additional_analysis/`, `SuNSTAMP-AD/plot_table.py`, `SuNSTAMP-MI/correlation-MI_3.py`: cross-replicate summary scripts (see [Additional scripts](#additional-scripts)).

In the script variables, the time point is stored in a variable named `months` for both datasets; for MI it denotes days post-MI (3 or 7).

### Folder Organization

Each replicate directory contains four main folders:
- `data/`: raw and processed data files (downloaded separately; see below)
- `src/`: source code and analysis scripts
- `result/`: analysis outputs and results
- `figure/`: generated plots and visualizations

### Data Access

**Important**: Due to file size limitations, the repository structure is split across two locations:

1. **This GitHub Repository**: contains the complete folder structure for the `src/`, `result/` and `figure/` directories of all replicates.
2. **Google Drive**: contains only the `data/` folders for each replicate.

**To use this repository**: download both the GitHub repository and the corresponding `data/` folders from the Google Drive link below, then merge them to recreate the complete directory structure.

📁 **Data folders available at**: https://drive.google.com/drive/folders/1jXY8URujCCVqC4V-0jcnMtVKxvr5HCgh?usp=sharing

#### Files expected in each `data/` folder

These are the file names read by the scripts. `{t}` is the time point and `{r}` the replicate number, e.g. `gene_exp_13_1.csv`.

| Dataset | File | Used by |
|---|---|---|
| AD | `gene_exp_{t}_{r}.csv` (spots × genes) | `nmf_final_2.py`, `new_supervised_final_run.py`, `Predictive_Scoring.py` |
| AD | `spatial_{t}_{r}.csv` (cell coordinates and metadata) | `spatial_decay_model.py` |
| AD | `plaque_{t}_{r}.csv` (plaque coordinates and radii) | `spatial_decay_model.py`, `Predictive_Scoring.py` |
| AD | `supervision_{t}_{r}.csv` (output of `spatial_decay_model.py`) | `nmf_final_2.py`, `new_supervised_final_run.py`, `Coefficients.py`, `Predictive_Scoring.py` |
| MI | `adata_{t}_{r}.h5ad` (AnnData object) | `process_adata.py` |
| MI | `gene_exp_{t}_{r}.csv`, `spatial_{t}_{r}.csv`, `supervision_{t}_{r}.csv` (outputs of `process_adata.py`) | `nmf_final_2.py`, `new_supervised_final_run.py`, `Coefficients.py` |

The raw AD expression matrix used by `extraction.py` (`expression_matrix_normalized.csv`) is obtained from Zenodo record 7332091 [1] instead (see Step 1).

### Example Structure
```
SuNSTAMP/SuNSTAMP-AD/replicate_8_1/
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

2. Download the data folders from the Google Drive link above.

3. Extract and place the `data/` folders into their corresponding replicate directories.

4. Verify that the structure matches the example above.

---

## Software environment

| File | Purpose |
|---|---|
| `requirements.txt` | **Recommended installation file.** Packages imported by the SuNSTAMP scripts, with the versions recorded in the original environment. Tested: installs cleanly with `pip check` passing on macOS (arm64), Python 3.12. |
| `requirements_windows.txt` | Archival snapshot of the original development environment (a Windows `pip freeze` export, previously named `requirement.txt`), kept unchanged, including its UTF-16 encoding. It contains Windows-only and unrelated packages and is not the recommended installation file. |

`scanpy` (used by the MI preprocessing script and the Figure 4 notebooks) and `statannotations` (Figure 4 notebooks) were not recorded in the original export. `requirements.txt` therefore gives a version range compatible with the other pinned packages instead of an exact version.

---

## Pipeline Steps

The analysis pipeline is organized as a sequence of scripts inside each replicate's `src/` folder. Below, we describe each script, its inputs/outputs, and its role in the workflow. **All commands are run from inside the replicate's `src/` folder.** Each copy of a script already contains the time point and replicate number of its folder (variables `months` and `replicate`).

---

<details>
  <summary> Step-1:  extraction.py</summary>

**Purpose:**
Extracts and preprocesses a subset of the gene expression matrix corresponding to a given replicate. It selects spots (cells) by index range, transposes the matrix, and saves the result as a new CSV file for downstream analysis. The gene expression csv file should be placed in the `data` folder of the corresponding replicate.

**Inputs:**
- `expression_matrix_normalized.csv`: full normalized gene expression matrix (rows = genes, columns = spot IDs).
- Parameters inside the script: `start_index_spot` and `end_index_spot` define the range of spot IDs to extract.

**Outputs:**
- `gene_exp_13_1.csv` (example for replicate `13_1`): transposed CSV file where rows = spots, columns = genes.

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

For MI, there is a similar file named `process_adata.py`. It reads `../data/adata_{t}_{r}.h5ad` and writes `gene_exp_{t}_{r}.csv`, `spatial_{t}_{r}.csv` and `supervision_{t}_{r}.csv` to `../data/`. If you want to run this file, you need to download the original dataset from [2]. Otherwise, you may use the processed datasets from the Google Drive.
```bash
python process_adata.py      # MI only
```

</details>

---

<details>
  <summary>Step-2: spatial_decay_model.py (for Alzheimer's analysis)</summary>

**Purpose:**
Calculates a spatially smoothened supervision value for each cell based on nearby plaques, normalizes the values and initial intensity, and generates visualizations. This step is part of the supervision pipeline for SuNSTAMP in AD analysis. For MI, we directly used pathological scores from the SPaSE [3] pipeline. The python file also creates two key plots:

1. Colormap of normalized attraction values for all cells.
2. Overlay of cells and plaques with plaque radii.

**Inputs:**
- `spatial_13_1.csv`: coordinates and metadata for cells (from the `data/` folder).
- `plaque_13_1.csv`: coordinates and radii of plaques (from the `data/` folder).

**Outputs:**
- `supervision_13_1.csv` (written to `data/`), containing:
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
  <summary>Step-3: new_conda_2.sh (runs nmf_final_2.py)</summary>

**Purpose:**
Performs **k-fold nested cross-validation** for SuNSTAMP using a wide range of hyperparameters (ranks and lambda values). This step is computationally intensive and was run on an HPC cluster.

**Inputs:**
- `nmf_final_2.py`: Python script performing the NMF decomposition with cross-validation.
- Preprocessed data files from the replicate `data/` folder (e.g., `gene_exp_13_1.csv`, `supervision_13_1.csv`).
- Hyperparameter configuration inside `nmf_final_2.py`.

**Outputs:**
- Cross-validation losses for each rank and lambda combination, one file per SLURM array task: `result/mou_errors_{task_id}_{t}_{r}.csv` (task IDs 0–20). Step 4 reads these files from the subfolder `result/mou_errors_{t}_{r}/`, so move them there after the run. The saved results in this repository are already organized this way.
- Due to the heavy computation (~23 hours runtime on 16 CPUs, 128GB RAM), **these results are already provided in the `result/` folder** for readers/users.

**Usage:**

- `new_conda_2.sh` is the SLURM array job script used for the paper (`#SBATCH --array=0-20`). It contains settings specific to our cluster (`#SBATCH --chdir=...`, `module load anaconda3/2020.11`, `source /home/.../.bashrc`, the conda environment name and the Python path). **Edit these lines for your own cluster** before submitting:

```bash
sbatch new_conda_2.sh
```

- `nmf_final_2.py` reads the array index from the environment variable `SLURM_ARRAY_TASK_ID`. Without SLURM, set it explicitly for each task (0–20), e.g.:

```bash
for i in $(seq 0 20); do SLURM_ARRAY_TASK_ID=$i python nmf_final_2.py; done
```

**Notes:**

- All scripts assume the folder structure is maintained (data, src, result, figure).
- Users can skip heavy computations because the precomputed results are available in the `result` folder.

</details>

---

<details>
  <summary>Step-4: errors_grouping.py</summary>

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
- A combined CSV file with summed losses across folds: `../result/summed_loss_by_rank_lambda_{months}_{replicate}.csv`

**Usage:**
```bash
python errors_grouping.py
```

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
- MI replicates only: plot of `-log10(p-value)` vs rank, with a threshold line at **p=0.05**, saved in `../figure/`. (The AD copies of `rank_test.py` write only the CSV.)

**Usage:**
```bash
python rank_test.py
```

With the rank fixed, the next step is to explore how different lambda values (the supervision regularization parameter) affect the model. This allows us to tune the strength of supervision applied in SuNSTAMP and evaluate its stability across a wide range of configurations.

</details>

---

<details>
  <summary>Step-6: automate_W_matrix.sh</summary>

**Purpose:**
Performs **supervised NMF runs** for the chosen rank while scanning across a list of lambda values. Each run produces the corresponding **W and H matrices**, which form the core of the SuNSTAMP interpretation.

**Inputs:**
- `new_supervised_final_run.py`: Python script implementing supervised NMF for a given `(rank, lambda)` pair.
- Data folder with the expression matrix (`gene_exp_[month]_[replicate].csv`) and supervision values (`supervision_[month]_[replicate].csv`).
- The rank selected in Step 5, set in the variable `RANK` at the top of `automate_W_matrix.sh`.

**Outputs:**
- `result/W_Matrices_{t}_{r}/output_W_matrix_rank_{rank}_{lambda}_{t}_{r}.csv`: W matrices for each `(rank, lambda)` configuration.
- `result/H_Matrices_{t}_{r}/output_H_matrix_first_row_rank_{rank}_{lambda}_{t}_{r}.csv`: first row of the H matrix for each `(rank, lambda)` configuration.

**Usage:**
```bash
bash automate_W_matrix.sh
# or a single configuration:
python new_supervised_final_run.py --rank 10 --lambda 3e-5
```

</details>

---


<details>
  <summary>Step-7: Coefficients.py </summary>

**Purpose:**
In the previous step, we stored the **W matrices** and the first row of the **H matrices** across different configurations.
Now, we want to assess how well the first factor (`Supervised Factor`) from the W matrix reflects the supervision signal.

**Inputs:**
- Stored W matrices from the previous step (`result/W_Matrices_{months}_{replicate}`).
- Supervision values' file (`data/supervision_{months}_{replicate}.csv`).
- The rank to evaluate, set in the variable `rank` at the top of the script.

**Outputs:**
- `result/Supervised_Factor_correlation_results_{months}_{replicate}.csv`

**Usage:**
```bash
python Coefficients.py
```

</details>

---

<details>
  <summary>Step-8: Predictive_Scoring.py (AD)</summary>


**Purpose:**
Apply gene weights learned from a **source replicate** to a **target replicate** and evaluate their predictive power for **plaque proximity**. This quantifies the cross-replicate generalizability of SuNSTAMP factors.

### Inputs
- **Command-line arguments**
  - `top_n_genes`: number of top-ranked genes from the H-matrix to use.
  - `output_file`: CSV file path to append results (Top-N genes vs. AUC).
- **Settings at the top of the script**: `replicate_now` (target replicate), `genes_from_replicate` (source replicate), and `which_rank_from` / `which_lambda_from` (the source replicate's final rank and lambda).

**Input files**
- `gene_exp_{months}_{replicate}.csv`: target replicate gene expression matrix (cells × genes).
- `supervision_{months}_{replicate}.csv`: target replicate spatial coordinates (`X-scaled`, `Y-scaled`).
- `output_H_matrix_first_row_rank_{rank}_{lambda}_{months}_{genes_from_replicate}.csv`: source replicate H-matrix scores (gene weights), read from the source replicate's `result/H_Matrices_*` folder.
- `plaque_{months}_{replicate}.csv`: plaque positions in the tissue (`m.cx`, `m.cy`).

### Outputs
- **Result CSV**: rows with AUC-ROC scores

**Usage** (example: top 20 genes):
```bash
python Predictive_Scoring.py 20 ../result/AUC_scores.csv
```

This step evaluates **SuNSTAMP's cross-replicate prediction**: gene weights learned from a **source replicate** are used to score plaque proximity in a **target replicate**, and performance is summarized by **AUC-ROC**.

</details>

---

<details>
  <summary>Step-9: predictive_score_automate.sh and AUC_incremental.py (AD, 13-month replicates)</summary>

  **Purpose:**
Evaluate how the predictive performance of SuNSTAMP-derived gene weights varies as the number of top-ranked genes increases (from 1 to 100). This identifies the number of genes needed for maximal plaque proximity prediction in the target replicate. The manuscript shows the first 50 values.

**Usage:**
```bash
bash predictive_score_automate.sh   # writes ../result/AUC_scores.csv (header GeneCount,AUC; top 1..100 genes)
python AUC_incremental.py           # reads ../result/AUC_scores.csv, writes ../figure/auc_curve_top50.png
```

The saved results of this step are `result/AUC_scores_13_1.csv` and `result/AUC_scores_13_2.csv`.

</details>

---

## Final parameters used in the manuscript

These are the final rank and lambda for each sample (Supplementary Table S2 of the manuscript), with the corresponding saved outputs.

| Dataset | Folder | Final rank | Final lambda | Saved W / H (first row) |
|---|---|---|---|---|
| AD 8 months, rep. 1 | `SuNSTAMP-AD/replicate_8_1` | 14 | 1.75e-5 | `result/W_Matrices_8_1/output_W_matrix_rank_14_1.75e-05_8_1.csv`, `result/H_Matrices_8_1/output_H_matrix_first_row_rank_14_1.75e-05_8_1.csv` |
| AD 8 months, rep. 2 | `SuNSTAMP-AD/replicate_8_2` | 14 | 1.75e-5 | `..._rank_14_1.75e-05_8_2.csv` |
| AD 13 months, rep. 1 | `SuNSTAMP-AD/replicate_13_1` | 10 | 3e-5 | `..._rank_10_3e-05_13_1.csv` |
| AD 13 months, rep. 2 | `SuNSTAMP-AD/replicate_13_2` | 14 | 2.75e-5 | `..._rank_14_2.75e-05_13_2.csv` (folders `W_matrices_13_2`, `H_matrices_13_2`) |
| MI day 3, rep. 1 | `SuNSTAMP-MI/replicate_3_1` | 9 | 500 | `..._rank_9_500.0_3_1.csv`; also `result/final_result/W_3_1.csv`, `H_3_1.csv` |
| MI day 3, rep. 3 | `SuNSTAMP-MI/replicate_3_3` | 12 | 500 | `..._rank_12_500.0_3_3.csv`; also `result/final_result/` |
| MI day 7, rep. 2 | `SuNSTAMP-MI/replicate_7_2` | 9 | 8 | `..._rank_9_8.0_7_2.csv`; also `result/final_result/` |
| MI day 7, rep. 3 | `SuNSTAMP-MI/replicate_7_3` | 11 | 50 | `..._rank_11_50.0_7_3.csv`; also `result/final_result/` |

**Note on MI day 3, replicate 1:** `src/automate_W_matrix.sh` and `src/Coefficients.py` in this folder are still set to rank 13 from an earlier run. The final rank is 9 (see `src/steps.txt`); the rank-9 outputs are saved in `result/W_Matrices_3_1/`, `result/H_Matrices_3_1/` and `result/final_result/`, and these are the files used for Figure 4. To rerun Steps 6–7 at the final rank, set `RANK=9` / `rank = 9`.

---

## Additional scripts

| Script | Run from | Purpose | Inputs → outputs |
|---|---|---|---|
| `SuNSTAMP-AD/additional_analysis/consistent_gene_count.py` | `SuNSTAMP-AD/` | Overlap and ordering consistency of the top genes between the 13-month replicates | `replicate_13_*/result/final_gene_list_13_*.csv` |
| `SuNSTAMP-AD/additional_analysis/plot_AUC_both.py` | `SuNSTAMP-AD/` | AUC vs number of top genes for both 13-month replicates | `replicate_13_*/result/AUC_scores_13_*.csv` → `auc_curve_top50_both.png` |
| `SuNSTAMP-AD/plot_table.py` | `SuNSTAMP-AD/` | Pairwise ordering table of shared genes (Supplementary Figure S3) | `replicate_13_*/result/final_gene_list_13_*.csv` → `table_genes.png` |
| `SuNSTAMP-MI/correlation-MI_3.py` | folder containing its inputs | Applies the day-3 replicate-1 gene weights to day-3 replicate-3 data | `gene_exp_3_3.csv`, `supervision_3_3.csv`, `H_3_1.csv` |

---

## Reproducing the figures

[`FIGURE_REPRODUCTION.md`](FIGURE_REPRODUCTION.md) lists, for every panel of Figures 1–4 and Supplementary Figures S1–S4, the script, the input files, the command, and whether the panel can be regenerated from the saved results in this repository, needs the Google Drive data, or depends on an external resource.

---


## General Guidelines for Applying SuNSTAMP to Any Dataset


SuNSTAMP is designed to be flexible and can be experimented with on **any spatial observation-feature dataset**. Below are the key guidelines for adapting it:

### Required Input Files
1. **Observation-feature matrix** (`gene_exp.csv` or equivalent)
   - Rows correspond to observations (e.g., cells, samples).
   - Columns correspond to features (e.g., genes, measurements).
   - Can be any non-negative matrix; does not need to be limited to gene expression.

2. **Supervision scores** (`supervision.csv`)
   - Contains supervision values for each observation.
   - For example, these could be spatial coordinates, pathology scores, or other quantitative labels.
   - Must have an identifier column (matching the observation IDs in the observation-feature matrix).



### Cross-Validation and Parameter Selection
- SuNSTAMP uses **nested k-fold cross-validation** to select the optimal factorization **rank** and regularization **lambda**.
- Users can adjust the number of folds, rank range, and lambda values according to their dataset and experimental goals.



### Quick Manual Experiments
- To bypass nested cross-validation and experiment with **one or two specific values** of rank and lambda:
  - Use the `new_supervised_final_run.py` script, e.g. `python new_supervised_final_run.py --rank 10 --lambda 3e-5`.
  - This allows for **fast exploration** of parameters without full cross-validation.



### Recommendations
- Ensure that the **observation IDs** in both files match exactly.
- Normalize or scale your observation-feature matrix if features have widely different ranges.
- Use cross-validation when possible to avoid overfitting and ensure generalizability.

---



## Citation ##

[1]: Zeng, H., Huang, J., Zhou, H., Meilandt, W.J., Dejanovic, B., Zhou, Y., Bohlen, C.J., Lee, S.-H., Ren, J., Liu, A., et al. *Integrative in situ mapping of single-cell transcriptional states and tissue histopathology in a mouse model of Alzheimer's disease.* Nature Neuroscience, 26(3):430–446, 2023.


[2]: Calcagno, M., Taghdiri, N., Ninh, V.K., Mesfin, J.M., Toomu, A., Sehgal, R., Lee, J., Liang, Y., Duran, J.M., Adler, E., et al. *Single-cell and spatial transcriptomics of the infarcted heart define the dynamic onset of the border zone in response to mechanical destabilization.* Nature Cardiovascular Research, 1(11):1039–1055, 2022.


[3] Rahman, M.N., Abrar, M.A., Shaw, V.R., Martin, J.F., Rahman, M.S., Samee, M.A.H. *SPASE: Spatially resolved pathology scores using optimal transport on spatial transcriptomics data.* Cell Systems, 2025.
