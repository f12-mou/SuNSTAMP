# Figure reproduction guide

This guide maps every panel of the manuscript's Figures 1–4 and Supplementary Figures S1–S4 to its script, input files and command, where these exist. Not every panel can be regenerated end to end from this repository: some depend on external inputs (SPaSE outputs, Enrichr, AlzMAP) or are provided as images without a plotting script. These cases are marked below and summarized under *Limitations*. **No panel requires retraining SuNSTAMP:** the final W and H matrices for every analysed sample, at the rank and lambda reported in Supplementary Table S2, are saved in the replicate `result/` folders (see the *Final parameters* table in `README.md`).

**Status labels**
- **Repo**: all inputs are included in this repository.
- **Drive**: needs the `data/` files from the Google Drive folder linked in `README.md`.
- **Re-plot**: the script plots values that are written into the script (copied from an earlier computation or a web tool); it does not recompute them.
- **External**: produced with an external tool or resource; no script.
- **Schematic**: illustration drawn in Figma; no code.

**Before running:** install the environment (`pip install -r requirements.txt`). Unless otherwise specified, run each figure-generation script from its corresponding panel directory (e.g. `cd Figure_Generation/Figure-3/Figure-3a && python Figure-3a.py`). Scripts write their figure into the same directory.

---

## Figure 1
| Panel | Status | Files |
|---|---|---|
| 1 (pipeline overview) | Schematic | `Figure_Generation/Figure-1/final_full_pipeline.svg` |

## Figure 2 (AD, 13 months, replicate 1)
| Panel | Status | Script | Inputs (place in `<panel folder>/data/` unless stated) | Command (from the panel folder) |
|---|---|---|---|---|
| 2a | Schematic | — | `Figure-2a.png` | — |
| 2b | Drive | `Figure-2b.py` | `data/spatial_13_1.csv`, `data/plaque_13_1.csv` (from `SuNSTAMP-AD/replicate_13_1/data/` on Drive) | `python Figure-2b.py` → `Figure-2b.png` |
| 2c | Re-plot | `Figure-2c.py` | Pearson correlations written in the script; the underlying values are listed in `Figure-2c-Pearson_results.txt` | `python Figure-2c.py` → `comparison_page_1_5.png` |
| 2d | Drive + Repo | `Figure-2d.py` | `data/supervision_13_1.csv` (Drive); `data/W_13_1.csv` = final W matrix of replicate 13_1 (rank 10, λ = 3×10⁻⁵), saved in this repository as `SuNSTAMP-AD/replicate_13_1/result/W_Matrices_13_1/output_W_matrix_rank_10_3e-05_13_1.csv`; copy and rename | `python Figure-2d.py` → `Figure-2d.png` |
| 2e | Drive | `Figure-2e.py` | `data/supervision_13_1.csv`, `data/spatial_13_1.csv` | `python Figure-2e.py` → `Figure-2e.png` |
| 2f | Drive + Repo | `Figure-2f.py` | `data/W_13_1.csv` (as for 2d), `data/spatial_13_1.csv` | `python Figure-2f.py` → `Figure-2f.png` |

`supervision_13_1.csv` can also be regenerated from `spatial_13_1.csv` and `plaque_13_1.csv` with `SuNSTAMP-AD/replicate_13_1/src/spatial_decay_model.py` (README Step 2).

## Figure 3 (AD)
| Panel | Status | Script | Inputs | Command |
|---|---|---|---|---|
| 3a (top-10 gene weights) | Repo | `Figure-3a/Figure-3a.py` | `Figure-3a/H_genes.csv` (first row of H for replicate 13_1, rank 10, λ = 3×10⁻⁵; identical to `SuNSTAMP-AD/replicate_13_1/result/H_Matrices_13_1/output_H_matrix_first_row_rank_10_3e-05_13_1.csv`) | `python Figure-3a.py` → `Figure-3a.png` |
| 3b (8 vs 13 months KDE) | Drive | `Figure-3b/Figure-3b.py` | `gene_exp_8_1.csv`, `gene_exp_13_1.csv` (Drive; place in the panel folder or edit `file_path_1`/`file_path_2`) | Set `gene_name` (Cst3, Ctsb, Gfap, Apoe), then `python Figure-3b.py` → `<gene>.jpg`; one gene per run |
| 3c (top-20 overlap, Venn) | Repo (lists) | no plotting script | Gene lists: `SuNSTAMP-AD/replicate_13_1/result/final_gene_list_13_1.csv`, `.../replicate_13_2/result/final_gene_list_13_2.csv`. Overlap and pair-ordering counts: `SuNSTAMP-AD/additional_analysis/consistent_gene_count.py` (run from `SuNSTAMP-AD/`) | Diagram image: `Figure-3c/Figure-3c.png` |
| 3d (cross-replicate AUC-ROC) | Drive (scores only) | no plotting script | Computed with `SuNSTAMP-AD/replicate_13_*/src/Predictive_Scoring.py` (README Step 8; needs `gene_exp`, `supervision`, `plaque` of the target replicate from Drive and the source replicate's saved H). Saved summary: `SuNSTAMP-AD/replicate_13_1/result/Predictive_Scores.txt` | Diagram image: `Figure-3d/Figure-3d.png` |
| 3e (AUC vs number of top genes) | Repo | `SuNSTAMP-AD/additional_analysis/plot_AUC_both.py` (both replicates); `Figure-3e/Figure-3e.py` (one replicate) | `SuNSTAMP-AD/replicate_13_1/result/AUC_scores_13_1.csv`, `.../replicate_13_2/result/AUC_scores_13_2.csv` | From `SuNSTAMP-AD/`: `python additional_analysis/plot_AUC_both.py` → `auc_curve_top50_both.png`. `Figure-3e.py` reads `AUC_scores.csv` in its own folder: copy `AUC_scores_13_1.csv` there under that name |
| 3f (AlzMAP expression) | External | — | Images from AlzMAP (https://alzmap.org) | `Figure-3f/Figure-3f.png` |
| 3g (Enrichr, AD) | Re-plot | `Figure-3g/Figure-3g.py` | Enrichr results written in the script; the gene list submitted to Enrichr (top 20) is `Figure-3g/Top_Genes_AD_13_1.csv` | `python Figure-3g.py` → `Figure-3g.png`. Enrichment itself: Enrichr web tool (https://maayanlab.cloud/Enrichr/) |

## Figure 4 (MI)
| Panel | Status | Script | Inputs | Command |
|---|---|---|---|---|
| 4a (SPaSE pathological scores) | External | `Figure-4a/Figure-4a-SPaSE_Figures.ipynb` | SPaSE output files read from absolute paths on the original analysis machine; not included in this repository | Not runnable as distributed (see *Limitations*) |
| 4b (supervised factor vs SPaSE score) | External | `Figure-4b/Figure-4b-SunSTAMP-MI_figures.ipynb` | SPaSE score files read from absolute paths on the original analysis machine; the final MI W matrices are in `SuNSTAMP-MI/replicate_*/result/final_result/` | Not runnable as distributed (see *Limitations*) |
| 4c (top-10 genes, day 3) | Repo | `Figure-4c/Figure-4c.py` | `Figure-4c/H_3_1.csv` (identical to `SuNSTAMP-MI/replicate_3_1/result/final_result/H_3_1.csv`, rank 9, λ = 500) | `python Figure-4c.py` → `Figure-4c.png` |
| 4d (top-10 genes, day 7) | Repo | `Figure-4d/Figure-4d.py` | `Figure-4d/H_7_2.csv` (identical to `SuNSTAMP-MI/replicate_7_2/result/final_result/H_7_2.csv`, rank 9, λ = 8) | `python Figure-4d.py` → `Figure-4d.png` |
| 4e (overlap, day 3, Venn) | Repo (lists) | no plotting script | `SuNSTAMP-MI/replicate_3_1/result/final_result/H_3_1.csv`, `.../replicate_3_3/result/final_result/H_3_3.csv` | Diagram image: `Figure-4e/Figure-4e.png` |
| 4f (overlap, day 7, Venn) | Repo (lists) | no plotting script | `SuNSTAMP-MI/replicate_7_2/result/final_result/H_7_2.csv`, `.../replicate_7_3/result/final_result/H_7_3.csv` | Diagram image: `Figure-4f/Figure-4f.png` |
| 4g (Enrichr, day 3) | Re-plot | `Figure-4g/Figure-4g.py` | Enrichr results written in the script; ranked genes `Figure-4g/H_3_1.csv` (top 15 submitted) | `python Figure-4g.py` → `Figure-4g.png` |
| 4h (Enrichr, day 7) | Re-plot | `Figure-4h/Figure-4h.py` | Enrichr results written in the script; ranked genes `Figure-4h/H_7_2.csv` (top 15 submitted) | `python Figure-4h.py` → `Figure-4h.png` |

## Supplementary figures and tables
| Item | Status | Script | Inputs | Command |
|---|---|---|---|---|
| S1, S2 | Schematic | — | `Figure-S1.png`, `Figure-S2.png` | — |
| S3 (pairwise ordering of shared genes) | Repo | `Supplementary_Figures/Figure-S3/Figure-S3.py` (same as `SuNSTAMP-AD/plot_table.py`) | `SuNSTAMP-AD/replicate_13_*/result/final_gene_list_13_*.csv` | Paths are relative to `SuNSTAMP-AD/`: `cd SuNSTAMP-AD && python plot_table.py` → `table_genes.png` |
| S4 (AlzMAP expression) | External | — | Images from AlzMAP (https://alzmap.org) | `Figure-S4.png` |
| Table S2 (final rank and lambda) | Repo | — | Rank and lambda appear in the saved file names; see the *Final parameters* table in `README.md` | — |

---

## Summary
| Category | Panels |
|---|---|
| Regenerable from files in this repository | 3a, 3e, 4c, 4d, S3 |
| Re-plot of values written in the script | 2c, 3g, 4g, 4h |
| Regenerable after downloading the Google Drive `data/` folders | 2b, 2d, 2e, 2f, 3b |
| Underlying values recomputable (3d: after the Drive download) but no plotting script; panel provided as an image only | 3c, 3d, 4e, 4f |
| Depend on external tools or resources; not regenerable from this repository | 3f, S4 (AlzMAP); 4a, 4b (SPaSE outputs not distributed) |
| Schematics | 1, 2a, S1, S2 |

## Limitations
- **Figure 4a/4b notebooks** read SPaSE output files through absolute paths on the machine where the analysis was run; these files are not distributed here. SPaSE [Rahman et al., Cell Systems, 2025] is a separate tool; the SPaSE settings used for each MI replicate appear in the score file names inside the notebooks (e.g. `scores_D3_1_alpha_0.0001_lambda_0.1.csv`). The MI supervision values used by SuNSTAMP are the `supervision_{t}_{r}.csv` files written by `process_adata.py`.
- **Venn diagrams (3c, 4e, 4f) and panel 3d** are provided as images; there is no plotting script, but the underlying gene lists or scores are listed above.
- **Enrichment panels (3g, 4g, 4h)** were computed with the Enrichr web tool; the scripts re-plot the reported values.
- The commands above were documented from the scripts' file paths. They were not re-executed when this guide was written.
