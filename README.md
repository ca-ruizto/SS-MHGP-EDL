# Surrogate modeling of nanowire stress–strain curves

Research code and saved analysis results accompanying **Heteroskedastic multivariate Gaussian process regression for surrogate modeling of stress-strain curves**.

The analyses use atomistic stress–strain data for Mo–Nb–Ta–V–W alloy nanowires under tension and compression. The R code represents curves using log strain at failure and six Legendre coefficients, fits Gaussian process surrogates, propagates predictive uncertainty to curves and mechanical properties, and evaluates active-learning strategies. Python notebooks provide neural-network comparisons, including the composition- and loading-aware mixture-of-experts model.

## Repository contents

| Location | Contents |
| --- | --- |
| `ss_data_all.RData` | Saved atomistic stress–strain data used by the R notebooks. |
| `dat_filtered.RData` | Additional saved, filtered data. |
| `Functions/` | Gaussian process fitting and prediction, covariance kernels, coordinate transforms, curve bases, Monte Carlo prediction, and IMSPE utilities. |
| `Model Fitting/` | R Markdown notebooks for fitting models and running active learning. |
| `Results/ss_all_ort_coefs_*_uns.csv` | Curve coefficients and derived properties for compression (`c`) and tension (`t`). |
| `Results/ss_ter_ort_coefs_*_uns.csv` | Additional higher-order-alloy coefficient data. |
| `Results/Fixed Data/` | Saved model fits and active-learning selection records. |
| `Results/Predictions_*.RData` | Cached heteroskedastic/homoskedastic predictions for diagonal and low-rank output covariance specifications. |
| `Results/NN_Experts/` | Neural-network prediction tables read by the R comparison notebook. |
| `Results/R1C2_*.csv` | Saved predictive-density, coverage, and interval-width comparisons. |
| `Results/AL_strategy_comparison/` | Saved active-learning strategy comparison results. |
| `Figures/` | R Markdown analysis and plotting notebooks. |
| `Python_code/Data/` | Neural-network input tables and results from several experiment variants. |
| `Python_code/NN Code/` | Neural-network training and prediction notebooks. |
| `Python_code/export_import_NN_Data.Rmd` | R/Python data exchange and prediction-import workflow. |

## Starting points

- **Main R analysis:** `Figures/Predictions_19-05 - Copy.Rmd`. This includes prediction caches, model comparisons, mechanical properties, representative curves, and later revision analyses. Its filename is historical.
- **Revision plots:** `Figures/rev_fig.Rmd` and `Figures/rev_fig2.Rmd`. These include composition summaries, curve-error summaries, uncertainty comparisons, and representative curves.
- **Curve preparation:** `Figures/idn_failure_point_all_dat.Rmd`.
- **Model fitting:** `Model Fitting/odd_bin_models_AL.Rmd` and `Model Fitting/odd_bin_models_gen_sig_AL.Rmd`.
- **Current four-group neural-network notebook:** `Python_code/NN Code/MVGPR_NN.ipynb`. It models the mean and variance of log strain at failure and of the six curve coefficients as four output groups.

Other dated notebooks and neural-network experiment folders are retained for provenance. They are not interchangeable configurations of the final analysis.

## Software

Use R with RStudio or another R Markdown-capable environment, and Python with Jupyter. The repository does not contain a package lockfile or a complete frozen environment. Consult the manuscript and associated submission archive for the reported software environment.

The R sources use the following non-base packages across the main and exploratory workflows:

```r
install.packages(c(
  "here", "rmarkdown", "knitr", "tidyverse", "scales", "ggnewscale",
  "AlgDesign", "mgcv", "gridExtra", "nloptr", "latex2exp", "collapse",
  "lme4", "Matrix", "mnormt", "rio", "corrplot", "DiceDesign",
  "compositions", "cubature", "fdapace", "ggh4x", "hetGP", "kableExtra",
  "matrixStats", "matrixcalc", "nortest", "rdist"
))
```

The Python notebooks import NumPy, pandas, PyTorch, scikit-learn, Optuna, and joblib. For a Python environment in which these are not already installed:

```sh
python -m pip install numpy pandas torch scikit-learn optuna joblib jupyter
```

These commands install available package versions; they do not reconstruct the original research environment.

## Working with the saved R analyses

1. Open `SS-MHGP-EDL.Rproj` in RStudio so that `here::here()` resolves paths from this repository.
2. Open `Figures/Predictions_19-05 - Copy.Rmd` and run the setup, function-loading, and data-preparation chunks before the analysis chunks that depend on them.
3. Use the supplied prediction caches and model files where the notebook provides a cache-loading branch. Some later revision chunks also load those caches directly.
4. Review the inputs and output paths of each section before running it. Notebook sections can refit models, regenerate predictions, or overwrite saved results. Full training and Monte Carlo evaluation can take substantial time.
5. Open the revision plotting notebooks for the corresponding additional figures and diagnostics. Supply their external revision files as described below before running those sections.

Saved filenames containing `2026_05_19` are still read by the main analysis. The low-rank comparisons also use `2026_09_21` fits. A date in a filename does not by itself indicate that the result is obsolete. Preserve the supplied model fits unless deliberately rebuilding the relevant experiment.

For a new GP training run, inspect the parameter block in the fitting notebook first. It sets the experiment name, covariance kernel, mean/covariance specification, and other run controls. The checked-in parameter block represents one configuration, rather than a command that automatically reproduces every model in the paper. Use a new experiment name/output location to preserve the supplied fits.

## Neural-network workflow

Open `Python_code/NN Code/MVGPR_NN.ipynb` in Jupyter and change `DATA_DIR` and `OUTPUT_DIR` in its configuration cell to your own locations. The input directory should contain:

```text
dat_bin_odd_c.csv
dat_bin_odd_t.csv
dat_test_bin_even_c.csv
dat_test_bin_even_t.csv
dat_test_qq_c.csv
dat_test_qq_t.csv
```

Copies with these names are provided in `Python_code/Data/`. The training tables include composition fractions, an alloy identifier, and the `m_*` and `var_*` response columns listed in the notebook. Loading type is added when the compression and tension tables are read.

Run the notebook cells in order. The notebook selects CPU or CUDA when available and writes predictions, metrics, tuning records, and fitted model/preprocessing files under `OUTPUT_DIR`. Review the checked-in trial counts, ensemble size, cross-validation folds, epoch limit, and early-stopping settings before running; these settings should not be assumed to reproduce a particular saved experiment without comparison to its training records.

The `M-Grid-1`, `V-Grid-1`, `M-Fixed-Branch-1`, and `V-Fixed-Branch-1` notebooks are earlier neural-network experiments. Several use different data layouts and machine-specific paths. Keep their results separate from the four-group notebook's outputs.

## External inputs and local paths

The repository is a research working tree, not a single-command reproduction package. Some revision sections point to a separate `SS-MHGP-EDL-revision` directory using absolute Windows paths. They require additional files such as:

- `ss_data_all_revision.RData`;
- revision coefficient tables;
- `Results/Revision Tables/` curve metrics and selected-case tables;
- `Results/Revision Evidence/` evaluation and integrated-variance records.

These files are not included in this checkout. Supply them from the associated revision/submission archive and update `revision_root`, `revision_dir`, or equivalent variables. Also change `figure_dir`, neural-network directory settings, and any other absolute paths before running those sections. The scripts can otherwise write to the original author's desktop or project folders.

Related atomistic simulation resources are available in [Binary-Wire_2026](https://github.com/shuozhixu/Binary-Wire_2026). The accompanying article and its Data S1 archive describe the submitted data and code collection.

## Output maintenance

Superseded rendered prediction reports, orphaned figure-export folders from earlier drafts, and unused preliminary prediction caches were removed during repository cleanup. Source notebooks, input data, model fits, current prediction caches, and revision result tables were retained. Regenerate figures from the appropriate notebook rather than treating an older rendered report as the final paper's source of truth.

## Citation

Please cite the accompanying article when using this code or the associated results. Use the final journal citation and article DOI once available.
