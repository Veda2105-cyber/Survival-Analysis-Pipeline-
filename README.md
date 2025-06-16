🧬 Survival Analysis Pipeline for Clinical and Multi-Omics Data

Welcome! This repository contains a highly modular and reproducible **Snakemake-based pipeline** that allows you to integrate, preprocess, and analyze **multi-omics and clinical survival datasets**. Whether you're a researcher working with TCGA datasets or a student diving into survival analysis for the first time — this pipeline is here to help you run end-to-end analysis smoothly and transparently.



🌟 What Can You Do With This?

✅ Preprocess large clinical datasets  
✅ Merge gene expression, mutation, and survival data  
✅ Normalize features & assign survival labels  
✅ Select top prognostic genes using **L1-regularized Cox models**  
✅ Visualize high-risk vs low-risk groups  
✅ Reuse or modify the workflow for your own dataset 🎯  

---

🗂 Folder Structure

```
Project scripts/
├── Snakefile                        # Master workflow
├── config.yaml                     # Editable config (file paths, parameters)
├── clean_data.py                   # Cleans raw clinical data
├── merge_data.py                   # Combines omics and survival data
├── normalize_expression.py         # Z-score scaling of expression data
├── select_top_cox_features.py     # Feature selection with L1-regularized Cox
├── train_cox_model.py             # Final model training
├── plot_heatmap.py                # Visualization
├── filter_mutations.py            # Mutation cleanup
├── label_survival.py              # Labeling survival bins
├── check_ids.py / check_final_matrix.py # Sanity checks
├── *.yaml                         # Additional plotting/model settings
```

---

 🛠️ How to Use This Pipeline

Step 1: Clone the Repo & Install Requirements

```bash
git clone https://github.com/yourusername/survival-analysis-pipeline.git
cd Project\ scripts
pip install -r requirements.txt  
```

Step 2: Add Your Data

Replace the following files in your working directory or update paths in `config.yaml`:
- `clinical.csv` — Patient data
- `expression.csv` — Gene expression matrix
- `mutation.csv` — Binary mutation matrix (optional)
- `survival.csv` — Survival times & status

Step 3: Configure Parameters

Open `config.yaml` and edit:

```yaml
input_data:
  clinical: "path/to/clinical.csv"
  expression: "path/to/expression.csv"
  mutation: "path/to/mutations.csv"
  survival: "path/to/survival.csv"

parameters:
  normalize: true
  top_k_features: 25
  cox_penalty: 0.8
```

Step 4: Run It

```bash
snakemake --cores 4
```

Use `--forceall` to overwrite all steps.

---

🧠 How It Works Behind the Scenes

 `clean_data.py`
- Handles missing values
- Converts categorical columns
- Drops irrelevant fields

`merge_data.py`, `merge_with_labels.py`
- Combines omics and clinical data
- Assigns survival status

`normalize_expression.py`
- Z-score normalizes the expression matrix

`select_top_cox_features.py`
- Trains a **Lasso-Cox (L1 penalty)** model
- Extracts top k predictors

`train_cox_model.py`
- Trains final Cox model
- Produces survival plots + evaluation metrics

`plot_heatmap.py`
- Visualizes expression of top genes across groups



🧪 Example Outputs

-  Cleaned dataset with survival labels
-  List of top survival-linked genes
-  Trained CoxPH model + coefficients
- Heatmaps and survival curves



🎯 Who Should Use This?

This is for you if you're:
- A **bioinformatics student** doing a survival project
- A **researcher** identifying survival markers
- A **data scientist** exploring omics integration methods



✨ Tips

- 🧬 Don’t have mutation data? Remove mutation-based rules from the Snakefile!
- 🧪 Want a different model? Plug in your method in `train_cox_model.py`.
- 📊 Add visualizations? Use `matplotlib` or `seaborn` at the final stage!



👩‍💻 Maintained by

[Veda Gunti](https://www.linkedin.com/in/yourlinkedin)  
Computational Biology | AI for Clinical Genomics  
📧 vedagunti@gmail.com



🌱 Built with curiosity, coffee, and the belief that clean code = kind code.
