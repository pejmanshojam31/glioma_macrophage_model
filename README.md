
# Glioma-macrophage-interaction-model

This Repository depicts the codes used for the manuscript below:

[(https://www.nature.com/articles/s41540-024-00478-7)](https://www.nature.com/articles/s41540-024-00478-7)

# Framework of the Study

This study employs a two-phase approach to investigate the role of tumor-associated macrophages (TAMs) in glioma recurrence after surgical resection. The aim is to understand the predictive power of each feature in a virtual reality environment. We want to investigate to what extent localized biopsies are essential and how much more data they can provide than conventional random biopsies. This study involves two phases. Phase one develops a novel spatiotemporal mathematical model based on the interplay between the gliomas and macrophage phenotypes based on the role of oxygen in the tumor microenvironment. Phase two investigates the generated dataset to predict the tumor recurrence after the resection. Here, we explicitly define these two steps.

<img src="https://github.com/user-attachments/assets/09626458-9041-4686-ac28-4b0c4c9dc0db" alt="Study Framework" width="600">



## Phase 1: Mathematical Modeling

- Developed a spatio-temporal mechanistic model based on partial differential equations (PDEs).
- Simulates interactions between glioma cells, TAMs, and oxygen levels in the tumor microenvironment (TME).
- Incorporates the "go-or-grow" dichotomy, representing the trade-off between glioma cell proliferation and migration.
- Derives two clinically relevant observables:
  - **Tumor Infiltration Width (IW):** Distance between the tumor core (80% tumor cell density) and tumor edge (2% tumor cell density).
  - **Tumor Size (TS):** Calculated by integrating the spatial profile of tumor density normalized by its maximum value.
 
<img src="https://github.com/user-attachments/assets/0ed49bc0-474c-4747-9df6-1306417c53b7" alt="Study Framework" width="400">

### Sensitivity analysis

We also performed a time-dependent sensitivity analysis of the model based on the LHS-PRCC method. Based on the results, we chose the most important parameters involved in the model.

## Phase 2: Machine Learning Analysis

- Generated a synthetic dataset of 20,000 virtual patient profiles using the mathematical model and Latin Hypercube Sampling (LHS).
- Simulated various biopsy techniques (localized and non-localized) and imaging data (T1Gd and FLAIR MRI) volume.
- Gaussian noise was added to replicate real-world variability in the dataset.
- Applied machine learning regression models to predict **IW** at 3 and 12 months post-resection.
- Evaluated the predictive power of multiple features:
  - **Ki67 index**
  - **Macrophage ratios** at the tumor core and edge
  - **Tumor volumes** derived from T1Gd and FLAIR MRI
- Investigated the impact of data quality (e.g., noise) on model accuracy and robustness.

## Study Goal

The primary objective is to:
1. Understand the contribution of TAMs to glioma recurrence.
2. Evaluate the predictive value of various diagnostic features for developing personalized treatment strategies.

## How to use this repository:

First we can check the mathematical folder to generate the data. We already have the generated data of 20,000 synthetic patients. Therefore, we can just run the main file to get the raw dataset for both initial tumor growth and the tumor recurrence post resection.
After obtaining the dataset, you can run the file Selected data extracted to generate the dataset for the machine learning predictions similar to what we have in [here](https://github.com/pejmanshojam31/glioma_macrophage_model/tree/main/data). 

For the machine learning approaches, you can simply run the notebooks to obtain the results. If you have questions, please [contact](pejman.shojaee@tu-dresden.de) me.
