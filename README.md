
# Glioma-macrophage-interaction-model

This Repository depicts the codes used for the manuscript below:
https://www.biorxiv.org/content/10.1101/2024.06.25.600613v1.abstract 

# Framework of the Study

This study employs a two-phase approach to investigate the role of tumor-associated macrophages (TAMs) in glioma recurrence after surgical resection. The aim is to understand the predictive power of each feature in a virtual reality environment. We want to investigate to what extent localized biopsies are essential and how much more data they can provide than conventional random biopsies. This study involves two phases. Phase one develops a novel spatiotemporal mathematical model based on the interplay between the gliomas and macrophage phenotypes based on the role of oxygen in the tumor microenvironment. Phase two investigates the generated dataset to predict the tumor recurrence after the resection. Here, we explicitly define these two steps.

[Figure 1.pdf](https://github.com/user-attachments/files/18009840/Figure.1.pdf)

## Phase 1: Mathematical Modeling

- Developed a spatio-temporal mechanistic model based on partial differential equations (PDEs).
- Simulates interactions between glioma cells, TAMs, and oxygen levels in the tumor microenvironment (TME).
- Incorporates the "go-or-grow" dichotomy, representing the trade-off between glioma cell proliferation and migration.
- Derives two clinically relevant observables:
  - **Tumor Infiltration Width (IW):** Distance between the tumor core (80% tumor cell density) and tumor edge (2% tumor cell density).
  - **Tumor Size (TS):** Calculated by integrating the spatial profile of tumor density normalized by its maximum value.
  - 
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


