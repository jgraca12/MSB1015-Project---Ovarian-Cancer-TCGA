# MSB1015-Project---Ovarian-Cancer-TCGA

## Aim of this project
This project aimed to investigate if adding gene-expression data improves prediction of disease-free survival in ovarian cancer beyond the clinical features alone.

The analysis uses clinical and mRNA gene-expression data from the TCGA ovarian serous cystadenocarcinoma cohort.

## Software
MATLAB (version R2026b (26.2.0.3386108)) using Live Scripts.

## Required toolboxes
- Statistics and Machine Learning Toolbox
- Bioinformatics Toolbox

## Repository structure
- 'scripts/' - MATLAB live scrips containing the analysis
- 'results/' - generated intermediate/result files
- 'figures/' - exported figures

## Running the analysis

Run the MATLAB Live Scripts in the following order:

1. `data_processing_01.mlx`
2. `clinical_eda_02.mlx`
3. `expression_qc_03.mlx`
4. `grade_classification_04.mlx`
5. `survival_models_05.mlx`

The first script imports and cleans the datasets and saves `results/cleaned_data.mat`, which is loaded by the subsequent scripts.

## Data

The analysis uses TCGA ovarian cancer clinical and unified mRNA expression data.

Data can be obtained from cBioPortal: https://www.cbioportal.org/study/summary?id=ov_tcga_pub

Place the following files in the `data/` folder:

- `data_clinical_patient.txt`
- `data_clinical_sample.txt`
- `data_mRNA_unified.txt`

Source:
TCGA ovarian serous cystadenocarcinoma dataset, Nature 2011, DOI: 10.1038/nature10166
