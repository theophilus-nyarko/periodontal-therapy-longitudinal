# periodontal-therapy-longitudinal
Longitudinal analysis of surgical vs. non-surgical periodontal therapy using GLS and Linear Mixed-Effects Models

# Surgical vs. Non-Surgical Periodontal Therapy: A Longitudinal Analysis

Comparing how periodontal pocket depth (PPD) changes over 14 months under surgical and non-surgical therapy, using generalized least squares and linear mixed-effects models for repeated measures.

**Course project, Western Michigan University (April 2026)**




## Problem

Periodontitis destroys the tissues that support the teeth and, if untreated, leads to tooth loss. Surgical therapy is more expensive than non-surgical therapy, and earlier studies suggest its short-term advantage may fade over time. This project asks three questions:

1. Does the rate of change in PPD differ between surgical and non-surgical therapy?
2. Is surgical therapy more or less effective in patients who are current smokers with moderate-to-heavy alcohol use?
3. Do periodontal conditions improve the prediction of PPD beyond demographic factors?

## Data

- Randomized clinical trial of **182 patients** with moderate to severe periodontitis
- Outcome: average PPD across each patient's teeth (mm)
- **8 visits:** baseline, then every 2 months up to 14 months
- Covariates: age, race, BMI, HbA1c, smoking, alcohol use, and periodontal conditions (number of teeth, plaque index, gingival index, bleeding on probing, CRP, rheumatoid factor)
- Balanced design with no missing visits

> The dataset is not included in this repository. Only the analysis code and report are shared.

## Methods

- **Exploratory analysis:** mean and individual (spaghetti) PPD trajectories, baseline distributions, and the residual autocorrelation function (ACF)
- **Time trend:** linear vs. quadratic trends compared with a likelihood ratio test and AIC
- **Covariance structure:** the residual ACF did not decay as AR(1) would predict, so an **unstructured** covariance was used in GLS
- **Model selection:** covariates added one at a time and kept only if they lowered the AIC; GLS compared with a random-intercept, random-slope LME model; models compared under ML, with final estimates under REML
- **Effect modification:** a treatment × smoking × alcohol interaction
- **Prediction:** nested models compared with a likelihood ratio test
- **Diagnostics:** residual histograms, Q-Q plots, and residuals vs. fitted values

## Key Findings

- PPD follows a **quadratic trajectory**: it falls from baseline to about 10 months and then rebounds in both groups.
- The **rate of change differs significantly between therapies** (p = 0.0013). Non-surgical therapy produced lower PPD than surgical therapy, after adjusting for age.
- **Age** was the only demographic factor that improved model fit, and adjusting for it changed the treatment-by-time estimate.
- There was **no evidence** that surgical therapy works differently for smokers with moderate-to-heavy alcohol use (p = 0.73).
- Adding the **plaque index** improved model fit (likelihood ratio test p = 0.027); other periodontal conditions did not.
- GLS with an unstructured covariance fit better than the LME model.



## Limitations

- Averaging PPD across teeth hides tooth-level variation, even though periodontitis is site-specific.
- Fourteen months may be too short to capture long-term differences, especially given the late rebound in PPD.
- Model selection by AIC and significance may exclude clinically relevant variables.

## Repository Structure

```
├── README.md
├── code/
│   └── periodontal_analysis.R
├── report/
│   └── periodontal_therapy_report.pdf
└── figures/
    └── mean_ppd_by_treatment.png
```

## Tools

R: `nlme` (GLS, LME), `lattice`, `reshape2`, `dplyr`, `ggplot2`
