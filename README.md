# Incidence-Based Mortality (IBM) Rate and Rate Ratio

This repository contains R functions to compute **Incidence-Based Mortality (IBM) rates** and **rate ratios**.  
The methods account for small-count bias and provide variance estimates suitable for rare cancers or small populations.

---

##  IBM Rate
The **IBM rate** is calculated by linking deaths to incident cancer cases in the registry.  
- Adjustments for small counts are available via:
  - **Fay–Feuer method** (recommended for rare events)
  - **Tiwari’s modification** (alternative small-count adjustment)

Rates can be age-adjusted as needed.

---

##  Rate Ratios
Rate ratios compare IBM rates across groups (e.g., sex, race, or calendar period).  
- Variance estimation is performed using the **Delta method**,  
  which approximates the standard error of the log rate ratio.  
- Confidence intervals are computed on the log scale and then exponentiated.

---

## Resources
- [Paper](TBA)  
- [R code](NCI/IBM.R)  

---

## ⚙️ Example Usage

```r
# Load functions
source("NCI/IBM.R")

# Example: Compute IBM rate with Fay–Feuer method
ibm_rate <- compute_dsr_and_rr_for_subset(df, idx1, idx2, "ER- & NHW & 30-54","ER-& NHB & 30-54","Subset", ci_method = "fayfeuer")





## Models for health disparity analysis 

#### Qauntile forward regression for high-dimensional survival data (2023)
 * Our main objective is to develop a quantile-specific prediction model using high-dimensional data. To demonstrate the effectiveness of this model, we present two examples. The first example focuses on the quantile-specific sequential selection of dietary factors associated with BMI. The second example highlights the quantile-specific sequential selection of risk factors for time-to-event analysis. Importantly, our model incorporates race and the interaction between race and risk factors to capture their impact on the predictions.
  * [paper](TBA) &nbsp; &nbsp;&nbsp; &nbsp; [R code](QFR/demo.R)

#### Varying-coefficients for regional quantile via KNN-based LASSO with applications to health outcome study (2023)
 * This novel framework for dynamic modeling of the associations
between health outcomes and risk factors can capture
the time-varying effects of age.
  * [paper](TBA) &nbsp; &nbsp;&nbsp; &nbsp; [MATLAB code](https://github.com/younghhk/software/tree/master/KNN)
  


#### A quantile regression decomposition estimation of disparities for complex survey data (2023+)
 *  The method proposes a quantile regression (QR) decomposition approach to disparity research using complex survey data. 
 * [paper](TBA) &nbsp; &nbsp;&nbsp; &nbsp; [R code](R/PBQR_github.R)
   


 
## COVID-19 infectious disease modeling
### Time-varying Poisson SIR model (2020)
* [paper](https://arxiv.org/abs/2004.05730)  &nbsp; &nbsp;&nbsp; &nbsp;   [R code](R/tvpSIR.R)

## Survival data screening in high-dimensional data 

 
#### Feature selection of ultrahigh-dimensional covariates with survival outcomes: a selective review (2017)
  * [paper](https://www.stt.msu.edu/users/hhong/review_survival_high.pdf)
  
 
#### The Lq-norm learning for ultrahigh-dimensional survival data: an integrated framework (2018)
* The Lq-norm learning is proposed to detect predictors with various levels of impact, such as short- or long-term impact, on censored
outcome.
 * [paper](https://www.stt.msu.edu/users/hhong/2018-CMC-0715-4p.pdf) &nbsp; &nbsp;&nbsp; &nbsp;   [R code](R/Lq.R)
  
#### Integrated powered density (IPOD): screening ultrahigh dimensional covariates with survival outcome (2018)
 *  With a flexible weighting scheme, Kolmogorov statistic as a special case,  IPOD method can detect early or late impact on censored outcome.
   * [paper](https://www.stt.msu.edu/users/hhong/Hong_et_al-2017-Biometrics.pdf)  &nbsp; &nbsp;&nbsp; &nbsp;  [R code](R/IPOD.R)
 
#### Conditional screening for survival data (2018)
 * The recently developed variable screening methods, though powerful in many practical setting,  are less powerful in detecting marginally weak while jointly important signals. A new conditional screening method for survival outcome data computes the marginal contribution of each biomarker given priorly known biological information.
  * [paper](https://www.stt.msu.edu/users/hhong/conditional_survival.pdf)  &nbsp; &nbsp;&nbsp; &nbsp; [R code](R/CS.R)
  
#### Quantile adaptive model-free variable screening for high-dimensional heterogeneous data (2013)
 * The proposed nonlinear independence screening procedure employs spline approximations to model the marginal effects at a quantile level of interest.
 * [paper](https://www.stt.msu.edu/users/hhong/screening.pdf) &nbsp; &nbsp;&nbsp; &nbsp; [R code](R/QA.R)

