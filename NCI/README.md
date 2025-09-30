# SOFTWARE for SEER Cancer Data Analysis

## Incidence-Based Mortality (IBM) Rate and Rate Ratio

This repository provides R functions to compute **Incidence-Based Mortality (IBM) rates** and **rate ratios**.  
The methods include adjustments for small-count bias and variance estimation, making them suitable for rare cancers and small populations.

---

##  IBM Rate
The **IBM rate** links deaths to incident cancer cases in the registry.  
Adjustments for small counts are available:
- **Fay–Feuer method** (recommended for rare events)  
- **Tiwari’s modification** (alternative adjustment)  

Rates can be age-adjusted if needed.

---

##  Rate Ratios
Rate ratios compare IBM rates across groups (e.g., sex, race, calendar period).  
- Variances are estimated with the **Delta method**, which approximates the standard error of the log rate ratio.  
- Confidence intervals are computed on the log scale and then exponentiated.

---

##  Resources

- [Age-adjusted Rate Confidence Intervals (SEER Documentation)](https://seer.cancer.gov/help/seerstat/equations-and-algorithms/rate-algorithms)

- **Rate Ratios**  
  - Confidence interval formula:  
    Fay MP. *Approximate confidence intervals for rate ratios from directly standardized rates with sparse data.*  
    Communications in Statistics: Theory and Methods. 1999; 28(9):2141–2160.  
  - P-value formula:  
    Fay MP, Tiwari RC, Feuer EJ, Zou Z. *Estimating average annual percent change for disease rates without assuming constant change.*  
    Biometrics. 2006; 62(3):847–854.

- [R code](NCI/IBM.R)


---

##  Example Usage

```r
# Load functions
source("NCI/IBM.R")

# Example: Compute IBM rate with Fay–Feuer method
ibm_rate <- compute_dsr_and_rr_for_subset(
  df, idx1, idx2,
  "ER- & NHW & 30-54",
  "ER- & NHB & 30-54",
  "Subset",
  ci_method = "fayfeuer"
)
```
