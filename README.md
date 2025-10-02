
## Table of Contents
- [Quantile Forward Regression for High-Dimensional Survival Data (2023)](#qfr-2023)
- [Varying-Coefficient Regional Quantile via KNN–LASSO (2023)](#vc-knn-lasso-2023)
- [Quantile Regression Decomposition for Complex Surveys (2024)](#qr-decomp-2023)
- [COVID-19 Infectious Disease Modeling](#covid19)
- [Selective Review (2017)](#review-2017)
- [Lq-Norm Learning (2018)](#lq-2018)
- [IPOD Screening (2018)](#ipod-2018)
- [Conditional Screening (2018)](#cs-2018)
- [Quantile-Adaptive Model-Free Screening (2013)](#qa-2013)

---

<a id="qfr-2023"></a>
- **Quantile Forward Regression for High-Dimensional Survival Data (2023)**  
  **Summary:** Quantile-specific prediction for censored outcomes with sequential feature selection; allows effects to differ across survival quantiles; includes race and interactions for disparity analysis.  
  **Published record:** Lee, E. R., Park, S., Lee, S. K., & Hong, H. G. (2023). *Lifetime Data Analysis*, 29(4), 769–806. [DOI](https://doi.org/10.1007/s10985-023-09603-w) · [PMID](https://pubmed.ncbi.nlm.nih.gov/37393569/) · [PMCID](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12128792/)  
  **Code:** [R](QFR/demo.R)

<a id="vc-knn-lasso-2023"></a>
- **Varying-Coefficient Regional Quantile via KNN–LASSO (2023)**  
  **Summary:** Combines KNN-based local smoothing with LASSO to estimate region-/age-varying effects at targeted quantiles; captures time-varying associations without a single global effect.  
  **Published record:** Park, S., Lee, E. R., & Hong, H. G. (2023). *Statistics in Medicine*, 42(22), 3903–3918. [DOI](https://doi.org/10.1002/sim.9839) · [PMID](https://pubmed.ncbi.nlm.nih.gov/37365909/) · [PMCID](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11370892/)  
  **Code:** [MATLAB](https://github.com/younghhk/software/tree/master/KNN)

<a id="qr-decomp-2023"></a>
- **Quantile Regression Decomposition for Complex Surveys (2024)**  
  **Summary:** Peters–Belson / Oaxaca–Blinder–style decomposition at outcome quantiles with design-based inference (weights/strata/clusters); partitions explained vs unexplained components.  
  **Published record:** Hong, H. G., Graubard, B. I., Gastwirth, J. L., & Kim, M.-O. (2024). *Annals of Applied Statistics*, 18(3), 2012–2033. [DOI](https://doi.org/10.1214/23-AOAS1868) · [PMID](https://pubmed.ncbi.nlm.nih.gov/40995413/) · [PMCID](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12456447/)  
  **Code:** [R](R/PBQR_github.R)

<a id="covid19"></a>
- **Estimation of Time-Varying Reproduction Numbers (COVID-19) (2020)**  
  **Summary:** Poisson SIR framework with time-varying transmission to estimate effective reproduction numbers and track epidemic dynamics.  
  **Published record:** Hong, H. G., & Li, Y. (2020). *PLOS ONE*, 15(7), e0236464. [DOI](https://doi.org/10.1371/journal.pone.0236464)  
  **Preprint & Code:** [arXiv:2004.05730](https://arxiv.org/abs/2004.05730) · [R](R/tvpSIR.R)

<a id="review-2017"></a>
- **Feature Selection of Ultrahigh-Dimensional Covariates with Survival Outcomes: A Selective Review (2017)**  
  **Summary:** Review of screening strategies for censored outcomes in ultrahigh-dimensional settings; guidance on use and pitfalls.  
  **Published record:** Hong, H. G., & Li, Y. (2017). *Statistical Methods in Medical Research*. [PMCID](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5906071/) · [PDF](https://www.stt.msu.edu/users/hhong/review_survival_high.pdf)

<a id="lq-2018"></a>
- **The Lq-Norm Learning for Ultrahigh-Dimensional Survival Data: An Integrative Framework (2018)**  
  **Summary:** Detects predictors with heterogeneous (short/long-term) impacts on censored outcomes via Lq-norm learning.  
  **Published record:** Hong, H. G., Chen, X., Kang, J., & Li, Y. (2018). *Statistica Sinica*, 30(3), 1213–1233 (2020 print). [DOI](https://doi.org/10.5705/ss.202017.0537) · [PMCID](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7394456/)  
  **Code:** [R](R/Lq.R)

<a id="ipod-2018"></a>
- **Integrated Powered Density (IPOD): Screening Ultrahigh-Dimensional Covariates with Survival Outcomes (2018)**  
  **Summary:** Weighted screening method (Kolmogorov as special case) sensitive to early vs late survival effects.  
  **Published record:** Hong, H. G., Chen, X., Christiani, D. C., & Li, Y. (2018). *Biometrics*, 74(2), 421–429. [DOI](https://doi.org/10.1111/biom.12820) · [PMCID](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6495533/)  
  **Code:** [R](R/IPOD.R)

<a id="cs-2018"></a>
- **Conditional Screening for Ultra-High Dimensional Covariates with Survival Outcomes (2018)**  
  **Summary:** Improves detection of marginally weak but jointly important biomarkers by conditioning on prior biological information.  
  **Published record:** Hong, H. G., Kang, J., & Li, Y. (2018). *Lifetime Data Analysis*, 24(1), 45–71. [DOI](https://doi.org/10.1007/s10985-016-9387-7) · [PMCID](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5494024/)  
  **Code:** [R](R/CS.R)

<a id="qa-2013"></a>
- **Quantile-Adaptive Model-Free Variable Screening for High-Dimensional Heterogeneous Data (2013)**  
  **Summary:** Nonlinear, model-free screening using spline-based marginal effects at user-chosen quantiles; retains variables with quantile-specific signal.  
  **Published record:** He, X., Wang, L., & Hong, H. G. (2013). *Quantile-adaptive model-free variable screening for high-dimensional heterogeneous data.* *Annals of Statistics*, 41(1), 342–369. [arXiv: 1304.2186](https://arxiv.org/abs/1304.2186) <br>
  **Code:** [R](R/QR.R)
