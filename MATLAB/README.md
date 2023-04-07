# Varying-coefficients for regional quantile via KNN-based LASSO with applications to health outcome study (2023+)



## Overview

The proposed KNN-based LASSO is a novel framework for dynamic modeling of the associations between  outcomes and  factors using varying-coefficients (VC) regional quantile regression via K-nearest neighbors (KNN) fused Lasso.  It can capture the time-varying effects of age.  We develop an alternating direction method of multipliers (ADMM) algorithm. 


## Main functions

- [demo_simulation.m](https://github.com/younghhk/software/blob/master/MATLAB/demo_simulation.m)
: Toy example running the proposed method.

- [VC_qt_knn_admm.m](https://github.com/younghhk/software/blob/master/MATLAB/KNN/VC_qt_knn_admm.m)
: Main ADMM algorithm solving optimization problem.

- [likelihood_knn.m](https://github.com/younghhk/software/blob/master/MATLAB/KNN/likelihood_knn.m)
: Computing likelihood value (used for BIC computation).

- [supplementary_code](https://github.com/younghhk/software/tree/master/MATLAB/KNN/supplementary_code)
: Directory including all the source files and functions related to KNN Lasso


## Authors

* [**Seyoung Park**](https://sites.google.com/view/seyoungpark/home),   [**Eun Ryung Lee**](https://sites.google.com/view/eunryunglee/home)
and [**Hyokyoung (Grace) Hong**](https://dceg.cancer.gov/about/staff-directory/hong-grace)

 Department of Statistics, Sungkyunkwan University and
 Division of Cancer Epidemiology & Genetics, National Cancer Institute, NIH.


## Contact

* grace.hong@nih.gov





