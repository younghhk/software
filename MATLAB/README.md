# Varying-coefficients for regional quantile via KNN-based LASSO with applications to health outcome study (2023+)



## Overview

Our KNN-based LASSO framework is a new approach for dynamic modeling of the associations between risk factors and outcomes using VC regional quantile regression via K-nearest neighbors fused Lasso. One of its strengths is its ability to capture the time-varying effects of age. To make it all work, we've created an ADMM algorithm that uses an alternating direction method of multipliers.

## Main functions

- [demo_simulation.m](https://github.com/younghhk/software/blob/master/MATLAB/demo_simulation.m)
: Toy example running the proposed method.

- [VC_qt_knn_admm.m](https://github.com/younghhk/software/blob/master/MATLAB/KNN/VC_qt_knn_admm.m)
: Main ADMM algorithm solving optimization problem.

- [likelihood_knn.m](https://github.com/younghhk/software/blob/master/MATLAB/KNN/likelihood_knn.m)
: Computing likelihood value (used for BIC computation).

- [supplementary_code](https://github.com/younghhk/software/tree/master/MATLAB/KNN/supplementary_code)
: Directory including all the source files and functions related to KNN Lasso

## Note
The  algorithms require the use of nearestneighbour.m by Richard Brown. See details in https://www.mathworks.com/matlabcentral/fileexchange/12574-nearestneighbour-m.

For KNN Lasso algorithm, we use the algorithm by Steven Siwei Ye and Oscar Hernan Madrid Padilla. See details in Non-parametric quantile regression via the K-NN fused lasso. Journal of Machine Learning Research, Vol. 22, No. 111, 1-38, 2021.

We use parametric max-flow algorithm in "On Total Variation Minimization and Surface Evolution Using Parametric Maximum Flows" by Antonin Chambolle and Jérôme Darbon (https://link.springer.com/article/10.1007/s11263-009-0238-9). 




## Authors

* [**Seyoung Park**](https://sites.google.com/view/seyoungpark/home),   [**Eun Ryung Lee**](https://sites.google.com/view/eunryunglee/home)
and [**Hyokyoung (Grace) Hong**](https://dceg.cancer.gov/about/staff-directory/hong-grace)

 Department of Statistics, Sungkyunkwan University and
 Division of Cancer Epidemiology & Genetics, National Cancer Institute, NIH.


## Contact

* ishspsy@skku.edu





