# Quantile forward regression for high-dimensional survival data (2023)



## Overview

Our main objective is to develop a quantile-specific prediction model using high-dimensional data. 

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

   [**Eun Ryung Lee**](https://sites.google.com/view/eunryunglee/home), [**Seyoung Park**](https://sites.google.com/view/seyoungpark/home),

 Department of Statistics, Sungkyunkwan University 
 
  [**Sang Kyu Lee**](),  and [**Hyokyoung (Grace) Hong**](https://dceg.cancer.gov/about/staff-directory/hong-grace)
 
 Division of Cancer Epidemiology & Genetics, National Cancer Institute, NIH


* For any questions or inquiries related to the code, please reach out to us at ishspsy@skku.edu.

