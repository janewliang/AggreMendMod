# AggreMendMod

**Aggre**gate **Mend**elian **Mod**els

This repository contains code to reproduce the results presented in the paper 

> **Liang, J. W., Shannon, K. M., Bear, L. M., McCarthy, A. M., Idos, G. E., Hong, C., Gruber, S. B., Braun, D., & Parmigiani, G. (2025). Evaluating Mendelian risk prediction models that aggregate across genes and cancers.**

## Dependencies
- `Fam3PRO` (formerly PanelPRO) R package v1.1.0, Lee, G., et al. (2021)<sup>[1](#myfootnote1)</sup>. Available [here](https://projects.iq.harvard.edu/bayesmendel/panelpro).
- Other R packages: [`abind`](https://cran.r-project.org/web/packages/abind/index.html), [`knitr`](https://cran.r-project.org/web/packages/knitr/index.html), [`pROC`](https://cran.r-project.org/web/packages/pROC/index.html), [`tidyverse`](https://cran.r-project.org/web/packages/tidyverse/index.html)
- Data from the USC-Stanford Hereditary Cancer Panel (HCP) Testing Study. Idos, G., et al. (2019)<sup>[2](#myfootnote2)</sup>. 
- Data from the Center for Cancer Risk Assessment at Massachusetts General Hospital (MGH). 

## Navigation
Additional details can be found in the sub-directory README files. 

- `simulate_families/`: Functions for simulating detailed pedigrees, including family history of cancer, genotypes, and tumor marker testing results. 
- `est_pen/`: Code to estimate aggregate penetrances. 
- `simulations/`: Code to run simulation studies. 
- `validation/`: Code to validate models on the HCP and MGH cohorts. 

---

<a name="myfootnote1">1</a>. Lee, G., Zhang, Q., Liang, J. W., Huang, T., Choirat, C., Parmigiani, G., & Braun, D. (2021). PanelPRO: A R package for multi-syndrome, multi-gene risk modeling for individuals with a family history of cancer. arXiv preprint arXiv:2010.13011.

<a name="myfootnote2">2</a>. Idos, G. E., Kurian, A. W., Ricker, C., Sturgeon, D., Culver, J. O., Kingham, K. E., ... & Levonian, P. (2019). Multicenter prospective cohort study of the diagnostic yield and patient experience of multiplex gene panel testing for hereditary cancer risk. JCO Precision Oncology, 3, 1-12.

<a name="myfootnote3">3</a>. Statistical Research and Applications Branch, National Cancer Institute. (2020). DevCan: Probability of Developing or Dying of Cancer Software. Version 6.7.8.
