# Package predsurv

This package includes relevant code for building relevant prognostic models in medicine.
The most usual modeling techniques are available as well as measures to asses model performance.
The main functions are wrapped in module blocks mostly from relevant packages such as prodlim, pec, and resample.
This package is perfect for visualising the relevant measures.
We will review the most relevant methods including adaptations to the use of high throughput genomic data.
There are a vast number of packages focused in survival analysis and our goal is make a summary of the most useful functions and the best ways to combine them to obtain a meaningfull analysis.  

PS:It is planned to develop also measures for model performance with time-dependent covariates. 

## Getting Started

This vignette serves as an approach to the modeling of survival data in medical statistics.
```
library(predsurv)
survdata <- surv_sim_data(N = 1000, features = 100, CenRate = 1/10)
plot_km(survdata, time = os_months, status = os_deceased)
plot_tte(survdata, time = os_months, status = os_deceased)
```

### Prerequisites

This package is build in R to install R follow R-CRAN guidelines.


Many functions wraps many important packages for survival analysis.
For learning more about survival analysis in R: R-CRAN Task survival

```
		snowball,
		purrr,
		caret,
		dplyr,
		ggplot2,
		survival,
		glmnet,
		pec
		
```
### Installing
Before installing the package make sure that you have all the dependencies shown at the previous section.
I would not recommend install the package without installing first all the dependencies because it can take long to install specially 
```
devtools::install_github("csetraynor/predsurv")
```
