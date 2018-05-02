# Package predsurv

This package includes relevant code and figures of various prognostic models in patient outcome analysis.
Advanced modeling techniques are available as well as measures to asses model performance.
The main functions are wrapped in module blocks mostly from relevant R packages such as glmnet, pec, or survival.
This package is perfect for fitting, assesing and visualizing this relevant measures.
We will review the most relevant methods including adaptations to the use of high throughput data in the vignette.


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


Many functions wraps many important packages available from R-CRAN. 
For learning more about survival analysis in R: R-CRAN Task survival

```
		purrr,
		caret,
		dplyr,
		rlang,
		magrittr,
		assertthat,
		ggplot2,
		ipred,
		survival,
		tdROC,
		glmnet,
		eply,
		stats,
		mgcv,
		party,
		partykit
		tidyposterior
		rsample
		
```

In addition, it's advisable to also install Stan, Installation and documentation can be found here: http://mc-stan.org/users/documentation/
### Installing
Before installing, make sure that you have all the dependencies shown at the previous section.
```
devtools::install_github("csetraynor/predsurv")
```
On the other hand it is also possible to use dependencies = TRUE.

## Acknowledgments
R Core Team (2017). R: A language and environment for statistical
  computing. R Foundation for Statistical Computing, Vienna, Austria.
  URL https://www.R-project.org/.

Noah Simon, Jerome Friedman, Trevor Hastie, Rob Tibshirani (2011).
  Regularization Paths for Cox's Proportional Hazards Model via
  Coordinate Descent. Journal of Statistical Software, 39(5), 1-13. URL
  http://www.jstatsoft.org/v39/i05/.

Ulla B. Mogensen, Hemant Ishwaran, Thomas A. Gerds (2012). Evaluating
  Random Forests for Survival Analysis Using Prediction Error Curves.
  Journal of Statistical Software, 50(11), 1-23. URL
  http://www.jstatsoft.org/v50/i11/.

Max Kuhn and Hadley Wickham (2017). rsample: General Resampling
  Infrastructure. R package version 0.0.2.
  https://CRAN.R-project.org/package=rsample
  
Stan Development Team (2018). RStan: the R interface to Stan. R
  package version 2.17.3. http://mc-stan.org/.




