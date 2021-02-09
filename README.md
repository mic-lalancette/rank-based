# rank-based
This repository contains the R code necessary to reproduce the numerical experiments found in the paper "Rank-based estimation under asymptotic dependence and independence, with applications to spatial extremes" by M. Lalancette, S. Engelke and S. Volgushev.

First run `BivEstimators.R`, `WfEstimators.R`, `SpEstimators.R` and `DataEstimators.R` to create different `.Rdata` files. Then running `PaperPlots.R` will retrieve the estimators in those data files to sequentially create all the plots and save them (in the `.eps` format) in the `Plots` folder.

Note that running `DataEstimators.R` necessitates the rainfall data set `RainData.Rdata` (currently not included here for confidentiality). The corresponding section in `PaperPlots.R` is commented out.
