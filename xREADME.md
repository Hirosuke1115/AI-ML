# LEGO as an investment: does the 2015 evidence predict out of sample?

Empirical project, AI and Machine Learning (EX), WS 26/27. Hiro, Hendrik, Stefan.

## How to run
1. Put `LEGO_Project.Rmd` in one folder with the four data files from the paper's Mendeley dataset
   (Dobrynskaya and Kishilova, v9hhs66vm3): `DiB whole sample Dec 2015 prices.xlsx` and
   `DiB whole sample Jan 2018 - Apr 2019 prices.xlsx` are used; the other two are included for reference.
2. Packages: `tidyverse`, `readxl`, `glmnet`, `rpart`, `rpart.plot`, `randomForest`, `sandwich`, `knitr`
   (installed automatically if missing). `xgboost` is optional: the XGBoost rows appear only if it is installed.
3. Open the Rmd in RStudio and click Knit. Run time is about 3 minutes (the random forests and the cross-validation).
   `LEGO_Project_code_only.R` is the same code as a plain script.

## Files
- `LEGO_Project.Rmd` / `LEGO_Project.html`: the project (replication, Experiment 1, Experiment 2, clean window, sorts, persistence check, limitations).
- `LegoProject.R`: Hiro's first version (replication and the 2018-2019 comparison); its replication block is reused in the Rmd.
