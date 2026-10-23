library(tidyverse)                                    # dplyr, tidyr, ggplot2
library(readxl)                                       # Read .xlsx files
library(glmnet)                                       # Ridge, LASSO, elastic net
library(rpart)                                        # Trees
library(rpart.plot)                                   # Tree plots

data_dir <- "."
test_month <- as.Date("2018-12-01")
T_train <- 2015
T_test <- as.numeric(format(test_month, "%Y"))

# 0. Data
raw_2015 <- read_excel(file.path(data_dir, "DiB whole sample Dec 2015 prices.xlsx"),
                       sheet = "DATA",
                       na = c("", "NA", "ND"))
names(raw_2015) <- c("id", "theme", "name", "year", "pieces", "minifigs", "p_new_2015", "p_used_2015", "p_retail", "age", "size")
sets <- raw_2015 %>% 
  mutate(across(c(minifigs, p_new_2015, p_used_2015), as.numeric),
         id = as.character(id),
         theme = relevel(factor(theme),ref = "Miscellaneous"),
         size = factor(size, levels = 1:4),
         year_f = factor(year),
         log_retail = log(p_retail)) %>% 
           filter(p_new_2015 > 0, p_retail > 0, !is.na(minifigs)) %>% # exclude $ = 0: no sale
           droplevels() # drop empty factor levels

# Target: yearly log return since release
training_sample <- sets %>% 
  mutate(h = T_train - year,
         y = log(p_new_2015 / p_retail) / h) # yearly log return

# Test file set up
raw_test <- read_excel(file.path(data_dir, "DiB whole sample Jan 2018 - Apr 2019 prices.xlsx"),
                       skip = 1,# Raw 1 is a title
                       na = c("", "NA", "ND")) # "ND" is no data
month_cols <- names(raw_test)[-(1:6)] # Month headers are Excel serial numbers
names(raw_test)[-(1:6)] <- as.character(as.Date(as.numeric(month_cols),
                                                origin = "1899-12-30")) # Convert to dates

prices_test <- raw_test %>% 
  mutate(id = as.character(id)) %>% 
  select(id, matches("^20")) %>% # keep the id and all 16 month columns
  pivot_longer(-id, names_to = "date", values_to = "p_test") %>%  # long: one row per set x month
  mutate(date = as.Date(date))

testing_sample <- sets %>% 
  inner_join(prices_test, by = "id") %>% 
  filter(!is.na(p_test), p_test > 0) %>% 
  mutate(h = (year(date) - year) + (month(date) - 12) / 12, # Years held (release = end of release year)
         y = log(p_test / p_retail) / h)

y_train <- training_sample$y
y_test <- testing_sample$y
stopifnot(length(y_test) == nrow(testing_sample))

# 1. OLS test
# 1.a Replicate the table 5
fit_paper <- lm(log(p_new_2015 / p_retail) ~ theme + size + year_f,
                data = training_sample)
print(summary(fit_paper))
print(round(c(coef(fit_paper)[c("themeBatman", "themeDiscovery", "size2")],
              R2 = summary(fit_paper)$r.squared), 4)) # the result is the same as in the paper: 0.4484, 0.6124, -0.2438, R2 0.4695

# 1.b OLS for forecasting
f_paper <- y ~ theme + size + year_f
f_full <- y ~ theme + size + year_f + pieces + I(pieces^2) + minifigs + log_retail

fit_ols_paper <- lm(f_paper, data = training_sample)
fit_ols_full <- lm(f_full, data = training_sample)

pred_ols_paper <- predict(fit_ols_paper, newdata = testing_sample)
pred_ols_full <- predict(fit_ols_full, newdata = testing_sample)
pred_naive <- rep(mean(y_train), nrow(testing_sample))

# 2. Ridge, Lasso, Elasticnet
x_formula <- ~ theme + size + year_f + pieces + I(pieces^2) + minifigs + log_retail
x_penalized <- model.matrix(x_formula, data = training_sample)[, -1] # Dummies; drop intercept
x_test <- model.matrix(x_formula, data = testing_sample)[, -1]
y_penalized <- y_train

# 2.a Elasticnet
alpha_seq <- (1:10) / 10
lambda_seq <- 10^seq(0, -5, by = -0.5)
enet_res <- lapply(alpha_seq, function(a){
  fit_enet <- glmnet(x_penalized, y_penalized, alpha = a, lambda = lambda_seq)
  pred <- predict(fit_enet, x_test)
  data.frame(alpha = a, 
             lambda = fit_enet$lambda,
             rmse = sqrt(colMeans((pred - y_test)^2)),
             nb_features = fit_enet$df) # Non-zero coef
}) %>% bind_rows()
ggplot(enet_res, aes(x = factor(alpha), y = factor(signif(lambda, 2)), fill = rmse)) +
  geom_tile() + 
  geom_text(aes(label = sprintf("%.4f", rmse)), size = 2) +
  scale_fill_gradient(low = "#0570EA", high = "white") + 
  labs(x = "alpha", y = "lambda", fill = "Test RMSE") +
  theme_light()

print(enet_res %>% arrange(rmse) %>% head(5))

# 2.b. Tuning
set.seed(42) # Reproducible folds
fold_id <- sample(rep(1:10, length.out = nrow(x_penalized))) # Same folds for every alpha

cv_fit <- function(a){
  cv <- cv.glmnet(x_penalized, y_penalized, alpha = a, foldid = fold_id) # cv = cross validation
  list(alpha = a, 
       lambda = cv$lambda.min,
       cv_err = min(cv$cvm),
       pred = as.numeric(predict(cv, x_test, s = "lambda.min")),
       fitted = as.numeric(predict(cv, x_penalized, s = "lambda.min")),
       nz = sum(coef(cv, s = "lambda.min")[-1] != 0))
}

ridge <- cv_fit(0)
lasso <- cv_fit(1)
enet_cv <- lapply(alpha_seq[2:10], cv_fit)
enet <- enet_cv[[which.min(sapply(enet_cv, `[[`, "cv_err"))]] # Best alpha by cv error

# 3. Regression tree
f_tree <- y ~ theme + size + year_f + pieces + minifigs + log_retail
fit_tree <- rpart(f_tree, data = training_sample,
                  minbucket = 50,
                  minsplit = 100,
                  cp = 0.001,
                  maxdepth = 4)
rpart.plot(fit_tree)
pred_tree <- predict(fit_tree, testing_sample)

# 4. Classification tree
training_sample <- training_sample %>% 
  mutate(y_C = as.factor(y > median(y))) #Above 2015 median?
testing_sample <- testing_sample %>% 
  group_by(date) %>% 
  mutate(y_C = as.factor(y > median(y))) %>% # Above that month median?
  ungroup()

fit_tree_C <- rpart(update(f_tree, y_C ~ .), data = training_sample, method = "class",
                    minbucket = 50,
                    minsplit = 100,
                    cp = 0.001, 
                    maxdepth = 4)
rpart.plot(fit_tree_C)
pred_C <- predict(fit_tree_C, testing_sample, type = "class") # Majority label of the leaf
print(table(Predicted = pred_C, Actual = testing_sample$y_C)) # Confusion Matrix

# 5. Comparison over all 16 test months
preds <- list(`OLS (paper features)` = pred_ols_paper,
              `OLS (all features)`   = pred_ols_full,
              Ridge                  = ridge$pred,
              LASSO                  = lasso$pred,
              `Elastic net`          = enet$pred,
              `Regression tree`      = pred_tree,
              `Naive (mean)`         = pred_naive)

# Long table: one row per model x set x month
long <- testing_sample %>%
  select(id, date, y) %>%                             # Keys and the actual return
  bind_cols(as_tibble(preds)) %>%                     # One column of predictions per model
  pivot_longer(-c(id, date, y), names_to = "Model", values_to = "pred") %>%  # Stack the models
  group_by(Model, date) %>%                           # Medians are taken within each month
  mutate(sign_ok = pred * y > 0,                      # Book hit: same sign?
         med_ok  = if (sd(pred) == 0) NA              # Undefined for a constant forecast
         else (pred > median(pred)) == (y > median(y))) %>%  # Same side of the month's median?
  ungroup()

long <- long %>%
  bind_rows(testing_sample %>%                        # Add the classification tree
              transmute(id, date, y, Model = "Classification tree",
                        pred = NA_real_, sign_ok = NA,
                        med_ok = pred_C == y_C))      # Correct class = median hit

# 5a. Pooled results (all set x month rows)
fitted_train <- tibble(Model = names(preds)[1:6],
                       `Train MSE` = c(mean((fitted(fit_ols_paper) - y_train)^2),
                                       mean((fitted(fit_ols_full)  - y_train)^2),
                                       mean((ridge$fitted - y_train)^2),
                                       mean((lasso$fitted - y_train)^2),
                                       mean((enet$fitted  - y_train)^2),
                                       mean((predict(fit_tree, training_sample) - y_train)^2)))

pooled <- long %>%
  group_by(Model) %>%
  summarise(Rows         = n(),                       # Set x month observations
            `Test MSE`   = mean((pred - y)^2),        # NA for the classification tree
            `Sign hit`   = mean(sign_ok),
            `Median hit` = mean(med_ok)) %>%
  left_join(fitted_train, by = "Model") %>%
  arrange(`Test MSE`)
print(pooled, width = Inf)
cat("Sign-hit baseline (share of positive test returns):", round(mean(y_test > 0), 3), "\n")

# 5b. Month by month
monthly <- long %>%
  group_by(Model, date) %>%
  summarise(MSE = mean((pred - y)^2), `Median hit` = mean(med_ok), .groups = "drop")

monthly %>%
  filter(!is.na(`Median hit`)) %>%
  ggplot(aes(x = date, y = `Median hit`, color = Model)) +
  geom_line() + geom_point(size = 1) +
  geom_hline(yintercept = 0.5, linetype = 2) +        # Coin-flip baseline
  labs(x = NULL, y = "Median hit ratio") + theme_light()

# 5c. Is the median hit ratio significant? Three ways, from naive to honest.
#   (1) Binomial test on all rows: treats every set x month as independent. Overstated,
#       because each set appears up to 16 times with nearly the same price.
#   (2) t-test across the 16 monthly hit ratios: months are not independent either.
#   (3) t-test across sets (cluster by set): each set's share of correct months is one
#       observation. This is the honest test; extra months add little new information.
sig <- long %>%
  filter(!is.na(med_ok)) %>%
  group_by(Model) %>%
  group_modify(~{
    by_month <- .x %>% group_by(date) %>% summarise(h = mean(med_ok))   # One hit ratio per month
    by_set   <- .x %>% group_by(id)   %>% summarise(h = mean(med_ok))   # One hit ratio per set
    tibble(`Hit`             = mean(.x$med_ok),
           `p binomial (rows)` = binom.test(sum(.x$med_ok), nrow(.x), p = 0.5)$p.value,
           `Worst month`     = min(by_month$h),
           `Best month`      = max(by_month$h),
           `t (months)`      = unname(t.test(by_month$h, mu = 0.5)$statistic),
           `t (sets)`        = unname(t.test(by_set$h,   mu = 0.5)$statistic),
           `p (sets)`        = t.test(by_set$h, mu = 0.5)$p.value,
           `Same answer all months` = mean(by_set$h %in% c(0, 1)))  # Share of sets never flipping
  }) %>%
  ungroup() %>%
  arrange(desc(Hit))
print(sig, width = Inf)

# 5d. Does the sign hit beat "always predict up"? Clustered by set.
sign_vs_up <- long %>%
  filter(!is.na(sign_ok)) %>%
  group_by(Model, id) %>%
  summarise(edge = mean(sign_ok) - mean(y > 0), .groups = "drop") %>%  # Model hit minus always-up hit
  group_by(Model) %>%
  summarise(`Edge vs always-up` = mean(edge),
            `p (sets)` = if (sd(edge) == 0) NA else t.test(edge)$p.value)