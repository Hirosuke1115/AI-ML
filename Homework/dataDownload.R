# =============================================================================
# LEGO_noise_study.R
# Does reducing noise in the target change what we learn about predictability?
#   Treatment A: average prices over k months at each end of the test window
#                (k = 1 is the current target, k = 3 and k = 5 are smoothed,
#                each compared with a single-month target of the same horizon)
#   Treatment B: cross-validation folds grouped by theme (no theme in both
#                training and test), compared with random folds by set
#   Treatment C: score the forecasts on theme averages as well as single sets
# Data: Dobrynskaya and Kishilova, Mendeley Data v9hhs66vm3 (same files as
#       xLEGO_Project.Rmd). Put this file next to the Excel files and run it.
# =============================================================================

# ---- 0. Packages and settings -------------------------------------------------
library(tidyverse)                                   # Data handling and plots
library(readxl)                                      # Reads the Excel files
library(glmnet)                                      # Ridge, LASSO, elastic net
library(randomForest)                                # Random forest

data_dir <- "."                                      # Folder with the DiB Excel files
out_dir  <- "noise_study_output"                     # Results are written here
dir.create(out_dir, showWarnings = FALSE)            # Create the folder if it is missing
N_FOLDS   <- 5                                       # Number of cross-validation folds
SEED      <- 42                                      # One seed for every random step

# ---- 1. December 2015 cross-section (features), as in xLEGO_Project.Rmd -------
raw_2015 <- read_excel(file.path(data_dir, "DiB whole sample Dec 2015 prices.xlsx"),
                       sheet = "DATA", na = c("", "NA", "ND"))           # One row per set
names(raw_2015) <- c("id", "theme", "name", "year", "pieces", "minifigs",
                     "p_new_2015", "p_used_2015", "p_retail", "age", "size")   # Short column names
sets <- raw_2015 %>%
  mutate(across(c(minifigs, p_new_2015, p_used_2015), as.numeric),       # Text columns to numbers
         id    = as.character(id),                                       # Set number as text key
         theme = factor(theme),                                          # Theme as a category
         size  = factor(size, levels = 1:4),                             # Size group 1 = biggest
         log_retail = log(p_retail)) %>%                                 # Retail price in logs
  filter(p_new_2015 > 0, p_retail > 0, !is.na(minifigs)) %>%             # Same filter as the Rmd
  droplevels() %>%                                                       # Remove empty theme levels
  mutate(y15     = log(p_new_2015 / p_retail) / (2015 - year),           # Return since release, p.a.
         prem15  = log(p_new_2015 / p_retail),                           # Premium over retail, Dec 2015
         log_p15 = log(p_new_2015),                                      # Price level, Dec 2015
         theme_ml  = fct_lump_min(theme, min = 15, other_level = "Other"),        # Pool rare themes
         year_ml   = factor(ifelse(year < 1995, "pre-1995", as.character(year))), # Pool early years
         pieces_sq = pieces^2)                                           # Allows a curved size effect

# ---- 2. Monthly prices January 2018 to April 2019 ----------------------------
raw_1819 <- read_excel(file.path(data_dir, "DiB whole sample Jan 2018 - Apr 2019 prices.xlsx"),
                       skip = 1, na = c("", "NA", "ND"))                 # Row 1 is a title
month_cols <- names(raw_1819)[-(1:6)]                                    # Month headers (Excel serials)
names(raw_1819)[-(1:6)] <- as.character(as.Date(as.numeric(month_cols),
                                                origin = "1899-12-30"))  # Serials to dates
prices <- raw_1819 %>%
  mutate(id = as.character(id)) %>%                                      # Same key as 'sets'
  select(id, matches("^20")) %>%                                         # id plus the 16 month columns
  pivot_longer(-id, names_to = "date", values_to = "p") %>%              # Long: one row per set-month
  mutate(date = as.Date(date), p = as.numeric(p)) %>%                    # Proper types
  filter(!is.na(p), p > 0) %>%                                           # Keep valid prices only
  mutate(m = (year(date) - 2018) * 12 + month(date))                     # Month index: Jan 2018 = 1 ... Apr 2019 = 16

# ---- 3. Noise diagnostic: Roll (1984) autocovariance of monthly changes -------
# Pure noise in a price makes this month's change and next month's change move
# in opposite directions, so their covariance is negative. Roll (1984): the
# noise standard deviation is about sqrt(-cov) when cov < 0.
roll <- prices %>%
  arrange(id, m) %>%                                                     # Months in order within each set
  group_by(id) %>%
  mutate(dp = log(p) - log(lag(p)),                                      # Monthly log change
         consecutive = m - lag(m) == 1) %>%                              # TRUE if the previous month exists
  mutate(dp_lag = lag(dp), ok = consecutive & lag(consecutive)) %>%      # Previous change, both months adjacent
  ungroup() %>%
  filter(ok %in% TRUE)                                                   # Pairs of adjacent changes only
roll_tab <- tibble(
  pairs           = nrow(roll),                                          # Number of change pairs used
  autocorrelation = cor(roll$dp, roll$dp_lag),                           # Negative = noise
  autocovariance  = cov(roll$dp, roll$dp_lag),
  noise_sd_roll   = ifelse(autocovariance < 0, sqrt(-autocovariance), NA),   # Implied noise SD (log units)
  sd_monthly_change = sd(roll$dp))                                       # Total SD of a monthly change
print(roll_tab)                                                          # Show on screen
write_csv(roll_tab, file.path(out_dir, "1_roll_noise_diagnostic.csv"))   # Save

# ---- 4. Targets: return over the window, with k-month averages at each end ---
# Smoothed target: start window = months 1..k, end window = months (17-k)..16.
#   Each end is the mean log price of the k months.
# Single target with the SAME horizon: one month at the centre of each window,
#   i.e. month (k+1)/2 and month 16-(k-1)/2. k must be odd so centres are whole months.
# Comparing the two at the same k isolates the effect of averaging, because
# wider averaging also shortens the gap between the two ends.
make_target <- function(k, smooth){
  if (smooth) { m_start <- 1:k;           m_end <- (17 - k):16 }         # k months at each end
  else        { m_start <- (k + 1) / 2;   m_end <- 16 - (k - 1) / 2 }    # One centre month at each end
  start <- prices %>% filter(m %in% m_start) %>%                         # Months in the start window
    group_by(id) %>% summarise(lp0 = mean(log(p)), n0 = n(), .groups = "drop")   # Mean log price
  end <- prices %>% filter(m %in% m_end) %>%                             # Months in the end window
    group_by(id) %>% summarise(lp1 = mean(log(p)), n1 = n(), .groups = "drop")
  h <- (mean(m_end) - mean(m_start)) / 12                                # Holding period in years (same for both)
  inner_join(start, end, by = "id") %>%                                  # Sets with both ends
    filter(n0 == length(m_start), n1 == length(m_end)) %>%               # Every month present, no partial averages
    transmute(id, y = (lp1 - lp0) / h,                                   # Annualised log return
              k = k, target = ifelse(smooth, "Averaged", "Single month"))
}
K_WINDOWS <- c(1, 3, 5)                                                  # Odd k only (see above)
targets <- bind_rows(make_target(1, smooth = FALSE),                     # k = 1: Jan 2018 -> Apr 2019
                     map_dfr(K_WINDOWS[-1], make_target, smooth = TRUE), # Averaged, k = 3 and 5
                     map_dfr(K_WINDOWS[-1], make_target, smooth = FALSE))# Single month, same horizons

# Same sets for every target, so differences are not driven by the sample
n_targets  <- n_distinct(paste(targets$k, targets$target))               # Number of target versions (5)
common_ids <- targets %>% count(id) %>% filter(n == n_targets) %>% pull(id)
targets    <- targets %>% filter(id %in% common_ids)
target_sd  <- targets %>% group_by(k, target) %>%
  summarise(sets = n(), mean = mean(y), sd = sd(y), .groups = "drop")    # Averaging should lower the SD at equal k
print(target_sd)
write_csv(target_sd, file.path(out_dir, "2_target_summary.csv"))

# Agreement between target versions (rank correlations)
target_cor <- targets %>% mutate(v = paste0("k", k, "_", target)) %>%
  select(id, v, y) %>% pivot_wider(names_from = v, values_from = y) %>%
  select(-id) %>% cor(method = "spearman")
print(round(target_cor, 3))
write.csv(round(target_cor, 3), file.path(out_dir, "2b_target_rank_correlations.csv"))

# ---- 5. Features (information as of December 2015) ---------------------------
f_char <- ~ theme_ml + size + year_ml + pieces + pieces_sq + minifigs + log_retail      # Characteristics
f_mom  <- ~ theme_ml + size + year_ml + pieces + pieces_sq + minifigs + log_retail +
  prem15 + log_p15 + y15                                       # + 2015 price history
rhs_char <- c("theme_ml", "size", "year_ml", "pieces", "minifigs", "log_retail")   # For the forest
rhs_mom  <- c(rhs_char, "prem15", "log_p15", "y15")
FEATURES <- list("Characteristics" = list(f = f_char, rhs = rhs_char),
                 "Characteristics + 2015 price history" = list(f = f_mom, rhs = rhs_mom))
MODELS <- c("Naive (mean)", "OLS", "LASSO", "Elastic net", "Random forest")

# ---- 6. Folds: random by set, or grouped by theme ----------------------------
make_folds <- function(df, scheme, k = N_FOLDS, seed = SEED){
  set.seed(seed)                                                         # Same folds every run
  if (scheme == "Random by set") {
    sample(rep(1:k, length.out = nrow(df)))                              # Deal sets into k piles
  } else {
    th <- unique(as.character(df$theme))                                 # Original themes, not lumped
    fold_th <- sample(rep(1:k, length.out = length(th)))                 # Deal whole themes into k piles
    fold_th[match(as.character(df$theme), th)]                           # Each set gets its theme's pile
  }
}

# ---- 7. Fit one model on training rows, predict test rows --------------------
fit_predict <- function(model, X_tr, y_tr, X_te, df_tr, df_te, rhs, seed){
  set.seed(seed)                                                         # Same seed for every model
  switch(model,
         "Naive (mean)" = rep(mean(y_tr), nrow(X_te)),                        # Training mean for everyone
         "OLS" = { b <- coef(lm.fit(cbind(1, X_tr), y_tr))                    # Least squares
         b[is.na(b)] <- 0                                           # Columns absent in training: weight 0
         as.numeric(cbind(1, X_te) %*% b) },
         "LASSO"       = as.numeric(predict(cv.glmnet(X_tr, y_tr, alpha = 1),   X_te, s = "lambda.min")),
         "Elastic net" = as.numeric(predict(cv.glmnet(X_tr, y_tr, alpha = 0.5), X_te, s = "lambda.min")),
         "Random forest" = predict(randomForest(reformulate(rhs, "y"), data = df_tr,
                                                ntree = 300, mtry = max(2, floor(length(rhs) / 3)),
                                                nodesize = 10), df_te))       # Same settings as the Rmd
}

# ---- 8. Run every combination and keep the out-of-fold predictions ----------
grid <- expand_grid(distinct(targets, k, target),                        # 5 target versions
                    scheme = c("Random by set", "Grouped by theme"),
                    features = names(FEATURES))                          # 5 x 2 x 2 = 20 designs
oof <- pmap_dfr(grid, function(k, target, scheme, features){
  df <- sets %>% inner_join(targets %>% filter(k == !!k, target == !!target), by = "id")   # Features + this target
  X  <- model.matrix(FEATURES[[features]]$f, data = df)[, -1, drop = FALSE]   # Dummies, no intercept
  fold <- make_folds(df, scheme)                                         # Fold label for each set
  map_dfr(MODELS, function(mdl){
    pred <- rep(NA_real_, nrow(df))                                      # Filled fold by fold
    for (f in 1:N_FOLDS) {
      tr <- fold != f                                                    # Training rows
      pred[!tr] <- fit_predict(mdl, X[tr, , drop = FALSE], df$y[tr], X[!tr, , drop = FALSE],
                               df[tr, ], df[!tr, ], FEATURES[[features]]$rhs, seed = SEED + f)
    }
    tibble(k = k, target = target, scheme = scheme, features = features, model = mdl,
           id = df$id, theme = as.character(df$theme), y = df$y, pred = pred)
  })
})
write_csv(oof, file.path(out_dir, "3_out_of_fold_predictions.csv"))      # Keep for later checks

# ---- 9. Scores on single sets ------------------------------------------------
score_sets <- oof %>%
  group_by(k, target, scheme, features, model) %>%
  summarise(sets       = n(),
            oos_r2     = 1 - mean((y - pred)^2) / mean((y - mean(y))^2),   # Out-of-sample R squared
            rmse       = sqrt(mean((y - pred)^2)),
            rmse_to_sd = rmse / sd(y),                                    # Comparable across k
            rank_corr  = if (sd(pred) > 0) cor(pred, y, method = "spearman") else NA_real_,
            median_hit = if (sd(pred) > 0) mean((pred > median(pred)) == (y > median(y))) else NA_real_,
            .groups = "drop") %>%
  arrange(features, scheme, model, k, target)
print(score_sets, n = Inf, width = Inf)
write_csv(score_sets, file.path(out_dir, "4_scores_single_sets.csv"))

# ---- 10. Scores on theme averages --------------------------------------------
# Average the realised and the predicted return within each theme, then score
# the themes. Set-level noise partly cancels inside a theme average.
MIN_SETS_PER_THEME <- 10                                                 # Themes with fewer sets are too noisy
score_themes <- oof %>%
  group_by(k, target, scheme, features, model, theme) %>%
  summarise(n = n(), y = mean(y), pred = mean(pred), .groups = "drop") %>%   # One row per theme
  filter(n >= MIN_SETS_PER_THEME) %>%
  group_by(k, target, scheme, features, model) %>%
  summarise(themes    = n(),
            oos_r2    = 1 - sum(n * (y - pred)^2) / sum(n * (y - weighted.mean(y, n))^2),   # Weighted by set count
            rank_corr = if (sd(pred) > 0) cor(pred, y, method = "spearman") else NA_real_,
            .groups = "drop") %>%
  arrange(features, scheme, model, k, target)
print(score_themes, n = Inf, width = Inf)
write_csv(score_themes, file.path(out_dir, "5_scores_theme_averages.csv"))

# ---- 11. One summary chart ---------------------------------------------------
p <- score_sets %>%
  filter(model != "Naive (mean)") %>%
  ggplot(aes(x = factor(k), y = rank_corr, colour = model, shape = target,
             group = interaction(model, target))) +                      # One line per model and target type
  geom_line() + geom_point(size = 2) +
  facet_grid(features ~ scheme) +                                        # Rows: features, columns: folds
  labs(x = "Window size k (averaged months, or centre month of that window)",
       y = "Rank correlation, held-out sets", colour = NULL, shape = NULL) +
  theme_light() + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "6_rank_corr_by_k.png"), p, width = 9, height = 6, dpi = 150)

# ---- 12. Session info for reproducibility -----------------------------------
writeLines(capture.output(sessionInfo()), file.path(out_dir, "7_session_info.txt"))