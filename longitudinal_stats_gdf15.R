# ── Packages ──────────────────────────────────────────────────────────────────
library(lme4); library(lmerTest); library(tidyverse); library(emmeans)
library(ggplot2); library(readxl); library(performance); library(nortest)
library(nlme); library(patchwork); library(forcats); library(pROC)

setwd("~/IMRB/Stats GDF15")

# ══════════════════════════════════════════════════════════════════════════════
# 1. IMPORT & CLEANING
# ══════════════════════════════════════════════════════════════════════════════

recode_visit <- function(df, visit_col = "visit") {
  df |> mutate(time_num = case_when(
    str_detect(.data[[visit_col]], "V0") ~ 0,
    str_detect(.data[[visit_col]], "V1") ~ 1,
    str_detect(.data[[visit_col]], "V2") ~ 2,
    TRUE ~ NA_real_
  ))
}

df_long    <- read_excel("data/Complete_data_Evanne_GDF15_BD_cleaned.xlsx", sheet = "V0 V1 V2") |> recode_visit("visit")
df_ymrs    <- read_excel("data/Patients_Evanne_cleaned.xlsx")                                   |> recode_visit("visit")
df_samples <- read_excel("data/samples_V0V1V2_Evanne_GDF15.xlsx")                              |> recode_visit("VISITE")

df_long <- df_long |>
  merge(df_ymrs[, c("fondacode", "time_num", "ymrs_num", "edulevel")], by = c("fondacode", "time_num"))

df_filtered <- df_long |>
  merge(df_samples[, c("fondacode", "time_num", "GDF15 pg/ml")], by = c("fondacode", "time_num")) |>
  rename(GDF15 = `GDF15 pg/ml`)

cat("Fondacodes with no match in clinical data:\n")
df_samples |> anti_join(df_filtered, by = c("fondacode", "time_num")) |>
  distinct(fondacode, time_num) |> arrange(fondacode) |> print()

df_filtered <- df_filtered |>
  mutate(corrected_age = age + time_num) |>
  group_by(fondacode) |> filter(n() == 3) |> ungroup()

cat("\nObservations per patient (should all be 3):\n")
df_filtered |> count(fondacode) |> summary() |> print()

# Clean string columns
na_strings <- c("BD not specified", "Ne sais pas", "ne sais pas",
                "BD NOT SPECIFIED", "Unknown", "unknown", "NA", "N/A", "n/a")

df_filtered <- df_filtered |>
  mutate(
    across(where(~ is.character(.x) && any(grepl(",",    .x, fixed = TRUE))), ~ as.numeric(gsub(",", ".", .x))),
    across(where(~ is.character(.x) && any(grepl("^[<>]", .x))),              ~ as.numeric(ifelse(grepl("^[<>]", .x), NA, gsub(",", ".", .x)))),
    across(where(~ is.factor(.) | is.character(.)),                            ~ ifelse(trimws(.x) %in% na_strings, NA, .x))
  )

# Lithium treatment
lithium_data <- read_excel("data/final_merged_data_samples_V0.xlsx") |>
  select(fondacode, Lithium_treat) |>
  mutate(
    lithium_treat = case_when(
      Lithium_treat %in% c("Yes", "yes", "YES") ~ "Yes",
      Lithium_treat %in% c("No",  "no",  "NO")  ~ "No",
      TRUE ~ NA_character_
    ) |> factor(levels = c("No", "Yes"))
  ) |> select(fondacode, lithium_treat)

n_before <- n_distinct(df_filtered$fondacode)
df_filtered <- df_filtered |> left_join(lithium_data, by = "fondacode")
n_na_lithium <- sum(is.na(df_filtered$lithium_treat[!duplicated(df_filtered$fondacode)]))
message(sprintf("Merge done. %d unique patients. %d with NA lithium_treat.", n_before, n_na_lithium))

# ══════════════════════════════════════════════════════════════════════════════
# 2. OUTLIER REMOVAL & AGE CATEGORIES
# ══════════════════════════════════════════════════════════════════════════════

df_time0 <- df_filtered |> filter(time_num == 0)
num_vars  <- df_time0 |> select(where(is.numeric)) |> names()
df_patient <- df_time0 |> filter(fondacode == "0401616545")

cat("\nZ-scores for patient 0401616545 at T0:\n")
tibble(variable = num_vars,
       z = sapply(num_vars, function(v)
         (df_patient[[v]] - mean(df_time0[[v]], na.rm = TRUE)) / sd(df_time0[[v]], na.rm = TRUE))) |>
  arrange(desc(abs(z))) |> print()

df_filtered <- df_filtered |>
  filter(fondacode != "0401616545") |>
  mutate(age_cat = factor(case_when(
    corrected_age <  30 ~ "<30", corrected_age <= 40 ~ "30–40", TRUE ~ ">40"),
    levels = c("<30", "30–40", ">40")))

# ══════════════════════════════════════════════════════════════════════════════
# 3. BASE MODELS
# ══════════════════════════════════════════════════════════════════════════════

df_filtered$log_GDF15 <- log10(df_filtered$GDF15)
shapiro.test(df_filtered$GDF15); shapiro.test(df_filtered$log_GDF15)

# Random intercept only
m1    <- lmer(log_GDF15 ~ time_num + age + lithium_treat + (1 | fondacode),            data = df_filtered, REML = TRUE)
# Random intercept + slope
m2    <- lmer(log_GDF15 ~ time_num + age + lithium_treat + (1 + time_num | fondacode), data = df_filtered, REML = TRUE)
# Time as factor (no linearity assumption)
m3    <- lmer(log_GDF15 ~ visit   + age + lithium_treat + (1 + time_num | fondacode),  data = df_filtered, REML = TRUE)

m1_ml <- update(m1, REML = FALSE); m2_ml <- update(m2, REML = FALSE); m3_ml <- update(m3, REML = FALSE)

cat("\nm1 vs m2 (random slope):\n");   print(anova(m1_ml, m2_ml))
cat("\nm1 vs m3 (time linearity):\n"); print(anova(m1_ml, m3_ml))

m1 <- lmer(log_GDF15 ~ time_num + age + lithium_treat + (1 | fondacode), data = df_filtered, REML = TRUE)
m1_bis <- lmer(log_GDF15 ~ time_num * lithium_treat + time_num * age        + (1 | fondacode), data = df_filtered, REML = TRUE)
summary(m1); summary(m1_bis)

# ══════════════════════════════════════════════════════════════════════════════
# 4. MAIN VISUALISATION
# ══════════════════════════════════════════════════════════════════════════════
facebd_colors <- list(
  dark_burgundy   = "#6B1A3A", medium_rose     = "#A83060",
  light_rose      = "#D4789A", very_light_pink = "#F2B8CC",
  blue_gray       = "#8899AA", dark_navy       = "#2D3B6E"
)
pal_age  <- c("<30" = "#274690", "30–40" = "#E59CB7", ">40" = "#5A102D")
pal_lith <- c("No" = "#274690", "Yes" = "#5A102D", "Unknown" = "#E59CB7")

newdata_base <- tibble(
  time_num      = 0:2,
  age           = mean(df_filtered$age),
  lithium_treat = factor("No", levels = levels(df_filtered$lithium_treat))
) |> mutate(pred = predict(m1, newdata = pick(everything()), re.form = NA))

coefs       <- coef(summary(m1_bis))
fmt_p       <- function(p) if (p < 0.001) "p < 0.001" else if (p < 0.01) sprintf("p = %.3f **", p) else if (p < 0.05) sprintf("p = %.3f *", p) else sprintf("p = %.3f", p)
label_age   <- paste0("Age : ",     fmt_p(coefs["age",                       "Pr(>|t|)"]), "     Time × Age : ",     fmt_p(coefs["time_num:age",              "Pr(>|t|)"]))
label_lith  <- paste0("Lithium : ", fmt_p(coefs["lithium_treatYes",          "Pr(>|t|)"]), "     Time × Lithium : ", fmt_p(coefs["time_num:lithium_treatYes", "Pr(>|t|)"]))

df_lith <- df_filtered |> 
  mutate(lithium_plot = fct_na_value_to_level(lithium_treat, "Unknown"))

# Compute n
n_age <- df_filtered |>
  filter(!is.na(age_cat)) |>
  distinct(fondacode, age_cat) |>
  count(age_cat) |>
  deframe()

n_lith <- df_lith |>
  distinct(fondacode, lithium_plot) |>
  count(lithium_plot) |>
  deframe()

# Relabel palettes
pal_age_n <- setNames(
  pal_age,
  paste0(names(pal_age), " (n=", n_age[names(pal_age)], ")")
)

pal_lith_n <- setNames(
  pal_lith,
  paste0(names(pal_lith), " (n=", n_lith[names(pal_lith)], ")")
)

df_age <- df_filtered |>
  filter(!is.na(age_cat)) |>                      # keep only rows with known age_cat
  mutate(age_cat_n = factor(
    paste0(age_cat, " (n=", n_age[as.character(age_cat)], ")"),
    levels = names(pal_age_n)
  ))

df_lith <- df_lith |>
  mutate(lithium_plot_n = factor(
    paste0(lithium_plot, " (n=", n_lith[as.character(lithium_plot)], ")"),
    levels = names(pal_lith_n)
  ))

make_traj_panel <- function(data, colour_var, palette, legend_title, label_txt,
                            newdata = NULL, visit_lty = "dashed") {
  p <- ggplot() +
    geom_line(data = data, aes(time_num, log_GDF15, group = fondacode, colour = .data[[colour_var]]),
              alpha = 0.4, linewidth = 0.5) +
    geom_point(data = data, aes(time_num, log_GDF15, colour = .data[[colour_var]]), alpha = 0.4, size = 1.2) +
    stat_summary(data = data |> filter(!is.na(.data[[colour_var]])),
                 aes(time_num, log_GDF15, colour = .data[[colour_var]]),
                 fun = mean, geom = "line",  linewidth = 1.4, linetype = visit_lty) +
    stat_summary(data = data |> filter(!is.na(.data[[colour_var]])),
                 aes(time_num, log_GDF15, colour = .data[[colour_var]]),
                 fun = mean, geom = "point", size = 3, shape = 18) +
    annotate("text", x = 0, y = Inf, hjust = 0, vjust = 0.8,
             label = label_txt, size = 3.2, colour = "grey20", lineheight = 1.4) +
    scale_x_continuous(breaks = 0:2, labels = c("V0", "V1", "V2")) +
    scale_colour_manual(values = palette) +
    coord_cartesian(clip = "off") +
    labs(x = "Visit", y = "GDF15 [log10(pg/mL)]", colour = legend_title) +
    theme_classic(base_size = 14) +
    theme(plot.margin = margin(t = 30, r = 5, b = 5, l = 5))
  if (!is.null(newdata)) {
    p <- p +
      geom_line(data  = newdata, aes(time_num, pred, colour = "Population mean"), linewidth = 1.8) +
      geom_point(data = newdata, aes(time_num, pred, colour = "Population mean"), size = 3)
  }
  p
}

# ── Plots ────────────────────────────────────────────────────────────────────
p1 <- make_traj_panel(df_age,  "age_cat_n",
                      c(pal_age_n, "Population mean" = "grey30"),
                      "Age category", label_age, newdata = newdata_base) +
  labs(title = "By age category")

p2 <- make_traj_panel(df_lith, "lithium_plot_n",
                      pal_lith_n,
                      "Lithium treatment", label_lith) +
  labs(title = "By lithium treatment")

plot_age_lithium <- (p1 | p2) + plot_annotation(
  title = "Longitudinal trajectories of GDF15",
  theme = theme(plot.title = element_text(size = 16, face = "bold")))

print(plot_age_lithium)
ggsave("outcome/mixed_model_age_lithium_GDF15_numbers.png", plot_age_lithium, width = 12, height = 8, dpi = 300, bg = "white")

# ══════════════════════════════════════════════════════════════════════════════
# 5. HELPER FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

sig_stars <- function(p) case_when(
  is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", p < 0.10 ~ ".", TRUE ~ "")

pull_coef <- function(model, rows, coefs) {
  if (length(rows) == 0 || is.na(rows[1]))
    return(list(beta = NA_real_, se = NA_real_, tval = NA_real_, pval = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_))
  row   <- rows[1]
  beta  <- coefs[row, "Estimate"];    se   <- coefs[row, "Std. Error"]
  tval  <- coefs[row, "t value"]
  pval  <- if ("Pr(>|t|)" %in% colnames(coefs)) coefs[row, "Pr(>|t|)"] else
    if ("Pr(>|z|)" %in% colnames(coefs)) coefs[row, "Pr(>|z|)"] else NA_real_
  df_v  <- if ("df" %in% colnames(coefs)) coefs[row, "df"] else Inf
  list(beta = beta, se = se, tval = tval, pval = pval,
       ci_lo = beta - qt(0.975, df_v) * se, ci_hi = beta + qt(0.975, df_v) * se)
}

dw_panel <- function(model, id_col, time_col) {
  df_r    <- model@frame
  df_r$.r <- tryCatch(residuals(model, type = "pearson"), error = function(e) residuals(model))
  df_r    <- df_r[order(df_r[[id_col]], df_r[[time_col]]), ]
  sq_diffs <- unlist(lapply(unique(df_r[[id_col]]), function(p) {
    r <- df_r$.r[df_r[[id_col]] == p]; if (length(r) < 2) NULL else diff(r)^2 }))
  num <- sum(sq_diffs, na.rm = TRUE); den <- sum(df_r$.r^2, na.rm = TRUE)
  if (den == 0 || length(sq_diffs) == 0) return(c(DW = NA_real_, p = NA_real_))
  dw <- num / den; z <- (1 - dw / 2) * sqrt(nrow(df_r) - 1)
  c(DW = round(dw, 3), p = round(2 * (1 - pnorm(abs(z))), 4))
}

make_lhs <- function(v, tr) switch(tr,
                                   sqrt = paste0("sqrt(", v, ")"), log = paste0("log(", v, ")"),
                                   log10 = paste0("log10(", v, ")"), log1p = paste0("log1p(", v, ")"), v)

fit_lmer <- function(formula_str, data, fam_str = "gaussian") {
  f <- as.formula(formula_str)
  if (fam_str == "gaussian") {
    tryCatch(lmer(f, data = data, REML = TRUE), error = function(e) NULL,
             warning = function(w) suppressWarnings(lmer(f, data = data, REML = TRUE)))
  } else {
    fam_obj <- switch(fam_str, poisson = poisson("log"), binomial = binomial("logit"),
                      Gamma = Gamma("log"), gaussian("identity"))
    tryCatch(glmer(f, data = data, family = fam_obj), error = function(e) NULL,
             warning = function(w) suppressWarnings(glmer(f, data = data, family = fam_obj)))
  }
}

model_diagnostics <- function(model, id, time, outlier_thr = 3) {
  vc_df     <- tryCatch(as.data.frame(VarCorr(model)), error = function(e) NULL)
  resid_vec <- tryCatch(residuals(model, type = "pearson"), error = function(e) residuals(model))
  re_int    <- tryCatch(ranef(model)[[id]][, 1], error = function(e) NULL)
  dw        <- tryCatch(dw_panel(model, id_col = id, time_col = time), error = function(e) c(DW = NA_real_, p = NA_real_))
  std_r     <- resid_vec / sd(resid_vec, na.rm = TRUE)
  list(
    AIC             = round(AIC(model), 1), BIC = round(BIC(model), 1),
    ICC             = if (!is.null(vc_df) && nrow(vc_df) >= 2) round(vc_df$vcov[1] / sum(vc_df$vcov), 3) else NA_real_,
    shapiro_resid_p = round(tryCatch(shapiro.test(resid_vec)$p.value, error = function(e) NA_real_), 4),
    shapiro_re_p    = round(tryCatch(if (!is.null(re_int) && length(re_int) >= 3) shapiro.test(re_int)$p.value else NA_real_, error = function(e) NA_real_), 4),
    DW_stat = dw[["DW"]], DW_p = dw[["p"]],
    n_outliers = sum(abs(std_r) > outlier_thr, na.rm = TRUE),
    converged  = length(model@optinfo$conv$lme4$messages) == 0
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# 6. run_clinical_on_GDF15()
# clinical_var ~ time_num * log10(GDF15) + age + lithium + (1 | fondacode)
# ══════════════════════════════════════════════════════════════════════════════

run_clinical_on_GDF15 <- function(data, vars, gdf15 = "GDF15", log_gdf15 = TRUE,
                                  id = "fondacode", time = "time_num", age = "age",
                                  lithium = "lithium_treat", transformations = list(),
                                  families = list(), fdr_method = "BH", outlier_thr = 3) {
  gdf15_pred  <- if (log_gdf15) paste0("log10(", gdf15, ")") else gdf15
  covariates  <- if (!is.null(lithium) && nzchar(lithium)) paste(age, lithium, sep = " + ") else age
  
  results <- lapply(vars, function(v) {
    fam_str <- families[[v]] %||% "gaussian"
    tr      <- transformations[[v]] %||% "identity"
    formula_str <- sprintf("%s ~ %s * %s + %s + (1 | %s)", make_lhs(v, tr), time, gdf15_pred, covariates, id)
    model <- fit_lmer(formula_str, data, fam_str)
    if (is.null(model)) { warning(sprintf("Model failed: %s", v)); return(NULL) }
    
    coefs    <- summary(model)$coefficients
    int_rows <- grep(sprintf("^(%s:%s|%s:%s)$", time, gsub("([()])", "\\\\\\1", gdf15_pred),
                             gsub("([()])", "\\\\\\1", gdf15_pred), time), rownames(coefs), value = TRUE)
    main_rows   <- grep(paste0("^", gsub("([()])", "\\\\\\1", gdf15_pred), "$"), rownames(coefs), value = TRUE)
    lith_rows   <- grep(paste0("^", lithium), rownames(coefs), value = TRUE)
    int_coef    <- pull_coef(model, int_rows,  coefs)
    main_coef   <- pull_coef(model, main_rows, coefs)
    lithium_coef <- pull_coef(model, lith_rows, coefs)
    diag        <- model_diagnostics(model, id, time, outlier_thr)
    n_skip_li   <- if (!is.null(lithium) && nzchar(lithium)) length(unique(data[[id]][is.na(data[[lithium]])])) else NA_integer_
    
    tibble(
      variable = v, transform = tr, family = fam_str, n_skipped_lithium = n_skip_li,
      beta_slope = round(int_coef$beta, 4),  SE_slope = round(int_coef$se, 4),
      CI_low_slope = round(int_coef$ci_lo, 4), CI_high_slope = round(int_coef$ci_hi, 4),
      t_slope = round(int_coef$tval, 3),     p_slope = round(int_coef$pval, 4),
      beta_level = round(main_coef$beta, 4), SE_level = round(main_coef$se, 4),
      CI_low_level = round(main_coef$ci_lo, 4), CI_high_level = round(main_coef$ci_hi, 4),
      t_level = round(main_coef$tval, 3),    p_level = round(main_coef$pval, 4),
      beta_lithium = round(lithium_coef$beta, 4), SE_lithium = round(lithium_coef$se, 4),
      CI_low_lithium = round(lithium_coef$ci_lo, 4), CI_high_lithium = round(lithium_coef$ci_hi, 4),
      t_lithium = round(lithium_coef$tval, 3), p_lithium = round(lithium_coef$pval, 4),
      AIC = diag$AIC, BIC = diag$BIC, ICC = diag$ICC,
      shapiro_resid_p = diag$shapiro_resid_p, shapiro_re_p = diag$shapiro_re_p,
      DW_stat = diag$DW_stat, DW_p = diag$DW_p,
      n_outliers = diag$n_outliers, converged = diag$converged
    )
  }) |> setNames(vars)
  
  bind_rows(results) |>
    mutate(
      p_adj_fdr_slope   = round(p.adjust(p_slope,   method = fdr_method), 4),
      p_adj_fdr_level   = round(p.adjust(p_level,   method = fdr_method), 4),
      p_adj_fdr_lithium = round(p.adjust(p_lithium, method = fdr_method), 4),
      sig_slope_raw     = sig_stars(p_slope),   sig_slope_fdr   = sig_stars(p_adj_fdr_slope),
      sig_level_raw     = sig_stars(p_level),   sig_level_fdr   = sig_stars(p_adj_fdr_level),
      sig_lithium_raw   = sig_stars(p_lithium), sig_lithium_fdr = sig_stars(p_adj_fdr_lithium)
    ) |> arrange(p_slope)
}

# ══════════════════════════════════════════════════════════════════════════════
# 7. run_clinical_on_GDF15_within() — Mundlak decomposition
# ══════════════════════════════════════════════════════════════════════════════

run_clinical_on_GDF15_mundlak <- function(data, vars, gdf15 = "GDF15", log_gdf15 = TRUE,
                                         id = "fondacode", time = "time_num", age = "age",
                                         time_x_within = FALSE, adjust_lithium = FALSE,
                                         transformations = list(), families = list(),
                                         fdr_method = "BH", outlier_thr = 3) {
  if (adjust_lithium) {
    if (!"lithium_treat" %in% names(data)) stop("Column 'lithium_treat' not found.")
    n_na <- sum(is.na(data$lithium_treat[!duplicated(data[[id]])]))
    if (n_na > 0) message(sprintf("%d patient(s) with NA lithium_treat will be excluded.", n_na))
  }
  
  gdf15_raw <- if (log_gdf15) {
    log10_col         <- paste0("log10_", gdf15)
    data[[log10_col]] <- log10(data[[gdf15]])
    log10_col
  } else gdf15
  
  data <- data |>
    group_by(.data[[id]]) |>
    mutate(GDF15_mean = mean(.data[[gdf15_raw]], na.rm = TRUE),
           GDF15_dev  = .data[[gdf15_raw]] - GDF15_mean) |>
    ungroup()
  
  covariates   <- c(age, if (adjust_lithium) "lithium_treat")
  within_term  <- if (time_x_within) sprintf("GDF15_dev + %s:GDF15_dev", time) else "GDF15_dev"
  fixed_part   <- sprintf("%s + GDF15_mean + %s + %s", time, within_term, paste(covariates, collapse = " + "))
  
  results <- lapply(vars, function(v) {
    fam_str     <- families[[v]] %||% "gaussian"
    tr          <- transformations[[v]] %||% "identity"
    formula_str <- sprintf("%s ~ %s + (1 | %s)", make_lhs(v, tr), fixed_part, id)
    model <- fit_lmer(formula_str, data, fam_str)
    n_obs <- n_distinct(model@frame[[id]])
    if (is.null(model)) { warning(sprintf("Model failed: %s", v)); return(NULL) }
    
    coefs   <- summary(model)$coefficients
    within_rows  <- if (time_x_within) grep(paste0("^(", time, ":GDF15_dev|GDF15_dev:", time, ")$"), rownames(coefs), value = TRUE) else grep("^GDF15_dev$", rownames(coefs), value = TRUE)
    between_rows <- grep("^GDF15_mean$",                       rownames(coefs), value = TRUE)
    time_rows    <- grep(paste0("^", time, "$"),               rownames(coefs), value = TRUE)
    lith_rows    <- grep("^lithium_treatYes$",                 rownames(coefs), value = TRUE)
    
    within_coef  <- pull_coef(model, within_rows,  coefs)
    between_coef <- pull_coef(model, between_rows, coefs)
    time_coef    <- pull_coef(model, time_rows,    coefs)
    lithium_coef <- pull_coef(model, lith_rows,    coefs)
    diag         <- model_diagnostics(model, id, time, outlier_thr)
    
    tibble(
      variable = v, n = n_obs, transform = tr, family = fam_str, lithium_adjusted = adjust_lithium,
      beta_within  = round(within_coef$beta, 4),   SE_within  = round(within_coef$se, 4),
      CI_low_within  = round(within_coef$ci_lo, 4), CI_high_within  = round(within_coef$ci_hi, 4),
      t_within  = round(within_coef$tval, 3),       p_within  = round(within_coef$pval, 4),
      beta_between = round(between_coef$beta, 4),   SE_between = round(between_coef$se, 4),
      CI_low_between = round(between_coef$ci_lo, 4), CI_high_between = round(between_coef$ci_hi, 4),
      t_between = round(between_coef$tval, 3),       p_between = round(between_coef$pval, 4),
      beta_time = round(time_coef$beta, 4),          SE_time  = round(time_coef$se, 4),
      p_time    = round(time_coef$pval, 4),
      beta_lithium = round(lithium_coef$beta, 4),   SE_lithium = round(lithium_coef$se, 4),
      CI_low_lithium = round(lithium_coef$ci_lo, 4), CI_high_lithium = round(lithium_coef$ci_hi, 4),
      p_lithium = round(lithium_coef$pval, 4),
      AIC = diag$AIC, BIC = diag$BIC, ICC = diag$ICC,
      shapiro_resid_p = diag$shapiro_resid_p, shapiro_re_p = diag$shapiro_re_p,
      DW_stat = diag$DW_stat, DW_p = diag$DW_p,
      n_outliers = diag$n_outliers, converged = diag$converged
    )
  }) |> setNames(vars)
  
  bind_rows(results) |>
    mutate(
      p_adj_fdr_within  = round(p.adjust(p_within,  method = fdr_method), 4),
      p_adj_fdr_between = round(p.adjust(p_between, method = fdr_method), 4),
      sig_within_raw    = sig_stars(p_within),  sig_within_fdr  = sig_stars(p_adj_fdr_within),
      sig_between_raw   = sig_stars(p_between), sig_between_fdr = sig_stars(p_adj_fdr_between)
    ) |> arrange(p_within)
}

# ══════════════════════════════════════════════════════════════════════════════
# 8. plot_mixed_results() — trajectory plots by GDF15 quantile
# ══════════════════════════════════════════════════════════════════════════════

plot_mixed_results <- function(data, results, outcome = "GDF15", id = "fondacode",
                               time = "time_num", age = "age", lithium = "lithium_treat",
                               log_outcome = TRUE, filter_by = "p_slope",
                               plot_alpha = 0.05, ncol = 2) {
  if (!filter_by %in% names(results)) stop(sprintf("Column '%s' not found.", filter_by))
  sig_vars <- results |> filter(.data[[filter_by]] < plot_alpha) |> pull(variable)
  if (length(sig_vars) == 0) { message(sprintf("No variable with %s < %.2f.", filter_by, plot_alpha)); return(invisible(NULL)) }
  message(sprintf("Plotting %d variable(s): %s", length(sig_vars), paste(sig_vars, collapse = ", ")))
  
  gdf15_term   <- if (log_outcome) paste0("log10(", outcome, ")") else outcome
  int_term     <- paste0(time, ":", gdf15_term)
  gdf15_q25    <- quantile(data[[outcome]], 0.25, na.rm = TRUE)
  gdf15_q75    <- quantile(data[[outcome]], 0.75, na.rm = TRUE)
  age_m        <- mean(data[[age]], na.rm = TRUE)
  lab_q25      <- sprintf("Q25 (%.0f pg/mL)", gdf15_q25)
  lab_q75      <- sprintf("Q75 (%.0f pg/mL)", gdf15_q75)
  covariates   <- if (!is.null(lithium) && nzchar(lithium)) paste(age, lithium, sep = " + ") else age
  trend_cols   <- c(Q25 = facebd_colors$dark_burgundy, Q75 = facebd_colors$dark_navy)
  lithium_mode <- if (!is.null(lithium) && nzchar(lithium)) names(sort(table(data[[lithium]]), decreasing = TRUE))[1] else NULL
  sig_symbol   <- function(p) case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", p < 0.1 ~ ".", TRUE ~ "")
  phantom_df   <- data.frame(x = NA_real_, y = NA_real_)
  
  rev_stats <- list(); newdata_list <- list(); df_augs <- list(); n_patients <- list()  
  
  for (v in sig_vars) {
    res_row   <- results |> filter(variable == v)
    transform <- if ("transform" %in% names(res_row) && !is.na(res_row$transform) &&
                     !res_row$transform %in% c("identity", "none", "")) res_row$transform else ""
    v_lhs       <- if (nchar(transform) > 0) paste0(transform, "(", v, ")") else v
    formula_str <- sprintf("%s ~ %s * %s + %s + (1 | %s)", v_lhs, time, gdf15_term, covariates, id)
    model <- fit_lmer(formula_str, data)
    if (is.null(model)) { message(sprintf("Skipping %s — model failed.", v)); next }
    
    coef_tbl   <- summary(model)$coefficients
    alt_term   <- paste0(gdf15_term, ":", time)
    found_term <- if (int_term %in% rownames(coef_tbl)) int_term else if (alt_term %in% rownames(coef_tbl)) alt_term else NA_character_
    if (is.na(found_term)) {
      beta_rev <- p_rev <- ci_lo <- ci_hi <- NA_real_
    } else {
      beta_rev <- coef_tbl[found_term, "Estimate"]; se_rev <- coef_tbl[found_term, "Std. Error"]
      df_rev   <- coef_tbl[found_term, "df"];       p_rev  <- coef_tbl[found_term, "Pr(>|t|)"]
      ci_lo    <- beta_rev - qt(0.975, df_rev) * se_rev; ci_hi <- beta_rev + qt(0.975, df_rev) * se_rev
    }
    rev_stats[[v]] <- list(v_lhs = v_lhs, beta = beta_rev, p = p_rev, ci_lo = ci_lo, ci_hi = ci_hi)
    
    nd <- bind_rows(
      tibble(!!time := 0:2, !!age := age_m, !!outcome := gdf15_q25, group = lab_q25),
      tibble(!!time := 0:2, !!age := age_m, !!outcome := gdf15_q75, group = lab_q75)
    )
    if (!is.null(lithium_mode)) nd[[lithium]] <- lithium_mode
    nd$pred <- predict(model, newdata = nd, re.form = NA)
    
    df_augs[[v]]      <- model@frame |> mutate(y_obs = .data[[v_lhs]])
    n_patients[[v]]   <- n_distinct(model@frame[[id]])  
    newdata_list[[v]] <- nd
  }
  
  vars_ok <- names(rev_stats)
  p_fdr   <- setNames(p.adjust(sapply(vars_ok, function(v) rev_stats[[v]]$p), method = "fdr"), vars_ok)
  
  panel_list <- lapply(vars_ok, function(v) {
    s      <- rev_stats[[v]]; nd <- newdata_list[[v]]
    sub_txt <- sprintf("N = %d  \u00b7  \u03b2 = %.4f  \u00b7  p = %.4f%s  \u00b7  pFDR = %.4f%s",
                       n_patients[[v]], s$beta, s$p, sig_symbol(s$p), p_fdr[[v]], sig_symbol(p_fdr[[v]]))
    
    ggplot() +
      geom_line(data  = df_augs[[v]], aes(.data[[time]], y_obs, group = .data[[id]]), colour = "#8899AA", alpha = 0.12, linewidth = 0.4) +
      geom_point(data = df_augs[[v]], aes(.data[[time]], y_obs),                      colour = "#8899AA", alpha = 0.20, size = 0.9) +
      geom_line(data  = nd |> filter(group == lab_q25), aes(.data[[time]], pred, group = group), colour = trend_cols["Q25"], linetype = "dashed", linewidth = 0.9) +
      geom_point(data = nd |> filter(group == lab_q25), aes(.data[[time]], pred), colour = trend_cols["Q25"], shape = 1,  size = 3.5, stroke = 1.2) +
      geom_line(data  = nd |> filter(group == lab_q75), aes(.data[[time]], pred, group = group), colour = trend_cols["Q75"], linetype = "dashed", linewidth = 0.9) +
      geom_point(data = nd |> filter(group == lab_q75), aes(.data[[time]], pred), colour = trend_cols["Q75"], shape = 16, size = 3.5) +
      geom_line(data = phantom_df, aes(x, y, linetype = lab_q25), colour = trend_cols["Q25"], linewidth = 0.9) +
      geom_line(data = phantom_df, aes(x, y, linetype = lab_q75), colour = trend_cols["Q75"], linewidth = 0.9) +
      scale_x_continuous(breaks = 0:2, labels = c("V0", "V1", "V2")) +
      scale_linetype_manual(values = c("dashed", "dashed") |> setNames(c(lab_q25, lab_q75)), name = NULL) +
      guides(linetype = guide_legend(nrow = 1,
                                     override.aes = list(colour = unname(trend_cols), linewidth = c(1.1, 1.1)))) +
      labs(title = sprintf("Effect of GDF15 [log10(pg/mL)] on trajectory of %s", v),
           subtitle = sub_txt, x = "Visit", y = s$v_lhs) +
      theme_classic(base_size = 11) +
      theme(plot.title    = element_text(face = "bold", size = 12, colour = facebd_colors$dark_navy),
            plot.subtitle = element_text(size = 10, colour = "grey35", lineheight = 1.4),
            legend.key.width = unit(1.1, "cm"), legend.key.height = unit(0.5, "cm"))
  }) |> setNames(vars_ok)
  
  wrap_plots(panel_list, ncol = ncol) + plot_layout(guides = "collect") +
    plot_annotation(
      title   = "Clinical variable trajectories — mixed-model (GDF15 as predictor)",
      caption = "Adjusted for age and lithium treatment",
      theme   = theme(plot.title   = element_text(face = "bold", size = 13, colour = facebd_colors$dark_burgundy, hjust = 0.5),
                      plot.caption = element_text(size = 10, colour = "grey40", face = "italic", hjust = 0.5))) &
    theme(legend.position = "bottom", legend.box = "horizontal", legend.margin = margin(t = 4, b = 4),
          legend.text = element_text(size = 12), legend.title = element_text(size = 12, face = "bold"),
          legend.key.width = unit(1.6, "cm"), legend.key.height = unit(0.6, "cm"))
}


# ══════════════════════════════════════════════════════════════════════════════
# 9. plot_mixed_gdf15_xaxis() — GDF15 on x-axis, outcome on y-axis
# ══════════════════════════════════════════════════════════════════════════════

plot_mixed_gdf15_xaxis <- function(data, results, outcome = "GDF15", id = "fondacode",
                                   time = "time_num", age = "age", lithium = "lithium_treat",
                                   log_outcome = TRUE, age_cut = 35, filter_by = "p_slope",
                                   plot_alpha = 0.05, ncol = 2, show_obs = TRUE, n_pred = 100) {
  if (!filter_by %in% names(results)) stop(sprintf("Column '%s' not found.", filter_by))
  sig_vars <- results |> filter(.data[[filter_by]] < plot_alpha) |> pull(variable)
  if (length(sig_vars) == 0) { message(sprintf("No variable with %s < %.2f.", filter_by, plot_alpha)); return(invisible(NULL)) }
  
  gdf15_term   <- if (log_outcome) paste0("log10(", outcome, ")") else outcome
  int_term     <- paste0(time, ":", gdf15_term); alt_term <- paste0(gdf15_term, ":", time)
  time_vals    <- sort(unique(data[[time]])); visit_levels <- paste0("V", time_vals)
  age_m        <- mean(data[[age]], na.rm = TRUE)
  covariates   <- if (!is.null(lithium) && nzchar(lithium)) paste(age, lithium, sep = " + ") else age
  lithium_mode <- if (!is.null(lithium) && nzchar(lithium)) names(sort(table(data[[lithium]]), decreasing = TRUE))[1] else NULL
  visit_cols   <- setNames(c(facebd_colors$dark_navy, facebd_colors$dark_burgundy, "#D4789A")[seq_along(time_vals)], visit_levels)
  sig_symbol   <- function(p) case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", p < 0.1 ~ ".", TRUE ~ "")
  
  gdf15_raw_seq <- seq(quantile(data[[outcome]], 0.05, na.rm = TRUE),
                       quantile(data[[outcome]], 0.95, na.rm = TRUE), length.out = n_pred)
  x_pred_seq    <- if (log_outcome) log10(gdf15_raw_seq) else gdf15_raw_seq
  
  apply_tr <- function(x, tr) switch(tr, sqrt = sqrt(x), log = log(x), log10 = log10(x), log2 = log2(x), x)
  
  model_stats <- list(); obs_data <- list(); pred_data <- list(); n_patients <- list()
  
  for (v in sig_vars) {
    res_row   <- results |> filter(variable == v)
    transform <- if ("transform" %in% names(res_row) && !is.na(res_row$transform) && !res_row$transform %in% c("identity", "none", "")) res_row$transform else ""
    if (transform == "sqrt" && sum(data[[v]] < 0, na.rm = TRUE) > 0) { warning(sprintf("sqrt of negatives in '%s'. Skipping.", v)); next }
    v_lhs       <- if (nchar(transform) > 0) paste0(transform, "(", v, ")") else v
    formula_str <- sprintf("%s ~ %s * %s + %s + (1 | %s)", v_lhs, time, gdf15_term, covariates, id)
    model <- fit_lmer(formula_str, data)
    if (is.null(model)) { message(sprintf("Skipping '%s' — model failed.", v)); next }
    
    coef_tbl   <- summary(model)$coefficients
    found_term <- if (int_term %in% rownames(coef_tbl)) int_term else if (alt_term %in% rownames(coef_tbl)) alt_term else NA_character_
    if (is.na(found_term)) { beta <- se <- p_val <- ci_lo <- ci_hi <- NA_real_ } else {
      beta  <- coef_tbl[found_term, "Estimate"]; se   <- coef_tbl[found_term, "Std. Error"]
      df_t  <- coef_tbl[found_term, "df"];       p_val <- coef_tbl[found_term, "Pr(>|t|)"]
      ci_lo <- beta - qt(0.975, df_t) * se;      ci_hi <- beta + qt(0.975, df_t) * se
    }
    model_stats[[v]] <- list(transform = transform, v_lhs = v_lhs, beta = beta, p = p_val, ci_lo = ci_lo, ci_hi = ci_hi)
    
    obs_data[[v]] <- data |>
      filter(!is.na(.data[[v]]), !is.na(.data[[outcome]])) |>
      mutate(.y_obs = apply_tr(.data[[v]], transform),
             x_obs  = if (log_outcome) log10(.data[[outcome]]) else .data[[outcome]],
             visit  = factor(paste0("V", .data[[time]]), levels = visit_levels)) |>
      filter(is.finite(.y_obs), is.finite(x_obs)) |>
      rename(!!v_lhs := .y_obs)
    
    n_patients[[v]] <- n_distinct(obs_data[[v]][[id]])
    
    nd <- do.call(rbind, lapply(time_vals, function(t_val) {
      nd_v <- data.frame(matrix(ncol = 0, nrow = n_pred))
      nd_v[[time]] <- t_val; nd_v[[age]] <- age_m; nd_v[[outcome]] <- gdf15_raw_seq
      nd_v[["visit"]] <- paste0("V", t_val); nd_v[["x_pred"]] <- x_pred_seq
      if (!is.null(lithium_mode)) nd_v[[lithium]] <- lithium_mode
      nd_v
    }))
    nd$pred <- predict(model, newdata = nd, re.form = NA)
    pred_data[[v]] <- nd |> mutate(visit = factor(visit, levels = visit_levels))
  }
  
  vars_ok <- names(model_stats)
  if (length(vars_ok) == 0) { message("No valid models."); return(invisible(NULL)) }
  p_fdr <- setNames(p.adjust(sapply(vars_ok, function(v) model_stats[[v]]$p), method = "fdr"), vars_ok)
  
  panel_list <- lapply(vars_ok, function(v) {
    s     <- model_stats[[v]]; nd <- pred_data[[v]]; v_lhs <- s$v_lhs
    sub_txt <- sprintf("N = %d  \u00b7  \u03b2(time\u00d7GDF15) = %.4f  \u00b7  p = %.4f%s  \u00b7  pFDR = %.4f%s",
                       n_patients[[v]], s$beta, s$p, sig_symbol(s$p), p_fdr[[v]], sig_symbol(p_fdr[[v]]))  
    
    ggplot() +
      { if (show_obs) geom_point(data = obs_data[[v]], aes(x_obs, .data[[v_lhs]], colour = visit, shape = visit), alpha = 0.25, size = 1.3) else list() } +
      geom_line(data = nd, aes(x_pred, pred, colour = visit, linetype = visit), linewidth = 1.1) +
      scale_colour_manual(values = visit_cols, name = "Visit") +
      scale_linetype_manual(values = setNames(rep("dashed", length(time_vals)), visit_levels), name = "Visit") +
      scale_shape_manual(values   = setNames(c(16L, 17L, 15L)[seq_along(time_vals)], visit_levels), name = "Visit") +
      labs(title    = sprintf("Effect of GDF15 [log10(pg/mL)] on trajectory of %s", v),
           subtitle = sub_txt,
           x = if (log_outcome) bquote(log[10](.(as.name(outcome)))) else outcome, y = v_lhs) +
      theme_classic(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 12, colour = facebd_colors$dark_navy),
            plot.subtitle = element_text(size = 9.5, colour = "grey35"))
  }) |> setNames(vars_ok)
  
  wrap_plots(panel_list, ncol = ncol) + plot_layout(guides = "collect") +
    plot_annotation(
      title   = "Clinical variable trajectories — mixed-model (GDF15 as predictor)",
      caption = "Adjusted for age and lithium treatment",
      theme   = theme(plot.title   = element_text(face = "bold", size = 13, colour = facebd_colors$dark_burgundy, hjust = 0.5),
                      plot.caption = element_text(size = 10, colour = "grey40", face = "italic", hjust = 0.5))) &
    theme(legend.position = "bottom", legend.text = element_text(size = 12),
          legend.key.width = unit(1.6, "cm"), legend.key.height = unit(0.6, "cm"))
}

# ══════════════════════════════════════════════════════════════════════════════
# 10. run_delta_regressions() — V1 clinical outcome prediction from baseline GDF15
# ══════════════════════════════════════════════════════════════════════════════

run_delta_regressions <- function(data, vars, gdf15 = "GDF15", log_gdf15 = TRUE,
                                  id = "fondacode", time = "time_num",
                                  age = "age", transform_vars = FALSE,
                                  adjust_lithium = FALSE, lithium_col = "lithium_treat",
                                  fdr_method = "BH") {
  
  gdf15_col <- if (log_gdf15) "log_GDF15" else gdf15
  
  # Baseline GDF15 + covariates at V0
  df_v0 <- data |>
    filter(.data[[time]] == 0) |>
    select(all_of(c(id, gdf15_col, age, if (adjust_lithium) lithium_col))) |>
    distinct()
  
  results <- lapply(vars, function(v) {
    
    # Extract V0 and V1 values per patient
    tmp <- data |>
      filter(.data[[time]] %in% c(0, 1), !is.na(.data[[v]])) |>
      select(all_of(c(id, time, v))) |>
      pivot_wider(id_cols = all_of(id), names_from = all_of(time),
                  values_from = all_of(v), names_prefix = "t") |>
      inner_join(df_v0, by = id) |>
      tidyr::drop_na()
    
    if (!all(c("t0", "t1") %in% names(tmp))) {
      warning(sprintf("Missing V0 or V1 for variable '%s'. Skipping.", v)); return(NULL)
    }
    
    tmp$delta <- tmp$t1 - tmp$t0
    if (transform_vars) tmp$delta <- sign(tmp$delta) * sqrt(abs(tmp$delta))
    
    n <- nrow(tmp)
    if (n < 10) { warning(sprintf("Not enough observations for '%s'.", v)); return(NULL) }
    
    covs        <- c(gdf15_col, age, if (adjust_lithium) lithium_col)
    formula_str <- paste("delta ~", paste(covs, collapse = " + "))
    fit <- tryCatch(lm(as.formula(formula_str), data = tmp), error = function(e) NULL)
    if (is.null(fit)) { warning(sprintf("Model failed for '%s'.", v)); return(NULL) }
    
    ct     <- summary(fit)$coefficients
    r2     <- summary(fit)$r.squared
    adj_r2 <- summary(fit)$adj.r.squared
    beta   <- ct[gdf15_col, "Estimate"];   se   <- ct[gdf15_col, "Std. Error"]
    tval   <- ct[gdf15_col, "t value"];    pval <- ct[gdf15_col, "Pr(>|t|)"]
    ci_lo  <- beta - qt(0.975, fit$df.residual) * se
    ci_hi  <- beta + qt(0.975, fit$df.residual) * se
    
    tibble(
      window    = "\u0394 (V1 - V0)",
      variable  = v,
      transform = if (transform_vars) "sqrt(|\u0394|)\u00b7sign(\u0394)" else "identity",
      n = n, beta = round(beta, 4), SE = round(se, 4),
      CI_low = round(ci_lo, 4), CI_high = round(ci_hi, 4),
      t = round(tval, 3), p_value = round(pval, 4),
      R2 = round(r2, 4), adj_R2 = round(adj_r2, 4)
    )
  })
  
  out <- bind_rows(results)
  if (nrow(out) == 0) { message("No valid results."); return(out) }
  
  out |>
    mutate(
      p_adj_fdr = round(p.adjust(p_value, method = fdr_method), 4),
      sig_raw   = sig_stars(p_value),
      sig_fdr   = sig_stars(p_adj_fdr)
    ) |>
    arrange(p_value)
}

# ══════════════════════════════════════════════════════════════════════════════
# 11. run_logistic_gdf15() — Binary outcome prediction
# ══════════════════════════════════════════════════════════════════════════════

run_logistic_gdf15 <- function(df, clinical_vars, visit_predictor = "V0", visit_outcome = "V1_an",
                               gdf15_col = "log_GDF15", visit_col = "visit", subject_col = "fondacode",
                               q = 0.75, adjust_baseline = TRUE, visit_baseline = NULL, covariates = NULL) {
  if (is.null(visit_baseline)) visit_baseline <- visit_predictor
  df_pred <- df |> filter(.data[[visit_col]] == visit_predictor) |> select(all_of(c(subject_col, gdf15_col))) |> rename(gdf15_pred = all_of(gdf15_col))
  df_out  <- df |> filter(.data[[visit_col]] == visit_outcome)  |> select(all_of(c(subject_col, clinical_vars)))
  df_base <- df |> filter(.data[[visit_col]] == visit_baseline) |> select(all_of(c(subject_col, clinical_vars))) |>
    rename_with(~ paste0(.x, "_baseline"), all_of(clinical_vars))
  
  results <- purrr::map_dfr(clinical_vars, function(var) {
    na_row <- function(note) tibble(variable = var, threshold = NA_real_, n_severe = NA_integer_, n_total = n_total,
                                    OR_sd = NA_real_, CI_low_sd = NA_real_, CI_high_sd = NA_real_,
                                    p_value = NA_real_, AUC = NA_real_, gdf15_sd_log10 = NA_real_,
                                    gdf15_1sd_pgml = NA_real_, converged = FALSE, note = note)
    tmp     <- df_out |> select(all_of(c(subject_col, var))) |> inner_join(df_pred, by = subject_col)
    base_col <- paste0(var, "_baseline")
    use_base <- adjust_baseline && base_col %in% names(df_base)
    if (use_base) tmp <- tmp |> left_join(df_base |> select(all_of(c(subject_col, base_col))), by = subject_col)
    tmp <- tmp |> tidyr::drop_na(); n_total <- nrow(tmp)
    if (n_total < 20) return(na_row("Not enough observations"))
    
    tmp[[var]] <- sqrt(tmp[[var]]); if (use_base) tmp[[base_col]] <- sqrt(tmp[[base_col]])
    thr_sqrt   <- quantile(tmp[[var]], probs = q, na.rm = TRUE); thr_orig <- thr_sqrt^2
    tmp$y      <- as.integer(tmp[[var]] > thr_sqrt); n_severe <- sum(tmp$y)
    if (n_severe < 10 || (n_total - n_severe) < 10) return(na_row(paste0("Too few events (n_severe = ", n_severe, ")")))
    
    gdf15_mean     <- mean(tmp$gdf15_pred, na.rm = TRUE); gdf15_sd <- sd(tmp$gdf15_pred, na.rm = TRUE)
    gdf15_1sd_pgml <- round(10^(gdf15_mean + gdf15_sd) - 10^gdf15_mean, 0)
    tmp$gdf15_pred_sd <- scale(tmp$gdf15_pred)[, 1]
    
    formula_sd <- as.formula(paste("y ~", paste(c("gdf15_pred_sd", if (use_base) base_col,
                                                  covariates[covariates %in% names(tmp)]), collapse = " + ")))
    fit <- tryCatch(glm(formula_sd, data = tmp, family = binomial()), error = function(e) NULL)
    if (is.null(fit)) return(na_row("Model did not converge"))
    
    ct      <- broom::tidy(fit, conf.int = TRUE, conf.method = "profile", exponentiate = TRUE)
    gdf_row <- ct |> filter(term == "gdf15_pred_sd")
    roc_obj <- pROC::roc(tmp$y, predict(fit, type = "response"), quiet = TRUE)
    
    tibble(variable = var, threshold = round(thr_orig, 2), n_severe = as.integer(n_severe), n_total = as.integer(n_total),
           OR_sd = round(gdf_row$estimate, 3), CI_low_sd = round(gdf_row$conf.low, 3), CI_high_sd = round(gdf_row$conf.high, 3),
           p_value = round(gdf_row$p.value, 4), AUC = round(as.numeric(pROC::auc(roc_obj)), 3),
           gdf15_sd_log10 = round(gdf15_sd, 3), gdf15_1sd_pgml = gdf15_1sd_pgml, converged = TRUE,
           note = paste0(if (use_base) paste0("adjusted for ", visit_baseline, " score") else "",
                         if (length(covariates[covariates %in% names(tmp)]) > 0) paste0(" + ", paste(covariates[covariates %in% names(tmp)], collapse = ", ")) else ""))
  })
  
  valid_p <- !is.na(results$p_value)
  if (sum(valid_p) > 1) {
    results$p_value_adj_BH <- NA_real_
    results$p_value_adj_BH[valid_p] <- round(p.adjust(results$p_value[valid_p], method = "BH"), 4)
  }
  results |> arrange(p_value)
}

# ══════════════════════════════════════════════════════════════════════════════
# 12. EXPORT TABLES
# ══════════════════════════════════════════════════════════════════════════════

make_gt_theme <- function(tbl, dark_burgundy = "#6E123B", light_pink = "#F2B8CC",
                          light_gray = "#F7F7F7", border_gray = "#D0D0D0") {
  tbl |>
    cols_align(align = "center", columns = everything()) |>
    tab_style(style = list(cell_fill(color = dark_burgundy), cell_text(color = "white", weight = "bold")),
              locations = cells_column_labels()) |>
    tab_style(style = list(cell_fill(color = light_pink), cell_text(color = dark_burgundy, weight = "bold")),
              locations = cells_row_groups()) |>
    opt_row_striping() |>
    tab_options(row.striping.background_color = light_gray, table.border.top.color = border_gray,
                column_labels.border.bottom.color = border_gray, table.width = pct(100), table.font.size = 13)
}

export_mundlak_table <- function(results_mundlak, filename,
                                title    = "GDF15 [log10(pg/mL)] within- and between-person effects on clinical scores",
                                subtitle = "Adjusted for age and lithium treatment") {
  library(gt); dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  
  sig_fn <- function(p) case_when(
    !is.na(p) & p < 0.001 ~ "***", !is.na(p) & p < 0.01 ~ "**",
    !is.na(p) & p < 0.05  ~ "*",   !is.na(p) & p < 0.10  ~ ".", TRUE ~ ""
  )
  
  make_block <- function(effect, label) {
    results_mundlak |>
      filter(!is.na(.data[[paste0("beta_", effect)]])) |>
      transmute(
        effect_label  = label,
        variable,
        n,
        `β [95% CI]`  = paste0(
          sprintf("%.4f", .data[[paste0("beta_",    effect)]]), "  [",
          sprintf("%.4f", .data[[paste0("CI_low_",  effect)]]), " ; ",
          sprintf("%.4f", .data[[paste0("CI_high_", effect)]]), "]"
        ),
        t         = .data[[paste0("t_", effect)]],
        p_value   = round(.data[[paste0("p_",          effect)]], 4),
        p_adj_fdr = round(.data[[paste0("p_adj_fdr_",  effect)]], 4),
        AIC, BIC, ICC,
        sig       = sig_fn(.data[[paste0("p_",         effect)]]),
        sig_fdr   = sig_fn(.data[[paste0("p_adj_fdr_", effect)]])
      ) |>
      arrange(p_value)
  }
  
  tbl_data <- bind_rows(
    make_block("within",  "Within-person effect  (GDF15_dev)"),
    make_block("between", "Between-person effect  (GDF15_mean)")
  )
  
  tbl <- gt(tbl_data, groupname_col = "effect_label") |>
    tab_header(
      title    = md(paste0("**", title, "**")),
      subtitle = md(paste0("<span style='color:#888888; font-size:13px;'>", subtitle, "</span>"))
    ) |>
    cols_label(
      variable    = "Clinical variable",
      n           = "N",
      `β [95% CI]` = "\u03b2 [95% CI]",
      t           = "t",
      p_value     = "p",
      p_adj_fdr   = "FDR p",
      AIC         = "AIC", BIC = "BIC", ICC = "ICC",
      sig         = "Sig.", sig_fdr = "Sig. (FDR)"
    ) |>
    cols_align(align = "left",   columns = c(variable, `β [95% CI]`)) |>
    cols_align(align = "center", columns = c(n, sig, sig_fdr)) |>
    # ── pink row-group headers ──────────────────────────────────────────────
    tab_style(
      style = list(
        cell_fill(color = "#FAD4E0"),
        cell_text(color = "#6E123B", weight = "bold", size = px(13))
      ),
      locations = cells_row_groups()
    ) |>
    make_gt_theme() |>
    
    tab_footnote(
      footnote  = "\u00b7 p < 0.10  * p < 0.05  ** p < 0.01  *** p < 0.001",
      locations = cells_column_labels(columns = sig)
    )
  
  print(tbl); gtsave(tbl, filename = filename)
  cat("Table exported to:", filename, "\n"); invisible(tbl)
}

export_logistic_table <- function(results_all, filename, title = "GDF15 prognostic models") {
  library(gt); dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  sig_fn <- function(p, v = p) case_when(!is.na(v) & v < 0.001 ~ "***", !is.na(v) & v < 0.01 ~ "**",
                                         !is.na(v) & v < 0.05 ~ "*",   !is.na(v) & v < 0.10 ~ "·", TRUE ~ "")
  tbl_data <- results_all |>
    mutate(group = "V0 → V2", sig = sig_fn(p_value), sig_adj = sig_fn(p_value_adj_BH),
           `OR [95% CI]`   = ifelse(is.na(OR_sd), "NA", paste0(sprintf("%.2f", OR_sd), " [", sprintf("%.2f", CI_low_sd), " ; ", sprintf("%.2f", CI_high_sd), "]")),
           `1 SD (pg/mL)`  = paste0("+", gdf15_1sd_pgml, " pg/mL"),
           threshold = round(threshold, 2), AUC = round(AUC, 3)) |>
    select(group, variable, threshold, n_total, `OR [95% CI]`, `1 SD (pg/mL)`, p_value, p_value_adj_BH, AUC, sig, sig_adj)
  
  tbl <- gt(tbl_data, groupname_col = "group") |>
    tab_header(title = md(paste0("**", title, "**"))) |>
    cols_label(variable = "Clinical variable (sqrt)", threshold = "Threshold (orig.)", n_total = "N",
               `OR [95% CI]` = "OR per SD [95% CI]", `1 SD (pg/mL)` = "1 SD GDF15",
               p_value = "p", p_value_adj_BH = "FDR p", AUC = "AUC", sig = md("Sig."), sig_adj = md("Sig. (FDR)")) |>
    make_gt_theme() |>
    tab_style(style = cell_text(color = "#6E123B", weight = "bold"),
              locations = cells_body(columns = c(sig, sig_adj), rows = sig %in% c("·","*","**","***") | sig_adj %in% c("·","*","**","***"))) |>
    tab_footnote(footnote = md("· p < 0.10 &nbsp; * p < 0.05 &nbsp; ** p < 0.01 &nbsp; *** p < 0.001"), locations = cells_column_labels(columns = sig))
  
  print(tbl); gtsave(tbl, filename = filename)
  cat("Table exported to:", filename, "\n"); invisible(tbl)
}

export_delta_table <- function(results_delta, filename,
                               title = "GDF15 [log10(pg/mL)] predicting change in clinical scores",
                               subtitle = "Adjusted for age and lithium treatment",
                               transform = "sqrt(|Δ|)·sign(Δ)") {
  library(gt); dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  has_sqrt  <- any(results_delta$transform == "sqrt(|Δ|)·sign(Δ)", na.rm = TRUE)
  sig_fn    <- function(p) case_when(!is.na(p) & p < 0.001 ~ "***", !is.na(p) & p < 0.01 ~ "**",
                                     !is.na(p) & p < 0.05 ~ "*",   !is.na(p) & p < 0.10 ~ ".", TRUE ~ "")
  tbl_data  <- results_delta |>
    filter(!is.na(beta)) |>
    mutate(window = "\u0394 (V1 - V0)",
           `β [95% CI]` = paste0(sprintf("%.3f", beta), "  [", sprintf("%.3f", CI_low), " ; ", sprintf("%.3f", CI_high), "]"),
           sig = sig_fn(p_value), sig_fdr = sig_fn(p_adj_fdr),
           across(c(p_value, p_adj_fdr, R2, adj_R2), round, 4)) |>
    arrange(p_value) |>
    select(window, variable, n, `β [95% CI]`, t, p_value, p_adj_fdr, R2, adj_R2, sig, sig_fdr)
  
  tbl <- gt(tbl_data, groupname_col = "window") |>
    tab_header(title = md(paste0("**", title, "**")),
               subtitle = md(paste0("<span style='color:#888888; font-size:13px;'>", subtitle, "</span>"))) |>
    cols_label(variable = if (has_sqrt) "Clinical variable (sqrt)" else "Clinical variable",
               n = "N", `β [95% CI]` = "\u03b2 [95% CI]", t = "t",
               p_value = "p", p_adj_fdr = "FDR p", R2 = "R\u00b2", adj_R2 = "Adj. R\u00b2",
               sig = "Sig.", sig_fdr = "Sig. (FDR)") |>
    cols_align(align = "left", columns = c(variable, `β [95% CI]`)) |>
    make_gt_theme() |>
    tab_style(style = cell_text(color = "#6E123B", weight = "bold"),
              locations = cells_body(columns = c(p_value, p_adj_fdr, sig, sig_fdr), rows = p_value < 0.05)) |>
    tab_style(style = cell_text(color = "#A83060"),
              locations = cells_body(columns = c(p_value, sig), rows = p_value >= 0.05 & p_value < 0.10)) |>
    tab_footnote(footnote = "\u00b7 p < 0.10  * p < 0.05  ** p < 0.01  *** p < 0.001",
                 locations = cells_column_labels(columns = sig))
  
  print(tbl); gtsave(tbl, filename = filename)
  cat("Table exported to:", filename, "\n"); invisible(tbl)
}

# ══════════════════════════════════════════════════════════════════════════════
# 13. ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

df_filtered <- df_filtered |>
  rename(QIDS = qidsr120, PSQI = psqi_, MADRS = madrs_, FAST = fast_, STAI = staya, MARS = mars_, ALS = als_) |>
  filter(!is.na(lithium_treat))

vars_clinical <- c("PSQI", "MADRS", "FAST", "STAI")

hist(df_filtered$log_GDF15[df_filtered$time_num == 0])

# Mixed models: clinical ~ GDF15 trajectory
sqrt_tr <- setNames(rep(list("sqrt"), length(vars_clinical)), vars_clinical)

results_clinical <- run_clinical_on_GDF15(data = df_filtered, vars = vars_clinical, transformations = sqrt_tr)

results_clinical |> select(variable, transform, n_skipped_lithium, family, beta_slope, SE_slope,
                           CI_low_slope, CI_high_slope, t_slope, p_slope, sig_slope_raw, p_adj_fdr_slope, sig_slope_fdr) |> print(n = Inf, width = Inf)
results_clinical |> select(variable, transform, beta_level, SE_level, p_level, sig_level_raw, p_adj_fdr_level, sig_level_fdr) |> print(n = Inf, width = Inf)
results_clinical |> select(variable, transform, shapiro_resid_p, shapiro_re_p, DW_stat, DW_p, n_outliers, converged) |> print(n = Inf, width = Inf)

combined_plot <- plot_mixed_results(data = df_filtered, results = results_clinical,
                                    filter_by = "p_slope", plot_alpha = 0.05, ncol = 2,
                                    lithium = "lithium_treat")
print(combined_plot)
ggsave("outcome/mixed_model_color_clinical_GDF15_numbers.png", combined_plot, width = 12, height = 8, dpi = 300, bg = "white")

plot_gdf15_xaxis <- plot_mixed_gdf15_xaxis(data = df_filtered, results = results_clinical, outcome = "GDF15", show_obs = TRUE, ncol = 2)
print(plot_gdf15_xaxis)
ggsave("outcome/mixed_model_color_clinical_GDF15_xaxis_numbers.png", plot_gdf15_xaxis, width = 12, height = 8, dpi = 300, bg = "white")

# Mundlak within/between decomposition
results_mundlak <- run_clinical_on_GDF15_mundlak(data = df_filtered, vars = vars_clinical,
                                               transformations = sqrt_tr, adjust_lithium = TRUE)
results_mundlak |> select(variable, transform, family, beta_within, SE_within, CI_low_within, CI_high_within,
                         t_within, p_within, sig_within_raw, p_adj_fdr_within, sig_within_fdr) |> print(n = Inf, width = Inf)
results_mundlak |> select(variable, transform, beta_between, SE_between, p_between, sig_between_raw, p_adj_fdr_between, sig_between_fdr) |> print(n = Inf, width = Inf)
results_mundlak |> select(variable, transform, shapiro_resid_p, shapiro_re_p, DW_stat, DW_p, n_outliers, converged) |> print(n = Inf, width = Inf)

export_mundlak_table(results_mundlak, "outcome/gdf15_mundlak_between.png")

# Delta regressions: GDF15 V0 predicting change V1-V0
results_delta_regression <- run_delta_regressions(data = df_filtered, vars = vars_clinical,
                                               transform_vars = TRUE, adjust_lithium = TRUE,
                                               lithium_col = "lithium_treat")
results_delta_regression |> select(window, variable, n, beta, SE, CI_low, CI_high, t, p_value, sig_raw, p_adj_fdr, sig_fdr, R2, adj_R2) |> print(n = Inf)
results_delta_regression |> filter(p_value < 0.05)

export_delta_table(results_delta = results_delta_regression,
                   filename      = "outcome/GDF15_delta_regressions.png",
                   title         = "GDF15 [log10(pg/mL)] predicting change in clinical scores")

# Logistic models: GDF15 at V0 predicting severity at V2
results_v0_v2 <- run_logistic_gdf15(df = df_filtered, clinical_vars = vars_clinical,
                                    visit_predictor = "V0", visit_outcome = "V2_ans",
                                    visit_baseline = "V0", adjust_baseline = TRUE)
print(results_v0_v2)
results_v0_v2 |> select(variable, p_value_adj_BH, AUC)

export_logistic_table(results_all = results_v0_v2,
                               filename    = "outcome/GDF15_logistic_combined.png",
                               title       = "GDF15 [log10(pg/mL)] prognostic models across time")

