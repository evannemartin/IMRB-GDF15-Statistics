# ── Libraries ─────────────────────────────────────────────
library(readxl); library(ggplot2); library(reshape2)
library(dplyr); library(ggrepel); library(tidyverse)
library(lmtest); library(car); library(MASS)
library(multcomp); library(multcompView)
library(effectsize); library(patchwork); library(cowplot)

setwd("~/IMRB/Stats GDF15/")

# ── FACE-BD Color palette ─────────────────────────────────
facebd_colors <- list(
  dark_burgundy  = "#6B1A3A",  # figures foncées
  medium_rose    = "#A83060",  # figures moyennes
  light_rose     = "#D4789A",  # figures claires
  very_light_pink = "#F2B8CC", # figures très claires
  blue_gray      = "#8899AA",  # figures contrôles
  dark_navy      = "#2D3B6E"   # texte FACE-BD
)

# Palette séquentielle pour heatmap (bleu → blanc → bordeaux)
facebd_gradient <- list(
  low  = "#2D3B6E",  # bleu navy (négatif)
  mid  = "white",
  high = "#A83060"   # rose moyen = couleur Inflammation (positif max)
)

# Palette catégorielle pour groupes
facebd_groups <- c(
  "#6B1A3A",  # Clinical      — bordeaux foncé
  "#8899AA",  # Blood cells   — bleu-gris
  "#A83060",  # Inflammation  — rose moyen
  "#2D3B6E",  # KYN pathway   — bleu navy
  "#D4789A"   # Metabolism    — rose clair
)

# ════════════════════════════════════════════════════════
# DATA IMPORT & CLEANING
# ════════════════════════════════════════════════════════

df_trans <- read_excel("data/final_merged_data_samples_V0.xlsx")

# Remove columns with same values (limit of detection)
df_trans <- subset(df_trans, select = -c(IL_4, IL_2, IL_1beta, IL_1alpha))

# Convert "," to "." and "</>X" to NA
df_trans <- df_trans %>%
  mutate(
    across(where(~ is.character(.x) && any(grepl(",",    .x, fixed = TRUE))),
           ~ as.numeric(gsub(",", ".", .x))),
    across(where(~ is.character(.x) && any(grepl("^[<>]", .x))),
           ~ as.numeric(ifelse(grepl("^[<>]", .x), NA, gsub(",", ".", .x))))
  )

# Replace ambiguous strings with NA
na_strings <- c("BD not specified", "Ne sais pas", "ne sais pas",
                "BD NOT SPECIFIED", "Unknown", "unknown", "NA", "N/A", "n/a")
df_trans <- df_trans %>%
  mutate(across(where(~ is.factor(.) | is.character(.)),
                ~ ifelse(trimws(.x) %in% na_strings, NA, .x)))

# Re-check character columns that could now be numeric (keep binary as categorical)
char_cols <- names(df_trans)[sapply(df_trans, is.character)]
converted <- c(); kept_binary <- c()

for (col in char_cols) {
  x <- df_trans[[col]]; non_na <- x[!is.na(x)]
  converted_vals <- suppressWarnings(as.numeric(gsub(",", ".", non_na)))
  if (length(non_na) > 0 && !any(is.na(converted_vals))) {
    if (all(unique(converted_vals) %in% c(0, 1))) {
      kept_binary <- c(kept_binary, col)
    } else {
      df_trans[[col]] <- suppressWarnings(as.numeric(gsub(",", ".", x)))
      converted <- c(converted, col)
    }
  }
}
cat("Converted to numeric:", paste(converted, collapse = ", "), "\n")
cat("Kept as binary:", paste(kept_binary, collapse = ", "), "\n")

table(sapply(df_trans, class))

# Smoking variables
df_trans$`current smokers`  <- ifelse(df_trans$suncf_cigarettes_lt == "Yes",      "Yes", "No")
df_trans$`remitted smokers` <- ifelse(df_trans$suncf_cigarettes_lt == "Remitted", "Yes", "No")

# ════════════════════════════════════════════════════════
# FUNCTIONS
# ════════════════════════════════════════════════════════

###################### SPEARMAN CORRELATION MATRIX ######################

test_spearman_pairwise <- function(df, alpha = 0.05, check_autocor = FALSE,
                                   maha_threshold = 0.975, id_col = "fondacode",
                                   p_adjust_method = "BH") {
  num_vars  <- df[, sapply(df, is.numeric)]
  var_names <- names(num_vars); n_vars <- length(var_names)
  n_tests   <- n_vars * (n_vars - 1) / 2
  cat("══════════════════════════════════════════\n")
  cat(n_vars, "numeric variables | Total pairs:", n_tests,
      "| Adjustment:", p_adjust_method, "\n")
  cat(paste("->", var_names, collapse = "\n"), "\n")
  cat("══════════════════════════════════════════\n\n")
  results <- list()
  for (i in 1:(n_vars - 1)) {
    for (j in (i + 1):n_vars) {
      x <- num_vars[[i]]; y <- num_vars[[j]]
      pair <- paste0(var_names[i], " ~ ", var_names[j])
      valid <- complete.cases(x, y)
      x <- as.numeric(x[valid]); y <- as.numeric(y[valid])
      df_valid <- df[valid, ]; n <- length(x)
      if (n < 8 || sd(x) == 0 || sd(y) == 0) {
        results[[paste0(var_names[i], "_", var_names[j])]] <- data.frame(
          pair = pair, n = n, n_outliers = NA, independence_ok = FALSE,
          spearman_rho = NA, spearman_p = NA, spearman_p_adj = NA,
          rho_without_out = NA, rho_difference = NA, outlier_influence = NA,
          spearman_ok = FALSE, significance = NA,
          advice = "ERROR: n too small or zero variance", error = "n too small or zero variance")
        next
      }
      rx <- rank(x); ry <- rank(y)
      n_out <- NA; rho_without_out <- NA; rho_difference <- NA; outlier_influence <- NA
      mat <- cbind(rx, ry); cov_mat <- cov(mat)
      if (det(cov_mat) >= 1e-10) {
        maha <- mahalanobis(mat, colMeans(mat), cov_mat)
        idx  <- which(maha > qchisq(maha_threshold, df = 2)); n_out <- length(idx)
        if (n_out > 0) {
          outlier_labels <- if (id_col %in% names(df)) as.character(df_valid[[id_col]][idx]) else
            paste0("row:", which(valid)[idx])
          rho_with   <- cor(x, y, method = "spearman")
          x_no <- x[-idx]; y_no <- y[-idx]
          if (length(x_no) < 3 || sd(x_no) == 0 || sd(y_no) == 0) {
            outlier_influence <- "Cannot compute — not enough obs after outlier removal"
          } else {
            rho_without_out <- cor(x_no, y_no, method = "spearman")
            rho_difference  <- abs(rho_with - rho_without_out)
            outlier_influence <- ifelse(rho_difference < 0.05, "Low — outlier(s) not influential",
                                        ifelse(rho_difference < 0.10, "Moderate — mention in methods",
                                               "High — investigate outlier(s)"))
          }
          cat("\n──────────────────────────────────────────\n")
          cat("Pair:", pair, "\nOutlier(s):", paste(outlier_labels, collapse = ", "), "\n")
          cat("rho with:", round(rho_with, 4), "| without:", round(rho_without_out, 4),
              "| diff:", round(rho_difference, 4), "| influence:", outlier_influence, "\n")
          cat("──────────────────────────────────────────\n")
        }
      } else cat("Warning: singular covariance matrix for:", pair, "\n")
      independ  <- if (check_autocor) dwtest(lm(ry ~ rx))$p.value > alpha else TRUE
      spearman  <- cor.test(x, y, method = "spearman", exact = FALSE)
      spearman_ok <- independ & !is.na(n_out) & n_out == 0
      advice <- if (is.na(n_out)) "Singular matrix — check pair manually" else
        if (spearman_ok) "Spearman valid" else
          if (!independ)   "Warning: autocorrelation detected" else
            paste0("Check ", n_out, " outlier(s) — influence: ", outlier_influence)
      results[[paste0(var_names[i], "_", var_names[j])]] <- data.frame(
        pair = pair, n = n, n_outliers = n_out, independence_ok = independ,
        spearman_rho = round(spearman$estimate, 4), spearman_p = round(spearman$p.value, 4),
        spearman_p_adj = NA, rho_without_out = round(rho_without_out, 4),
        rho_difference = round(rho_difference, 4), outlier_influence = outlier_influence,
        spearman_ok = spearman_ok, significance = NA, advice = advice, error = NA_character_)
    }
  }
  df_results <- do.call(rbind, results)
  valid_p <- !is.na(df_results$spearman_p)
  df_results$spearman_p_adj[valid_p] <- p.adjust(df_results$spearman_p[valid_p], method = p_adjust_method)
  df_results$significance <- ifelse(is.na(df_results$spearman_p_adj), NA,
                                    ifelse(df_results$spearman_p_adj < 0.001, "***",
                                           ifelse(df_results$spearman_p_adj < 0.01,  "**",
                                                  ifelse(df_results$spearman_p_adj < 0.05,  "*", "ns"))))
  df_results <- df_results %>% arrange(spearman_p_adj)
  cat("\n══════════════════════════════════════════\n")
  cat("Significant (raw p):", sum(df_results$spearman_p < alpha, na.rm = TRUE),
      "| Significant (adj p):", sum(df_results$spearman_p_adj < alpha, na.rm = TRUE),
      "| Pairs with outliers:", sum(df_results$n_outliers > 0, na.rm = TRUE), "\n")
  cat("══════════════════════════════════════════\n")
  return(df_results)
}

###################### SPEARMAN CORRELATION HEATMAP ######################

plot_spearman_heatmap <- function(df, results, alpha = 0.05, focus_var = NULL,
                                  focus_pairs = NULL, use_adjusted_p = TRUE,
                                  strip_plot = FALSE, strip_var = NULL,
                                  rho_limits = NULL) {
  p_col <- if (use_adjusted_p && "spearman_p_adj" %in% names(results)) "spearman_p_adj" else "spearman_p"
  if (strip_plot) {
    if (is.null(strip_var)) stop("strip_var must be specified for strip_plot mode.")
    strip_df <- results %>%
      filter(sapply(strsplit(pair, " ~ "), `[`, 1) == strip_var |
               sapply(strsplit(pair, " ~ "), `[`, 2) == strip_var) %>%
      mutate(
        other_var = ifelse(sapply(strsplit(pair, " ~ "), `[`, 1) == strip_var,
                           sapply(strsplit(pair, " ~ "), `[`, 2),
                           sapply(strsplit(pair, " ~ "), `[`, 1)),
        p_use = .data[[p_col]],
        stars = ifelse(p_use < 0.001, "***", ifelse(p_use < 0.01, "**", ifelse(p_use < 0.05, "*", ""))),
        label = paste0(ifelse(round(spearman_rho, 3) == 0,
                              format(round(spearman_rho, 4), nsmall = 4),
                              round(spearman_rho, 3)), stars),
        sig       = p_use < alpha,
        other_var = reorder(other_var, spearman_rho), 
        x_fixed   = factor("x"))
    if (nrow(strip_df) == 0) { cat("No pairs found involving", strip_var, "\n"); return(invisible(NULL)) }
    p_label <- ifelse(use_adjusted_p & "spearman_p_adj" %in% names(results), "FDR-adjusted p", "raw p")
    p <- ggplot(strip_df, aes(x = x_fixed, y = other_var)) +
      geom_tile(aes(fill = spearman_rho), color = "white", linewidth = 0.5) +
      geom_text(aes(label = label, fontface = ifelse(sig, "bold", "plain")), size = 6.2, color = "black") +
      scale_fill_gradient2(
        low      = facebd_gradient$low,   # "#2D3B6E"
        mid      = facebd_gradient$mid,   # "white"
        high     = facebd_gradient$high,  # "#6B1A3A"
        midpoint = 0,
        limits   = rho_limits,
        name     = "Spearman ρ"
      ) +
      labs(title = paste("Spearman correlations with", strip_var),
           subtitle = paste0("Stars based on ", p_label, ": * p<.05  ** p<.01  *** p<.001 | Biological variables in pg/mL"),
           x = strip_var, y = NULL) +
      theme_minimal(base_size = 16) +
      theme(plot.title = element_text(face = "bold", size = 14),
            plot.subtitle = element_text(size = 9, color = "grey40"),
            axis.text.y = element_text(face = "bold", size = 16, hjust = 1),
            axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            panel.grid = element_blank(), legend.position = "right",
            legend.title = element_text(size = 14, face = "bold"),
            legend.text = element_text(size = 14), legend.key.size = unit(1, "cm"))
    print(p); return(invisible(p))
  }
  # Standard heatmap
  if (!is.null(focus_pairs)) {
    sig_pairs <- results %>%
      filter(sapply(strsplit(pair, " ~ "), `[`, 1) %in% focus_pairs &
               sapply(strsplit(pair, " ~ "), `[`, 2) %in% focus_pairs)
    sig_vars <- focus_pairs
  } else {
    sig_pairs <- results %>% filter(.data[[p_col]] < alpha)
    if (!is.null(focus_var)) sig_pairs <- sig_pairs %>%
        filter(sapply(strsplit(pair, " ~ "), `[`, 1) == focus_var |
                 sapply(strsplit(pair, " ~ "), `[`, 2) == focus_var)
    sig_vars <- unique(c(sapply(strsplit(sig_pairs$pair, " ~ "), `[`, 1),
                         sapply(strsplit(sig_pairs$pair, " ~ "), `[`, 2)))
  }
  if (nrow(sig_pairs) == 0) { cat("No significant pairs found.\n"); return(invisible(NULL)) }
  n_sig <- length(sig_vars)
  rho_mat   <- matrix(NA,    n_sig, n_sig, dimnames = list(sig_vars, sig_vars))
  valid_mat <- matrix(FALSE, n_sig, n_sig, dimnames = list(sig_vars, sig_vars))
  diag(rho_mat) <- 1; diag(valid_mat) <- TRUE
  for (k in seq_len(nrow(results))) {
    parts <- strsplit(results$pair[k], " ~ ")[[1]]; vi <- parts[1]; vj <- parts[2]
    if (!(vi %in% sig_vars) || !(vj %in% sig_vars)) next
    rho_mat[vi, vj] <- rho_mat[vj, vi] <- results$spearman_rho[k]
    valid_mat[vi, vj] <- valid_mat[vj, vi] <- results$spearman_ok[k]
  }
  rho_melt   <- melt(rho_mat,   varnames = c("Var1", "Var2"), value.name = "rho")
  valid_melt <- melt(valid_mat, varnames = c("Var1", "Var2"), value.name = "valid")
  plot_df    <- merge(rho_melt, valid_melt, by = c("Var1", "Var2"))
  sig_df <- results %>%
    filter(sapply(strsplit(pair, " ~ "), `[`, 1) %in% sig_vars &
             sapply(strsplit(pair, " ~ "), `[`, 2) %in% sig_vars) %>%
    mutate(Var1 = sapply(strsplit(pair, " ~ "), `[`, 1),
           Var2 = sapply(strsplit(pair, " ~ "), `[`, 2),
           p_use = .data[[p_col]],
           stars = ifelse(p_use < 0.001, "***", ifelse(p_use < 0.01, "**", ifelse(p_use < 0.05, "*", "")))) %>%
    dplyr::select(Var1, Var2, stars)
  sig_all <- rbind(sig_df, data.frame(Var1 = sig_df$Var2, Var2 = sig_df$Var1, stars = sig_df$stars),
                   data.frame(Var1 = sig_vars, Var2 = sig_vars, stars = ""))
  plot_df <- merge(plot_df, sig_all, by = c("Var1", "Var2"), all.x = TRUE)
  plot_df$stars[is.na(plot_df$stars)] <- ""
  plot_df$label <- ifelse(plot_df$Var1 == plot_df$Var2, "", paste0(round(plot_df$rho, 2), plot_df$stars))
  plot_df$border_col <- ifelse(plot_df$Var1 == plot_df$Var2, "grey80",
                               ifelse(plot_df$valid, "#2E7D32", "#C62828"))
  p_label <- ifelse(use_adjusted_p & "spearman_p_adj" %in% names(results), "FDR-adjusted p", "raw p")
  ggplot(plot_df, aes(x = Var2, y = Var1)) +
    geom_tile(aes(fill = rho), color = "white", linewidth = 0.6) +
    geom_tile(aes(color = border_col), fill = NA, linewidth = 1.1) +
    geom_text(aes(label = label), size = 3, fontface = "bold", color = "black") +
    scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#B71C1C",
                         midpoint = 0, limits = rho_limits, name = "Spearman ρ") +
    scale_color_identity() + scale_x_discrete(position = "top") +
    labs(title = "Pairwise Spearman correlation heatmap", x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, face = "bold"),
          axis.text.y = element_text(face = "bold"), panel.grid = element_blank(),
          plot.title = element_text(face = "bold", size = 14), legend.position = "right")
}

combine_strip_plots <- function(plot_list, group_labels, results_list, strip_var,
                                title = "Spearman correlations with GDF15 pg/mL",
                                label_fill = facebd_groups,
                                label_color = "white", rel_widths = NULL,
                                label_text_size = 5, legend_width = 1) {
  all_rhos <- lapply(results_list, function(res) {
    res %>% filter(sapply(strsplit(pair, " ~ "), `[`, 1) == strip_var |
                     sapply(strsplit(pair, " ~ "), `[`, 2) == strip_var) %>% pull(spearman_rho)
  }) %>% unlist()
  global_max <- max(abs(all_rhos), na.rm = TRUE)
  rho_limits <- c(-global_max, global_max)
  shared_legend <- cowplot::get_legend(
    plot_list[[1]] +
      scale_fill_gradient2(
        low      = facebd_gradient$low,   # "#2D3B6E"
        mid      = facebd_gradient$mid,   # "white"
        high     = facebd_gradient$high,  # "#6B1A3A"
        midpoint = 0,
        limits   = rho_limits,
        name     = "Spearman ρ"
      ) +
      theme(legend.position = "right", legend.title = element_text(size = 14, face = "bold"),
            legend.text = element_text(size = 12), legend.key.size = unit(0.6, "cm")))
  plots_with_labels <- lapply(seq_len(length(plot_list)), function(i) {
    fill_col <- label_fill[((i - 1) %% length(label_fill)) + 1]
    label_banner <- ggplot() +
      annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = fill_col) +
      annotate("text", x = 0.5, y = 0.5, label = group_labels[i], color = label_color,
               fontface = "bold", size = label_text_size) +
      theme_void() + theme(plot.margin = margin(0, 2, 0, 2))
    p_clean <- plot_list[[i]] +
      scale_fill_gradient2(
        low = facebd_gradient$low, mid = facebd_gradient$mid,
        high = facebd_gradient$high, midpoint = 0,
        limits = rho_limits, name = "Spearman ρ") +
      labs(title = NULL, subtitle = NULL) +  # ← doit être après scale_fill
      theme(legend.position = "none",
            plot.margin     = margin(0, 2, 5, 2),
            plot.title      = element_blank(),   # ← forcer en plus
            plot.subtitle   = element_blank())   # ← forcer en plus
    label_banner / p_clean + plot_layout(heights = c(1, 20))
  })
  combined_with_legend <- cowplot::plot_grid(
    wrap_plots(plots_with_labels, nrow = 1, widths = rel_widths), shared_legend,
    nrow = 1, rel_widths = c(sum(rel_widths), legend_width))
  cowplot::plot_grid(
    cowplot::plot_grid(
      cowplot::ggdraw() + cowplot::draw_label(title, fontface = "bold", size = 18, x = 0.01, hjust = 0),
      cowplot::ggdraw() + cowplot::draw_label(
        "Stars based on FDR-adjusted p: * p<.05  ** p<.01  *** p<.001 | Biological variables in pg/mL",
        size = 16, color = "grey40", x = 0.01, hjust = 0),
      ncol = 1, rel_heights = c(1, 0.6)),
    combined_with_legend, ncol = 1, rel_heights = c(0.08, 1))
}

# ════════════════════════════════════════════════════════
# SPEARMAN CORRELATIONS
# ════════════════════════════════════════════════════════

# ── Subgroup dataframes ───────────────────────────────────
df_infla <- df_trans[c(names(df_trans)[grepl("^IL_", names(df_trans))],
                       "TNF_beta", "TNF_alpha", "TNF_R1", "TNF_R2",
                       "IFN_gamma", "ICAM_1", "VCAM_1", "crp_lbstresc", "GDF15 pg/ml")]
names(df_infla) <- c("IL-12/IL-23p40", "IL-15", "IL-16", "IL-17A", "IL-5", "IL-7",
                     "IL-10", "IL-12p70", "IL-13", "IL-6", "IL-8", "IL-1RA", "IL-6R",
                     "TNF-β", "TNF-α", "TNF-R1", "TNF-R2", "IFN-γ", "ICAM-1", "VCAM-1", "CRP", "GDF15")

df_metabo <- df_trans[c("trig_lbstresc", "gluc_lbstresc", "chol_lbstresc", "hdl_lbstresc",
                        "ldl_lbstresc", "mtCN", "CCF_MTDNA_ND1", "creat_lbstresc", "Lactate", "GDF15 pg/ml")]
names(df_metabo) <- c("Triglycerides", "Glycemia", "Total Cholesterol", "HDL Cholesterol",
                      "LDL Cholesterol", "mtDNA Copy Number", "ccf-mtDNA",
                      "Creatinine", "Lactate", "GDF15")

df_kyn <- df_trans[c("TRP", "KYN", "OHKYN", "KA", "QUINA", "XA", "AA", "QUINO", "PICO", "GDF15 pg/ml")]
names(df_kyn) <- c("Tryptophan", "Kynurenine", "3-Hydroxykynurenine", "Kynurenic acid",
                   "Quinaldic acid", "Xanthurenic acid", "Anthranilic acid",
                   "Quinolinic acid", "Picolinic acid", "GDF15")

df_blood <- df_trans[c("rbc_lbstresc", "plat_lbstresc", "mono_lbstresc", "lym_lbstresc",
                       "eos_lbstresc", "baso_lbstresc", "neut_lbstresc", "wbc_lbstresc", "GDF15 pg/ml")]
names(df_blood) <- c("Red Blood Cells", "Platelets", "Monocytes", "Lymphocytes",
                     "Eosinophils", "Basophils", "Neutrophils", "White Blood Cells", "GDF15")

df_clinical <- df_trans[, c("bmi", "age", "madrs_", "ymrs_num", "fagers", "fast_", "bis10",
                             "staya", "mars_", "mathys_", "psqi_", "als_", "ctq39", "qidsr120", "GDF15 pg/ml")]
names(df_clinical) <- c("BMI", "Age", "MADRS", "YMRS", "FAGERS", "FAST", "BIS",
                        "STAY-A", "MARS", "MATHYS", "PSQI", "ALS", "CTQ", "QIDS", "GDF15")

# ── Run pairwise Spearman ─────────────────────────────────
run_and_plot_spearman <- function(df_sub, strip_var = "GDF15",
                                  filename, width = 8, height = 7) {
  results <- test_spearman_pairwise(df = df_sub, alpha = 0.05,
                                    check_autocor = FALSE, maha_threshold = 0.99)
  p <- plot_spearman_heatmap(df = df_sub, results = results, strip_plot = TRUE,
                             strip_var = strip_var, use_adjusted_p = TRUE, alpha = 0.05)
  ggsave(filename, plot = p, width = width, height = height, dpi = 300, bg = "white")
  return(list(results = results, plot = p))
}

out_clinical <- run_and_plot_spearman(df_clinical, filename = "outcome/stripplot_clinical_GDF15.png")
out_blood    <- run_and_plot_spearman(df_blood,    filename = "outcome/stripplot_blood_cells_GDF15.png")
out_infla    <- run_and_plot_spearman(df_infla,    filename = "outcome/stripplot_infla_GDF15.png")
out_kyn      <- run_and_plot_spearman(df_kyn,      filename = "outcome/stripplot_kyn_GDF15.png")
out_metabo   <- run_and_plot_spearman(df_metabo,   filename = "outcome/stripplot_metabo_GDF15.png")

# ── All variables heatmap ─────────────────────────────────
results_spearman <- test_spearman_pairwise(df = df_trans, alpha = 0.05,
                                           check_autocor = FALSE, maha_threshold = 0.99,
                                           id_col = "fondacode")

plot_spearman_heatmap(df = df_trans, results = results_spearman,
                             alpha = 0.05, focus_var = "GDF15 pg/ml")

# ── Combined strip plot ───────────────────────────────────
ggsave("outcome/combined_stripplot_color_GDF15.png",
       combine_strip_plots(
         plot_list    = list(out_clinical$plot, out_blood$plot, out_infla$plot,
                             out_kyn$plot, out_metabo$plot),
         results_list = list(out_clinical$results, out_blood$results, out_infla$results,
                             out_kyn$results, out_metabo$results),
         strip_var    = "GDF15",
         group_labels = c("Clinical", "Blood cells", "Inflammation", "KYN pathway", "Metabolism"),
         rel_widths   = c(1, 1, 1, 1, 1)),
       width = 22, height = 12, dpi = 300, bg = "white")