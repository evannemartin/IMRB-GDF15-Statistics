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
df_mito <- read_excel("data/MitoBDPRS.xlsx")

# Add MitoPRS
df_mito <- df_mito %>%
  mutate(ID_padded = stringr::str_pad(as.character(ID), width = 10, side = "left", pad = "0"))

df_trans <- df_trans %>%
  left_join(df_mito, by = c("fondacode" = "ID_padded"))

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
                                label_text_size = 5.5, legend_width = 1) {
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
        "Stars based on FDR-adjusted p: * p<.05  ** p<.01  *** p<.001",
        size = 18, color = "grey40", x = 0.01, hjust = 0),
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
names(df_kyn) <- c("Tryptophan", "Kynurenine", "3-HK", "Kynurenic acid",
                   "Quinaldic acid", "Xanthurenic acid", "Anthranilic acid",
                   "Quinolinic acid", "Picolinic acid", "GDF15")

df_blood <- df_trans[c("rbc_lbstresc", "plat_lbstresc", "mono_lbstresc", "lym_lbstresc",
                       "eos_lbstresc", "baso_lbstresc", "neut_lbstresc", "wbc_lbstresc", "GDF15 pg/ml")]
names(df_blood) <- c("Red Blood Cells", "Platelets", "Monocytes", "Lymphocytes",
                     "Eosinophils", "Basophils", "Neutrophils", "White Blood Cells", "GDF15")

df_clinical <- df_trans[, c("bmi", "age", "madrs_", "ymrs_num", "fagers", "fast_", "bis10",
                             "staya", "mars_", "mathys_", "psqi_", "als_", "ctq39", "GDF15 pg/ml")]
names(df_clinical) <- c("BMI", "Age", "MADRS", "YMRS", "FAGERS", "FAST", "BIS",
                        "STAI", "MARS", "MATHYS", "PSQI", "ALS", "CTQ", "GDF15")

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

###################### LINEAR REGRESSION ######################


# ════════════════════════════════════════════════════════
# FUNCTIONS
# ════════════════════════════════════════════════════════

run_robust_lm <- function(df, outcome = "log_GDF15") {
  predictors <- setdiff(names(df)[sapply(df, is.numeric)], outcome)
  results <- list()
  for (pred in predictors) {
    x <- df[[pred]]; y <- df[[outcome]]
    valid <- complete.cases(x, y); x <- x[valid]; y <- y[valid]; n <- length(x)
    mod_rob <- rlm(y ~ x, method = "MM"); mod_ols <- lm(y ~ x); s <- summary(mod_rob)
    beta <- coef(mod_rob)[["x"]]; se <- s$coefficients["x", "Std. Error"]
    p_val <- 2 * pt(abs(beta / se), df = n - 2, lower.tail = FALSE)
    r2 <- 1 - sum(mod_rob$residuals^2) / sum((y - mean(y))^2)
    results[[pred]] <- data.frame(
      predictor = pred, n = n, beta_robust = round(beta, 4),
      beta_ols = round(coef(mod_ols)[["x"]], 4), se = round(se, 4),
      p_value = round(p_val, 4), r_squared = round(r2, 4),
      n_downweighted = sum(mod_rob$w < 0.5),
      significance = ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**",
                                                         ifelse(p_val < 0.05, "*", "ns"))))
  }
  do.call(rbind, results) %>% arrange(p_value)
}

run_anova_against_outcome <- function(df, outcome = "log_GDF15", alpha = 0.05) {
  predictors <- df %>% dplyr::select(where(~ is.factor(.) | is.character(.))) %>% names()
  results <- list()
  for (pred in predictors) {
    x <- as.factor(df[[pred]]); y <- df[[outcome]]
    valid <- complete.cases(x, y); x <- droplevels(x[valid]); y <- y[valid]
    n <- length(x); n_groups <- nlevels(x)
    if (n_groups < 2 || any(table(x) < 3) ||
        any(tapply(y, x, var, na.rm = TRUE) == 0, na.rm = TRUE)) next
    mod <- aov(y ~ x); s <- summary(mod)
    f_stat <- s[[1]]["x", "F value"]; p_val <- s[[1]]["x", "Pr(>F)"]
    eta2 <- s[[1]]["x", "Sum Sq"] / sum(s[[1]][, "Sum Sq"])
    lm_coefs <- coef(lm(y ~ x))[-1]
    estimates_str <- paste0("ref=", levels(x)[1], " | ",
                            paste(paste0(gsub("^x", "", names(lm_coefs)), ": ",
                                         ifelse(lm_coefs > 0, paste0("+", round(lm_coefs, 3)), round(lm_coefs, 3))),
                                  collapse = " | "))
    sw  <- tryCatch(shapiro.test(residuals(mod)), error = function(e) list(p.value = NA))
    lev <- tryCatch(leveneTest(y ~ x), error = function(e) data.frame(`Pr(>F)` = NA))
    sw_p <- ifelse(is.na(sw$p.value), 0, sw$p.value); lev_p <- lev$`Pr(>F)`[1]
    assumptions_ok <- sw_p > alpha & !is.na(lev_p) & lev_p > alpha
    if (assumptions_ok) {
      posthoc_method <- ifelse(n_groups == 2, "t-test", "Tukey HSD")
      if (n_groups == 2) p_val <- t.test(y ~ x)$p.value
    } else {
      posthoc_method <- ifelse(n_groups == 2, "Wilcoxon", "Kruskal-Wallis")
      p_val <- if (n_groups == 2) wilcox.test(y ~ x)$p.value else kruskal.test(y ~ x)$p.value
    }
    results[[pred]] <- data.frame(
      predictor = pred, n = n, n_groups = n_groups, estimates = estimates_str,
      f_statistic = round(f_stat, 4), p_value = round(p_val, 4), eta_squared = round(eta2, 4),
      shapiro_resid_p = round(sw_p, 4), levene_p = round(lev_p, 4),
      assumptions_ok = assumptions_ok, posthoc_method = posthoc_method, posthoc_sig = p_val < alpha,
      significance = ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**",
                                                         ifelse(p_val < 0.05, "*", "ns"))))
  }
  do.call(rbind, results) %>% arrange(p_value)
}

run_ancova <- function(df, outcome = "log_GDF15",
                       covariates = "age", alpha = 0.05) {
  
  all_vars   <- setdiff(names(df), c(outcome, covariates))
  cat_preds  <- all_vars[sapply(df[all_vars], function(v) is.factor(v) | is.character(v))]
  cont_preds <- setdiff(all_vars[sapply(df[all_vars], is.numeric)],
                        c(toupper(covariates),
                          paste0(toupper(substring(covariates, 1, 1)),
                                 substring(covariates, 2))))
  
  cov_formula_str <- paste(covariates, collapse = " + ")
  results <- list()
  
  for (pred in c(cat_preds, cont_preds)) {
    
    is_cat   <- is.factor(df[[pred]]) | is.character(df[[pred]])
    x        <- if (is_cat) as.factor(df[[pred]]) else df[[pred]]
    y        <- df[[outcome]]
    cov_data <- df[, covariates, drop = FALSE]
    
    valid    <- complete.cases(x, y, cov_data)
    x        <- if (is_cat) droplevels(x[valid]) else x[valid]
    y        <- y[valid]; cov_data <- cov_data[valid, , drop = FALSE]
    n        <- length(x); n_groups <- if (is_cat) nlevels(x) else NA
    
    # ── Skip checks ───────────────────────────────────────
    if (is_cat) {
      if (n_groups < 2) { message("Skipping ", pred, " — < 2 groups"); next }
      if (any(table(x) < 3)) { message("Skipping ", pred, " — group with < 3 obs"); next }
      if (any(tapply(y, x, var, na.rm = TRUE) == 0, na.rm = TRUE)) {
        message("Skipping ", pred, " — zero variance group"); next }
    } else {
      if (n < 8 || sd(x) == 0) {
        message("Skipping ", pred, " — n too small or zero variance"); next }
    }
    
    df_mod  <- data.frame(y = y, x = x, cov_data)
    
    # ── Fit models ────────────────────────────────────────
    mod <- tryCatch(
      lm(as.formula(paste("y ~ x +", cov_formula_str)), data = df_mod),
      error = function(e) { message("Model failed for ", pred); NULL })
    if (is.null(mod)) next
    
    mod_cov <- lm(as.formula(paste("y ~", cov_formula_str)), data = df_mod)
    s       <- summary(mod)
    
    coef_rows <- rownames(coef(s))[startsWith(rownames(coef(s)), "x")]
    if (length(coef_rows) == 0 || all(is.na(coef(s)[coef_rows, "Estimate"]))) {
      message("Skipping ", pred, " — coefficient aliased or NA"); next
    }
    
    # ── Effect size ───────────────────────────────────────
    eta  <- tryCatch(effectsize::eta_squared(mod, partial = TRUE), error = function(e) NULL)
    eta2 <- tryCatch({
      if (!is.null(eta)) {
        val <- eta$Eta2_partial[eta$Parameter == "x"]
        if (length(val) == 0 || !is.numeric(val)) NA else as.numeric(val)
      } else NA
    }, error = function(e) NA)
    
    if (is_cat) {
      
      # ── p-value via Anova type III — same as old function ─
      aov3   <- tryCatch(car::Anova(mod, type = "III"), error = function(e) NULL)
      f_stat <- if (!is.null(aov3) && "x" %in% rownames(aov3)) aov3["x", "F value"] else NA
      p_val  <- if (!is.null(aov3) && "x" %in% rownames(aov3)) aov3["x", "Pr(>F)"]  else NA
      
      # Fallback to standard F if Anova type III fails
      if (is.na(p_val)) {
        aov_s  <- tryCatch(summary(aov(as.formula(paste("y ~ x +", cov_formula_str)), data = df_mod)),
                           error = function(e) NULL)
        f_stat <- if (!is.null(aov_s)) aov_s[[1]]["x", "F value"] else NA
        p_val  <- if (!is.null(aov_s)) aov_s[[1]]["x", "Pr(>F)"]  else NA
      }
      
      # ── Assumption checks ─────────────────────────────────
      sw  <- tryCatch(shapiro.test(residuals(mod)), error = function(e) list(p.value = NA))
      bp  <- tryCatch(bptest(mod),                  error = function(e) list(p.value = NA))
      lev <- tryCatch(leveneTest(y ~ x),             error = function(e) data.frame(`Pr(>F)` = NA))
      sw_p <- ifelse(is.na(sw$p.value), 0, sw$p.value)
      bp_p <- ifelse(is.na(bp$p.value), 0, bp$p.value)
      lev_p <- lev$`Pr(>F)`[1]
      assumptions_ok <- sw_p > alpha & bp_p > alpha & !is.na(lev_p) & lev_p > alpha
      
      # ── Slope homogeneity ─────────────────────────────────
      cont_covs  <- covariates[sapply(covariates, function(v) is.numeric(df_mod[[v]]))]
      if (length(cont_covs) > 0) {
        other_covs  <- setdiff(covariates, cont_covs[1])
        int_formula <- if (length(other_covs) == 0) {
          as.formula(paste("y ~ x *", cont_covs[1]))
        } else {
          as.formula(paste("y ~ x *", cont_covs[1], "+",
                           paste(other_covs, collapse = " + ")))
        }
        mod_int   <- tryCatch(lm(int_formula, data = df_mod), error = function(e) NULL)
        p_slopes  <- if (!is.null(mod_int)) anova(mod, mod_int)$`Pr(>F)`[2] else NA
        slopes_ok <- !is.na(p_slopes) & p_slopes > alpha
      } else {
        p_slopes <- NA; slopes_ok <- NA
      }
      
      # ── LM estimates ──────────────────────────────────────
      lm_coefs <- coef(mod)[coef_rows]
      cov_ests <- sapply(covariates, function(v) {
        est <- round(coef(mod)[v], 4)
        paste0(v, ": ", ifelse(est > 0, paste0("+", est), est))
      })
      estimates_str <- paste0(
        paste(cov_ests, collapse = " | "), " | ref=", levels(x)[1], " | ",
        paste(paste0(gsub("^x", "", names(lm_coefs)), ": ",
                     ifelse(lm_coefs > 0, paste0("+", round(lm_coefs, 3)),
                            round(lm_coefs, 3))), collapse = " | "))
      
      results[[pred]] <- data.frame(
        predictor            = pred, type = "categorical",
        covariates           = paste(covariates, collapse = " + "),
        n = n, n_groups = n_groups, estimates = estimates_str,
        f_statistic          = round(f_stat, 4),
        p_value              = round(p_val, 4),
        partial_eta2         = ifelse(is.numeric(eta2) & length(eta2) == 1, round(eta2, 4), NA),
        shapiro_p            = round(sw_p, 4),
        breusch_pagan_p      = round(bp_p, 4),
        levene_p             = round(lev_p, 4),
        slopes_homogeneity_p = round(p_slopes, 4),
        slopes_ok            = slopes_ok,
        assumptions_ok       = assumptions_ok,
        r2_full              = round(s$r.squared, 4),
        r2_increment         = round(s$r.squared - summary(mod_cov)$r.squared, 4),
        significance         = ifelse(p_val < 0.001, "***",
                                      ifelse(p_val < 0.01,  "**",
                                             ifelse(p_val < 0.05,  "*", "ns")))
      )
      
    } else {
      
      # ── Continuous predictor ──────────────────────────────
      sw   <- tryCatch(shapiro.test(residuals(mod)), error = function(e) list(p.value = NA))
      bp   <- tryCatch(bptest(mod),                  error = function(e) list(p.value = NA))
      sw_p <- ifelse(is.na(sw$p.value), 0, sw$p.value)
      bp_p <- ifelse(is.na(bp$p.value), 0, bp$p.value)
      p_val <- coef(s)["x", "Pr(>|t|)"]
      
      results[[pred]] <- data.frame(
        predictor            = pred, type = "continuous",
        covariates           = paste(covariates, collapse = " + "),
        n = n, n_groups = NA,
        estimates            = paste0("β = ", round(coef(mod)["x"], 4)),
        f_statistic          = NA,
        p_value              = round(p_val, 4),
        partial_eta2         = ifelse(is.numeric(eta2) & length(eta2) == 1, round(eta2, 4), NA),
        shapiro_p            = round(sw_p, 4),
        breusch_pagan_p      = round(bp_p, 4),
        levene_p             = NA,
        slopes_homogeneity_p = NA,
        slopes_ok            = NA,
        assumptions_ok       = sw_p > alpha & bp_p > alpha,
        r2_full              = round(s$r.squared, 4),
        r2_increment         = round(s$r.squared - summary(mod_cov)$r.squared, 4),
        significance         = ifelse(p_val < 0.001, "***",
                                      ifelse(p_val < 0.01,  "**",
                                             ifelse(p_val < 0.05,  "*", "ns")))
      )
    }
  }
  
  df_results <- do.call(rbind, results) %>% arrange(p_value)
  
  # ── FDR correction ────────────────────────────────────────
  valid_p <- !is.na(df_results$p_value)
  df_results$p_value_adj <- NA
  df_results$p_value_adj[valid_p] <- p.adjust(
    df_results$p_value[valid_p], method = "BH")
  
  df_results$significance_adj <- ifelse(
    is.na(df_results$p_value_adj), NA,
    ifelse(df_results$p_value_adj < 0.001, "***",
           ifelse(df_results$p_value_adj < 0.01,  "**",
                  ifelse(df_results$p_value_adj < 0.05,  "*", "ns"))))
  
  return(df_results)
}


extract_robust_pvalues <- function(mod) {
  s <- summary(mod); coefs <- s$coefficients
  n <- nrow(mod$model); k <- length(coef(mod))
  t_vals <- coefs[, "t value"]
  p_vals <- 2 * pt(abs(t_vals), df = n - k, lower.tail = FALSE)
  data.frame(term = rownames(coefs), estimate = round(coefs[, "Value"], 4),
             std_error = round(coefs[, "Std. Error"], 4), t_value = round(t_vals, 4),
             p_value = round(p_vals, 4),
             sig = ifelse(p_vals < 0.001, "***", ifelse(p_vals < 0.01, "**",
                                                        ifelse(p_vals < 0.05, "*", ifelse(p_vals < 0.1, ".", "ns"))))) %>%
    `rownames<-`(NULL)
}

plot_ancova_results <- function(results, df, outcome = "log_GDF15",
                                forced_covariates = "age",
                                alpha = 0.05,
                                boxplot = FALSE,
                                force_predictors = NULL,
                                drop_levels = NULL) {
  library(patchwork)
  
  # ── Filter results ────────────────────────────────────────
  if (!is.null(force_predictors)) {
    results_sig <- results %>% filter(predictor %in% force_predictors)
  } else {
    results_sig <- results %>% filter(p_value < alpha)
  }
  
  if (nrow(results_sig) == 0) {
    cat("No significant results to plot.\n"); return(invisible(NULL))
  }
  
  cont_covar <- forced_covariates[sapply(forced_covariates,
                                         function(v) is.numeric(df[[v]]))]
  plots <- list()
  
  for (i in seq_len(nrow(results_sig))) {
    
    pred    <- results_sig$predictor[i]
    p_val   <- results_sig$p_value[i]
    p_adj   <- if ("p_value_adj" %in% names(results_sig))
      results_sig$p_value_adj[i] else NA
    sig     <- results_sig$significance[i]
    sig_adj <- if ("significance_adj" %in% names(results_sig))
      results_sig$significance_adj[i] else NA
    eta2    <- results_sig$partial_eta2[i]
    est     <- results_sig$estimates[i]
    is_cat  <- results_sig$type[i] == "categorical"
    
    x        <- df[[pred]]
    y        <- df[[outcome]]
    cov_data <- df[, forced_covariates, drop = FALSE]
    
    valid    <- complete.cases(x, y, cov_data)
    x        <- if (is_cat) droplevels(as.factor(x[valid])) else x[valid]
    y        <- y[valid]
    cov_data <- cov_data[valid, , drop = FALSE]
    cov_x    <- cov_data[[cont_covar[1]]]
    
    if (is_cat && pred %in% names(drop_levels)) {
      levels_to_drop <- drop_levels[[pred]]
      keep     <- !x %in% levels_to_drop
      x        <- droplevels(x[keep])
      y        <- y[keep]
      cov_data <- cov_data[keep, , drop = FALSE]
      cov_x    <- cov_data[[cont_covar[1]]]
    }
    
    df_plot  <- data.frame(x = x, y = y, cov_x = cov_x, cov_data)
    n_groups <- if (is_cat) nlevels(x) else NA
    
    if (is_cat) {
      n_per_group <- table(x)
      # Rename x-axis labels to include n per group
      x_labels <- setNames(
        paste0(names(n_per_group), "\n(n=", n_per_group, ")"),
        names(n_per_group)
      )
    } else {
      n_total     <- length(y)
    }
    
    # ── Extract covariate betas ───────────────────────────
    cov_betas_str <- {
      parts     <- strsplit(est, " \\| ")[[1]]
      cov_parts <- parts[sapply(forced_covariates, function(cv)
        startsWith(parts, paste0(cv, ":")))]
      if (length(cov_parts) > 0) {
        cov_parts_rounded <- sapply(cov_parts, function(cp) {
          colon_split <- strsplit(cp, ": ")[[1]]
          name_part   <- colon_split[1]
          num         <- as.numeric(colon_split[2])
          rounded     <- round(num, 3)
          formatted   <- formatC(abs(rounded), format = "f", digits = 3)
          paste0(name_part, " = ", ifelse(rounded >= 0, "+", "-"), formatted)
        })
        paste(cov_parts_rounded, collapse = "\n")
      } else ""
    }
    
    beta_display <- if (is_cat) {
      parts   <- strsplit(est, " \\| ")[[1]]
      non_ref <- parts[!startsWith(parts, "ref=") &
                         !grepl(paste(forced_covariates, collapse = "|"), parts)]
      if (pred %in% names(drop_levels)) {
        levels_to_drop <- drop_levels[[pred]]
        non_ref <- non_ref[!grepl(paste(levels_to_drop, collapse = "|"), non_ref)]
      }
      paste(paste0(pred, " = ", non_ref), collapse = "\n")
    } else {
      paste0(pred, " = ", est)
    }
    
    p_label     <- ifelse(p_val < 0.001, "p < 0.001",
                          paste0("p = ", round(p_val, 3)))
    p_adj_label <- if (!is.na(p_adj))
      ifelse(p_adj < 0.001, "p(FDR) < 0.001",
             paste0("p(FDR) = ", round(p_adj, 3))) else ""
    
    sig_str     <- paste0(" ", sig)
    sig_adj_str <- if (!is.na(p_adj)) paste0(" ", sig_adj) else ""
    
    annot_label <- paste0(
      if (nchar(cov_betas_str) > 0) paste0(cov_betas_str, "\n") else "",
      beta_display,
      "\n", p_label, sig_str,
      if (p_adj_label != "") paste0("\n", p_adj_label, sig_adj_str) else ""
    )
    
    group_colors <- if (is_cat) facebd_groups[seq_len(n_groups)] else NULL
    
    if (boxplot & is_cat) {
      
      cov_formula   <- as.formula(paste("y ~", paste(forced_covariates, collapse = " + ")))
      mod_cov       <- lm(cov_formula, data = df_plot)
      df_plot$y_adj <- residuals(mod_cov) + mean(y)
      
      p <- ggplot(df_plot, aes(x = x, y = y_adj, fill = x)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5) +
        geom_jitter(aes(color = x), width = 0.15, size = 2, alpha = 0.6) +
        scale_fill_manual(values  = group_colors, guide = "none") +
        scale_color_manual(values = group_colors, guide = "none") +
        # ── N per group shown directly on x-axis tick labels
        scale_x_discrete(labels = x_labels) +
        annotate("text", x = -Inf, y = Inf, label = annot_label,
                 hjust = 0, vjust = 1.3, size = 3.5,
                 fontface = "italic", color = "grey20") +
        labs(
          title    = paste0("Effect of ", pred),
          subtitle = paste0("Adjusted for ", paste(forced_covariates, collapse = " + ")),
          x = pred,
          y = outcome
        ) +
        theme_minimal(base_size = 11) +
        theme(plot.title       = element_text(face = "bold", size = 12),
              plot.subtitle    = element_text(size = 10, color = "grey40"),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size = 12),
              axis.text.x      = element_text(angle = 20, hjust = 1, size = 11))
      
    } else if (!is_cat) {
      
      p <- ggplot(df_plot, aes(x = x, y = y)) +
        geom_point(alpha = 0.6, size = 2.5, color = "grey40") +
        geom_smooth(method = "lm", se = TRUE,
                    color = "#1565C0", fill = "#1565C0", alpha = 0.15) +
        annotate("text", x = -Inf, y = Inf, label = annot_label,
                 hjust = 0, vjust = 1.3, size = 3.2,
                 fontface = "italic", family = "mono", color = "grey20") +
        labs(
          title = paste0(pred, " → ", outcome, " | adjusted for ",
                         paste(forced_covariates, collapse = " + ")),
          x = pred,
          y = outcome
        ) +
        theme_minimal(base_size = 11) +
        theme(plot.title       = element_text(face = "bold", size = 12),
              panel.grid.minor = element_blank())
      
    } else {
      
      p <- ggplot(df_plot, aes(x = cov_x, y = y, color = x, fill = x)) +
        geom_point(alpha = 0.6, size = 2.5) +
        geom_smooth(method = "lm", se = TRUE, alpha = 0.15) +
        scale_color_manual(values = group_colors, name = pred) +
        scale_fill_manual(values  = group_colors, name = pred) +
        # ── N per group shown in legend labels
        scale_color_manual(values = group_colors, name = pred,
                           labels = x_labels) +
        scale_fill_manual(values  = group_colors, name = pred,
                          labels = x_labels) +
        annotate("text", x = -Inf, y = Inf, label = annot_label,
                 hjust = 0, vjust = 1.3, size = 3.2,
                 fontface = "italic", family = "mono", color = "grey20") +
        labs(
          title = paste0(pred, " → ", outcome, " | adjusted for ",
                         paste(forced_covariates, collapse = " + ")),
          x = cont_covar[1],
          y = outcome
        ) +
        theme_minimal(base_size = 11) +
        theme(plot.title       = element_text(face = "bold", size = 12),
              panel.grid.minor = element_blank(),
              legend.position  = "top")
    }
    
    plots[[pred]] <- p
  }
  
  # ── Shared Y scale across all panels ─────────────────────
  wrap_plots(plots, ncol = 2) & scale_y_continuous(limits = range(df[[outcome]], na.rm = TRUE))
}

to_superscript <- function(n) {
  supers <- c("\u00b9", "\u00b2", "\u00b3", "\u2074", "\u2075",
              "\u2076", "\u2077", "\u2078", "\u2079")
  supers[n]
}

build_unit_labels <- function(vars, units_map) {
  relevant    <- units_map[names(units_map) %in% vars]
  unique_units <- unique(unname(relevant))
  idx          <- seq_along(unique_units)
  unit_index   <- setNames(idx, unique_units)
  
  label_map <- setNames(vars, vars)
  for (v in names(relevant)) {
    u   <- relevant[[v]]
    sup <- to_superscript(unit_index[[u]])
    label_map[[v]] <- paste0(v, sup)
  }
  
  legend_lines <- paste0(to_superscript(idx), " ", unique_units)
  list(label_map = label_map, legend_lines = legend_lines)
}


export_results_table <- function(results, filename, title = "ANCOVA Results",
                                 alpha = 0.05,
                                 units_map = NULL) {         
  library(gt)
  
  is_grouped <- "group" %in% names(results)
  
  results_sig <- results %>%
    filter(p_value < alpha) %>%
    mutate(
      p_value      = ifelse(p_value < 0.001, "< 0.001",
                            as.character(round(p_value, 3))),
      p_value_adj  = ifelse(!is.na(p_value_adj) & as.numeric(p_value_adj) < 0.001,
                            "< 0.001",
                            as.character(round(as.numeric(p_value_adj), 3))),
      partial_eta2 = round(partial_eta2, 3),
      r2_increment = round(r2_increment, 3)
    )
  
  if (nrow(results_sig) == 0) {
    cat("No significant results to export.\n"); return(invisible(NULL))
  }
  
  if (!is.null(units_map)) {
    vars_in_table <- unique(results_sig$predictor)
    ul <- build_unit_labels(vars_in_table, units_map)
    results_sig <- results_sig %>%
      mutate(predictor = dplyr::recode(predictor, !!!ul$label_map))
    unit_footnote <- paste(ul$legend_lines, collapse = "   ")
  } else {
    unit_footnote <- NULL
  }
  # ───────────────────────────────────────────────────────────────────────────
  
  cols_to_show <- intersect(
    c("group", "predictor", "covariates", "n", "estimates",
      "p_value", "p_value_adj", "partial_eta2", "r2_increment",
      "significance", "significance_adj"),
    names(results_sig)
  )
  results_sig <- results_sig %>% dplyr::select(all_of(cols_to_show))
  
  tbl <- if (is_grouped) {
    results_sig %>% gt(groupname_col = "group")
  } else {
    results_sig %>% gt()
  }
  
  tbl <- tbl %>%
    tab_header(title    = title,
               subtitle = paste0("Significant results only (p < ", alpha, ")")) %>%
    cols_label(
      predictor        = "Predictor",         
      covariates       = "Adjusted for",
      n                = "N",
      estimates        = "Estimates",
      p_value          = "p-value",
      p_value_adj      = "p-value (FDR)",
      partial_eta2     = "Partial \u03b7\u00b2",
      r2_increment     = "\u0394R\u00b2",
      significance     = "Sig.",
      significance_adj = "Sig. (FDR)"
    ) %>%
    cols_hide(columns = any_of(c("type", "n_groups", "f_statistic", "shapiro_p",
                                 "levene_p", "breusch_pagan_p", "slopes_homogeneity_p",
                                 "slopes_ok", "assumptions_ok", "r2_full"))) %>%
    # Header colonnes
    tab_style(
      style = list(cell_fill(color = facebd_colors$dark_burgundy),
                   cell_text(color = "white", weight = "bold")),
      locations = cells_column_labels()
    ) %>%
    # Header groupes
    tab_style(
      style = list(cell_fill(color = "#F2B8CC"),
                   cell_text(color = facebd_colors$dark_burgundy, weight = "bold")),
      locations = cells_row_groups()
    ) %>%
    # Stars significatives
    tab_style(
      style = cell_text(weight = "bold", color = facebd_colors$dark_burgundy),
      locations = cells_body(columns = significance,
                             rows = significance %in% c("*", "**", "***"))
    ) %>%
    tab_style(
      style = cell_text(weight = "bold", color = facebd_colors$dark_burgundy),
      locations = cells_body(columns = significance_adj,
                             rows = significance_adj %in% c("*", "**", "***"))
    ) %>%
    tab_style(
      style = cell_text(weight = "bold", size = 16),
      locations = cells_title(groups = "title")) %>%
    opt_row_striping() %>%
    tab_options(table.font.size = 12, table.width = pct(100)) %>%
    tab_footnote(
      footnote  = "\u00b7 p < 0.10  * p < 0.05  ** p < 0.01  *** p < 0.001",
      locations = cells_column_labels(columns = significance)
    ) 
  
  if (!is.null(unit_footnote)) {
    tbl <- tbl %>%
      tab_footnote(
        footnote  = unit_footnote,
        locations = cells_column_labels(columns = predictor)
      )
  }
  # ───────────────────────────────────────────────────────────────────────────
  
  gtsave(tbl, filename = filename)
  cat("Table exported to:", filename, "\n")
  return(invisible(tbl))
}


# ════════════════════════════════════════════════════════
# REGRESSION - CATEGORICAL AND CONTINUOUS
# ════════════════════════════════════════════════════════

# ── Log transformation ──────────────────────────────────
df_trans$log_GDF15 <- log10(df_trans$`GDF15 pg/ml`)

###################### USAGE FOR GDF15 with categorical data ######################

df_clinical_cat <- df_trans[c("sex", "arm", "edulevel", "Antidepressants_treat", "Anxiolytics_treat", "Lithium_treat",
                              "Antipsychotics_treat", "Thymoregulators_treat", "cyclerap", "saison", "rad_tb_subst", "suoccur_alcool", "suoccur_cannabis", "current smokers", "remitted smokers", "Mito_PRS_group", "log_GDF15", "age")]
names(df_clinical_cat) <- c("Sex", "BD subtype", "Education level", "Antidepressants treatment", "Anxiolytics treatment", "Lithium treatment",
                            "Antipsychotics treatment", "Thymoregulators treatment", "Rapid cycling", "Season", "Substance use disorder", "Alcohol use disorder", "Cannabis use disorder", "Current smokers", "Remitted smokers", "Mito PRS group", "GDF15 [log10(pg/mL)]", "Age")


# ── Replace groups < 3 obs with NA ───────
df_clinical_cat <- df_clinical_cat %>%
  mutate(
    `Education level` = {
      x <- as.factor(`Education level`)
      group_sizes <- table(x)
      small_groups <- names(group_sizes[group_sizes < 3])
      ifelse(`Education level` %in% small_groups, NA, `Education level`)
    }
  )

results_ancova <- run_ancova(
  df        = df_clinical_cat,
  outcome   = "GDF15 [log10(pg/mL)]",
  covariate = "Age"
)
cat("\n=== Significant ANCOVA results (p < 0.05) ===\n")
print(results_ancova %>% filter(p_value < 0.05))

df_anc <- df_clinical_cat %>%
  dplyr::select(all_of(c("GDF15 [log10(pg/mL)]", "BD subtype", "Lithium treatment", "Age"))) %>%
  filter(complete.cases(.))

best_formula_clinical_ancova <- formula(step(
  lm(`GDF15 [log10(pg/mL)]` ~ `BD subtype` + `Lithium treatment`+ Age, data = df_anc),
  direction = "both"
))
print(best_formula_clinical_ancova)

mod_ancova <- lm(best_formula_clinical_ancova, data = df_clinical_cat)
summary(mod_ancova)

ggsave("outcome/boxplot_color_ancova_clinical_GDF15_ns_numbers.png",
      plot_ancova_results(results_ancova, df_clinical_cat, "GDF15 [log10(pg/mL)]", "Age", force_predictors = c("Sex", "Mito PRS group"), boxplot = TRUE, drop_levels = list(`Mito PRS group` = "Mid")),
      width = 12, height = 8, dpi = 300, bg = "white")

ggsave("outcome/scatterplot_color_ancova_clinical_GDF15_numbers.png",
       plot_ancova_results(results_ancova, df_clinical_cat, "GDF15 [log10(pg/mL)]", "Age"),
       width = 12, height = 8, dpi = 300, bg = "white")

ggsave("outcome/boxplot_color_ancova_clinical_GDF15_numbers.png",
       plot_ancova_results(results_ancova, df_clinical_cat, "GDF15 [log10(pg/mL)]", "Age", boxplot = TRUE),
       width = 12, height = 8, dpi = 300, bg = "white")



###################### USAGE FOR GDF15 with continuous data ######################


# Adjusted regression on all continuous variables

# ── Run ANCOVA on each subgroup for both adjustments ──────
df_list <- list(
  "Blood cells"  = df_blood,
  "Inflammation" = df_infla,
  "KYN pathway"  = df_kyn,
  "Metabolism"   = df_metabo,
  "Clinical"     = df_clinical
)

rename_map <- c(
  "age"           = "Age",
  "Lithium_treat" = "Lithium treatment",
  "log_GDF15"     = "Log(GDF15)"
)

for (covs in list(c("age"), c("age", "Lithium_treat"))) {
  
  cov_label <- paste(covs, collapse = " + ")
  cat("\n══════════════════════════════════════════\n")
  cat("Running ANCOVA adjusted for:", cov_label, "\n")
  cat("══════════════════════════════════════════\n")
  
  results_by_group <- do.call(rbind, lapply(names(df_list), function(group_name) {
    df_sub <- df_list[[group_name]] %>%
      dplyr::select(-any_of(c("GDF15", "GDF15 pg/ml"))) 
    df_sub$log_GDF15 <- df_trans$log_GDF15
    
    for (cov in covs) df_sub[[cov]] <- df_trans[[cov]]
    
    res <- tryCatch(
      run_ancova(df = df_sub, outcome = "log_GDF15", covariates = covs),
      error = function(e) { message("Error for ", group_name, ": ", e$message); NULL }
    )
    if (!is.null(res) && nrow(res) > 0) { res$group <- group_name; res }
  })) %>%
    filter(!predictor %in% covs) %>%
    distinct(predictor, group, .keep_all = TRUE) %>%
    # ── Apply display names ────────────────────────────────
    mutate(
      predictor  = ifelse(predictor %in% names(rename_map),
                          rename_map[predictor], predictor),
      covariates = sapply(covariates, function(cov_str) {
        for (old in names(rename_map)) {
          cov_str <- gsub(old, rename_map[old], cov_str, fixed = TRUE)
        }
        cov_str
      }),
      estimates = sapply(estimates, function(est_str) {
        for (old in names(rename_map)) {
          est_str <- gsub(old, rename_map[old], est_str, fixed = TRUE)
        }
        est_str
      })
    ) %>%
    arrange(group, p_value)
  
  cat("\n=== Significant results (p < 0.05) ===\n")
  print(results_by_group %>% filter(p_value < 0.05) %>%
          dplyr::select(group, predictor, p_value, p_value_adj, significance_adj))
  
  # ── Display names for covariates in title ────────────────
  cov_label_display <- paste(
    sapply(covs, function(c) ifelse(c %in% names(rename_map), rename_map[c], c)),
    collapse = " + "
  )
  cov_filename <- gsub(" \\+ ", "_", cov_label)
  
  group_order <- c("Clinical", "Blood cells", "Inflammation",
                   "KYN pathway", "Metabolism")
  
  results_by_group <- results_by_group %>%
    mutate(group = factor(group, levels = group_order)) %>%
    arrange(group, p_value)

  units_map <- c(
    # Blood cells
    "Neutrophils" = "10\u2079/L",
    "Basophils"   = "10\u2079/L",
    "White Blood Cells" = "10\u2079/L",
    # Inflammation
    "IL-10"       = "pg/mL",
    # KYN pathway
    "Kynurenine"  = "\u00b5M"
  )
  
  export_results_table(
    results  = results_by_group,
    filename = paste0("outcome/table_color_ancova_adj_", cov_filename, "_GDF15.png"),
    title    = paste0("ANCOVA — GDF15 [log10(pg/mL)] adjusted for ", cov_label_display),
    units_map = units_map
  )
}
