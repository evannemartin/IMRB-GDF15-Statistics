
library(readxl)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggrepel)

setwd("~/IMRB/Stats GDF15/")

# Import data
df=read_excel("data/merged_data_samples_V0.xlsx") ; head(df[,1:27])


###################### SPEARMAN CORRELATION MATRIX ######################

test_spearman_pairwise <- function(df, alpha = 0.05, check_autocor = FALSE,
                                   maha_threshold = 0.975,
                                   id_col = "fondacode") {
  
  num_vars  <- df[, sapply(df, is.numeric)]
  var_names <- names(num_vars)
  n_vars    <- length(var_names)
  
  results       <- list()
  outlier_plots <- list()
  
  for (i in 1:(n_vars - 1)) {
    for (j in (i + 1):n_vars) {
      
      x    <- num_vars[[i]]
      y    <- num_vars[[j]]
      pair <- paste0(var_names[i], " ~ ", var_names[j])
      
      # Remove NAs
      valid    <- complete.cases(x, y)
      x        <- as.numeric(x[valid])
      y        <- as.numeric(y[valid])
      df_valid <- df[valid, ]
      n        <- length(x)
      
      # Minimum checks
      if (n < 8 || sd(x) == 0 || sd(y) == 0) {
        results[[paste0(var_names[i], "_", var_names[j])]] <- data.frame(
          pair              = pair,
          n                 = n,
          n_outliers        = NA, independence_ok   = FALSE,
          spearman_rho      = NA, spearman_p        = NA,
          rho_without_out   = NA, rho_difference    = NA,
          outlier_influence = NA,
          spearman_ok       = FALSE, significance   = NA,
          advice            = "ERROR: n too small or zero variance",
          error             = "n too small or zero variance"
        )
        next
      }
      
      rx <- rank(x); ry <- rank(y)
      
      # --- Bivariate outliers (Mahalanobis on ranks) ---
      n_out            <- NA
      rho_without_out  <- NA
      rho_difference   <- NA
      outlier_influence <- NA
      mat     <- cbind(rx, ry)
      cov_mat <- cov(mat)
      
      if (det(cov_mat) >= 1e-10) {
        maha  <- mahalanobis(mat, colMeans(mat), cov_mat)
        idx   <- which(maha > qchisq(maha_threshold, df = 2))
        n_out <- length(idx)
        
        if (n_out > 0) {
          
          # Outlier labels
          if (id_col %in% names(df)) {
            outlier_ids    <- as.character(df_valid[[id_col]][idx])
            outlier_labels <- outlier_ids
          } else {
            outlier_ids    <- paste0("row:", which(valid)[idx])
            outlier_labels <- outlier_ids
            warning(paste("Column", id_col, "not found — using row index instead"))
          }
          
          # --- Rho comparison with vs without outliers ---
          rho_with       <- cor(x, y, method = "spearman")
          x_no_out       <- x[-idx]
          y_no_out       <- y[-idx]
          rho_without_out <- cor(x_no_out, y_no_out, method = "spearman")
          rho_difference  <- abs(rho_with - rho_without_out)
          
          outlier_influence <- ifelse(
            rho_difference < 0.05,  "Low — outlier(s) not influential",
            ifelse(
              rho_difference < 0.10, "Moderate — mention in methods",
              "High — investigate outlier(s)"
            )
          )
          
          # Console summary
          cat("\n──────────────────────────────────────────\n")
          cat("Pair:", pair, "\n")
          cat("Outlier(s) detected:", paste(outlier_labels, collapse = ", "), "\n")
          cat("rho with outlier(s)   :", round(rho_with, 4), "\n")
          cat("rho without outlier(s):", round(rho_without_out, 4), "\n")
          cat("difference            :", round(rho_difference, 4), "\n")
          cat("influence             :", outlier_influence, "\n")
          cat("──────────────────────────────────────────\n")
          
          # Plot
          plot_data <- data.frame(
            rx      = rx,
            ry      = ry,
            outlier = seq_along(rx) %in% idx,
            label   = ifelse(seq_along(rx) %in% idx, outlier_labels, "")
          )
          
          outlier_plots[[pair]] <- ggplot(plot_data, aes(x = rx, y = ry)) +
            geom_point(aes(color = outlier), size = 2.5) +
            geom_text_repel(aes(label = label), size = 3, color = "firebrick") +
            scale_color_manual(
              values = c("FALSE" = "steelblue", "TRUE" = "firebrick"),
              labels = c("Normal", "Outlier"),
              name   = ""
            ) +
            labs(
              title    = paste("Bivariate outliers:", pair),
              subtitle = paste0(
                n_out, " outlier(s) detected — Mahalanobis threshold ",
                maha_threshold * 100, "%\n",
                "rho with: ", round(rho_with, 3),
                " | rho without: ", round(rho_without_out, 3),
                " | diff: ", round(rho_difference, 3),
                " → ", outlier_influence
              ),
              x = paste0("rank(", var_names[i], ")"),
              y = paste0("rank(", var_names[j], ")")
            ) +
            theme_minimal(base_size = 11) +
            theme(
              plot.title    = element_text(face = "bold"),
              plot.subtitle = element_text(color = "firebrick", size = 9)
            )
        }
      } else {
        cat("Warning: singular covariance matrix for:", pair, "\n")
      }
      
      # --- Independence check (optional Durbin-Watson on ranks) ---
      independ <- if (check_autocor) {
        dwtest(lm(ry ~ rx))$p.value > alpha
      } else TRUE
      
      # --- Spearman correlation ---
      spearman    <- cor.test(x, y, method = "spearman", exact = FALSE)
      no_outlier  <- !is.na(n_out) & n_out == 0
      spearman_ok <- independ & no_outlier
      
      advice <- if (is.na(n_out)) {
        "Singular matrix — check pair manually"
      } else if (spearman_ok) {
        "Spearman valid"
      } else if (!independ) {
        "Warning: autocorrelation detected"
      } else {
        paste0("Check ", n_out, " outlier(s) — influence: ", outlier_influence)
      }
      
      results[[paste0(var_names[i], "_", var_names[j])]] <- data.frame(
        pair              = pair,
        n                 = n,
        n_outliers        = n_out,
        independence_ok   = independ,
        spearman_rho      = round(spearman$estimate, 4),
        spearman_p        = round(spearman$p.value, 4),
        rho_without_out   = round(rho_without_out, 4),
        rho_difference    = round(rho_difference, 4),
        outlier_influence = outlier_influence,
        spearman_ok       = spearman_ok,
        significance      = ifelse(spearman$p.value < 0.001, "***",
                                   ifelse(spearman$p.value < 0.01,  "**",
                                          ifelse(spearman$p.value < 0.05,  "*", "ns"))),
        advice            = advice,
        error             = NA_character_
      )
    }
  }
  
  df_results <- do.call(rbind, results) %>% arrange(spearman_p)
  
  # Display outlier plots after the loop
  if (length(outlier_plots) > 0) {
    cat("\n══════════════════════════════════════════\n")
    cat(length(outlier_plots), "pair(s) with bivariate outliers detected:\n")
    cat("══════════════════════════════════════════\n")
    for (nm in names(outlier_plots)) cat("->", nm, "\n")
    cat("\n")
    for (p in outlier_plots) print(p)
  } else {
    cat("\n✓ No bivariate outliers detected in any pair.\n")
  }
  
  return(df_results)
}



###################### SPEARMAN CORRELATION HEATMAP ######################

plot_spearman_heatmap <- function(df, results, alpha = 0.05,
                                  focus_var = NULL) {
  
  num_vars  <- df[, sapply(df, is.numeric)]
  var_names <- names(num_vars)
  
  # --- Filter significant pairs ---
  sig_pairs <- results %>% filter(spearman_p < alpha)
  
  if (nrow(sig_pairs) == 0) {
    cat("No significant pairs found at alpha =", alpha, "\n")
    return(invisible(NULL))
  }
  
  # --- If focus_var is specified, keep only pairs involving that variable ---
  if (!is.null(focus_var)) {
    
    if (!focus_var %in% c(
      sapply(strsplit(results$pair, " ~ "), `[`, 1),
      sapply(strsplit(results$pair, " ~ "), `[`, 2)
    )) {
      stop(paste("Variable", focus_var, "not found in results pairs."))
    }
    
    sig_pairs <- sig_pairs %>%
      filter(
        sapply(strsplit(pair, " ~ "), `[`, 1) == focus_var |
          sapply(strsplit(pair, " ~ "), `[`, 2) == focus_var
      )
    
    if (nrow(sig_pairs) == 0) {
      cat("No significant pairs found involving", focus_var, "at alpha =", alpha, "\n")
      return(invisible(NULL))
    }
    
    cat("Filtering pairs significantly associated with:", focus_var, "\n")
  }
  
  # --- Extract relevant variable names ---
  sig_vars <- unique(c(
    sapply(strsplit(sig_pairs$pair, " ~ "), `[`, 1),
    sapply(strsplit(sig_pairs$pair, " ~ "), `[`, 2)
  ))
  
  cat(length(sig_vars), "variables retained:\n")
  cat(paste("->", sig_vars, collapse = "\n"), "\n\n")
  
  n_sig <- length(sig_vars)
  
  # --- Build symmetric rho matrix ---
  rho_mat   <- matrix(NA,    n_sig, n_sig, dimnames = list(sig_vars, sig_vars))
  valid_mat <- matrix(FALSE, n_sig, n_sig, dimnames = list(sig_vars, sig_vars))
  diag(rho_mat)   <- 1
  diag(valid_mat) <- TRUE
  
  for (k in seq_len(nrow(results))) {
    parts <- strsplit(results$pair[k], " ~ ")[[1]]
    vi <- parts[1]; vj <- parts[2]
    if (!(vi %in% sig_vars) || !(vj %in% sig_vars)) next
    rho_mat[vi, vj]   <- results$spearman_rho[k]
    rho_mat[vj, vi]   <- results$spearman_rho[k]
    valid_mat[vi, vj] <- results$spearman_ok[k]
    valid_mat[vj, vi] <- results$spearman_ok[k]
  }
  
  # --- Melt for ggplot ---
  rho_melt   <- melt(rho_mat,   varnames = c("Var1", "Var2"), value.name = "rho")
  valid_melt <- melt(valid_mat, varnames = c("Var1", "Var2"), value.name = "valid")
  plot_df    <- merge(rho_melt, valid_melt, by = c("Var1", "Var2"))
  
  # --- Significance stars ---
  sig_df <- results %>%
    filter(
      sapply(strsplit(pair, " ~ "), `[`, 1) %in% sig_vars &
        sapply(strsplit(pair, " ~ "), `[`, 2) %in% sig_vars
    ) %>%
    mutate(
      Var1  = sapply(strsplit(pair, " ~ "), `[`, 1),
      Var2  = sapply(strsplit(pair, " ~ "), `[`, 2),
      stars = ifelse(spearman_p < 0.001, "***",
                     ifelse(spearman_p < 0.01,  "**",
                            ifelse(spearman_p < 0.05,  "*", "")))
    ) %>%
    select(Var1, Var2, stars)
  
  sig_mirror      <- sig_df
  sig_mirror$Var1 <- sig_df$Var2
  sig_mirror$Var2 <- sig_df$Var1
  diag_df <- data.frame(Var1 = sig_vars, Var2 = sig_vars, stars = "")
  sig_all <- rbind(sig_df, sig_mirror, diag_df)
  
  plot_df <- merge(plot_df, sig_all, by = c("Var1", "Var2"), all.x = TRUE)
  plot_df$stars[is.na(plot_df$stars)] <- ""
  plot_df$label <- ifelse(
    plot_df$Var1 == plot_df$Var2, "",
    paste0(round(plot_df$rho, 2), plot_df$stars)
  )
  plot_df$border_col <- ifelse(
    plot_df$Var1 == plot_df$Var2, "grey80",
    ifelse(plot_df$valid, "#2E7D32", "#C62828")
  )
  
  # --- Subtitle ---
  subtitle <- if (!is.null(focus_var)) {
    paste0(
      "Only variables significantly associated with ", focus_var,
      " (p < ", alpha, ")\n",
      "Border: green = assumptions met | red = assumptions violated\n",
      "Stars: * p<.05  ** p<.01  *** p<.001"
    )
  } else {
    paste0(
      "Only variables involved in at least one significant pair (p < ", alpha, ")\n",
      "Border: green = assumptions met | red = assumptions violated\n",
      "Stars: * p<.05  ** p<.01  *** p<.001"
    )
  }
  
  # --- Plot ---
  ggplot(plot_df, aes(x = Var2, y = Var1)) +
    geom_tile(aes(fill = rho), color = "white", linewidth = 0.6) +
    geom_tile(aes(color = border_col), fill = NA, linewidth = 1.1) +
    geom_text(aes(label = label), size = 3, fontface = "bold", color = "black") +
    scale_fill_gradient2(
      low      = "#1565C0",
      mid      = "white",
      high     = "#B71C1C",
      midpoint = 0,
      limits   = c(-1, 1),
      name     = "Spearman ρ"
    ) +
    scale_color_identity() +
    scale_x_discrete(position = "top") +
    labs(title = "Pairwise Spearman correlation heatmap",
         subtitle = subtitle, x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x      = element_text(angle = 90, hjust = 0, vjust = 0.5,
                                      face = "bold"),
      axis.text.y      = element_text(face = "bold"),
      panel.grid       = element_blank(),
      plot.title       = element_text(face = "bold", size = 14),
      plot.subtitle    = element_text(size = 9, color = "grey40"),
      legend.position  = "right"
    )
}


###################### USAGE FOR GDF15 WITH BIOMARKERS ######################


df_biomarkers_GDF15=cbind(df[,1:49], df["GDF15 pg/ml"])

df_biomarkers_GDF15 <- df_biomarkers_GDF15 %>%
  filter(rowSums(!is.na(.)) > 4)

# Remove columns with same values because of limit of detection (otherwise create outliers!)
df_biomarkers_GDF15=subset(df_biomarkers_GDF15, select=-c(IL_4, IL_2, IL_1beta, IL_1alpha))

df_biomarkers_GDF15 <- df_biomarkers_GDF15 %>%
  mutate(
    across(
      where(~ is.character(.x) && any(grepl(",", .x, fixed = TRUE))),
      ~ as.numeric(gsub(",", ".", .x))
    )
  )

# NAs introduced are "<LOQ" values so OK
df_biomarkers_GDF15 %>%
  select(where(is.character)) %>%
  summarise(across(everything(), ~ sum(!grepl("^[0-9,\\.]+$", .x) & !is.na(.x))))

unique(df$AA[!grepl("^[0-9,\\.]+$", df$AA)])

# TEST

results_biomarkers <- test_spearman_pairwise(
  df             = df_biomarkers_GDF15,
  alpha          = 0.05,
  check_autocor  = FALSE,
  maha_threshold = 0.99,
  id_col         = "fondacode"
)
# outliers not influential

# Complete significant tab
print(results_biomarkers) 

# Only suspect pairs
results_biomarkers %>% filter(advice != "Spearman valid")

# Significant and valid pairs
results_biomarkers %>% filter(spearman_p < 0.05, spearman_ok)

# Heatmap with all significant pairs
plot_spearman_heatmap(df = df, results = results_biomarkers, alpha = 0.05)

# Only pairs significantly associated with GDF15
heatmap_biomarkers_GDF15 <- plot_spearman_heatmap(
  df        = df,
  results   = results_biomarkers,
  alpha     = 0.05,
  focus_var = "GDF15 pg/ml"
) ; heatmap_biomarkers_GDF15

# Export with ggsave
ggsave(
  filename = "outcome/heatmap_biomarkers_GDF15.png",
  plot     = heatmap_biomarkers_GDF15,
  width    = 8,      # adjust based on number of variables
  height   = 7,
  dpi      = 300,    # high resolution for publication
  bg       = "white" # white background
)

###################### USAGE FOR GDF15 WITH CLINICAL DATA ######################

# Selected variables for the study
#df_clinical_GDF15=df[,c("age","madrs_","ymrs_num", "fast_","bis10","staya","mars_","mathys_", "psqi_", "als_", "ctq29", "ctq31", "ctq33", "ctq35", "ctq37", "ctq39","GDF15 pg/ml")]

# Exploratory variables
df_clinical_GDF15=cbind(df[,56:81], df[,c("fondacode", "age", "bmi", "ymrs_num", "GDF15 pg/ml")])

df_clinical_GDF15 <- df_clinical_GDF15 %>%
  filter(rowSums(!is.na(.)) > 1)
# no missing values for a complete row

# TEST

results_clinical <- test_spearman_pairwise(
  df             = df_clinical_GDF15,
  alpha          = 0.05,
  check_autocor  = FALSE,
  maha_threshold = 0.99,
  id_col         = "fondacode"
)
# outliers not influential

# Complete significant tab
print(results_clinical) 

# Only suspect pairs
results_clinical %>% filter(advice != "Spearman valid")

# Significant and valid pairs
results_clinical %>% filter(spearman_p < 0.05, spearman_ok)

# Heatmap with all significant pairs
plot_spearman_heatmap(df = df, results = results_clinical, alpha = 0.05)

# Only pairs significantly associated with GDF15
heatmap_clinical_GDF15 <- plot_spearman_heatmap(
  df        = df,
  results   = results_clinical,
  alpha     = 0.05,
  focus_var = "GDF15 pg/ml"
) ; heatmap_clinical_GDF15

# Export with ggsave
ggsave(
  filename = "outcome/heatmap_clinical_GDF15.png",
  plot     = heatmap_clinical_GDF15,
  width    = 8,      # adjust based on number of variables
  height   = 7,
  dpi      = 300,    # high resolution for publication
  bg       = "white" # white background
)
