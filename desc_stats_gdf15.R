
install.packages(c("dplyr", "gtsummary", "stringr", "openxlsx", "readxl"))

setwd("~/IMRB/Stats GDF15")

library(dplyr)
library(gtsummary)
library(stringr)
library(openxlsx)
library(readxl)

######################### DESCRIPTIVE TABLE ############################

# 1) DATA IMPORT

# Import data
data=read_excel("data/Complete_data_Evanne_GDF15_BD_cleaned.xlsx", sheet = "V0") ; head(data)
data["mdq_"]

# Clean empty columns
data_noblanks=data[, apply(data, 2, function(x) !all(is.na(x)))]

# Import available samples
data_samples=read_excel("data/samples_V0V1V2_Evanne_GDF15.xlsx") ; head(data_samples)
data_samples_V0=data_samples[data_samples$VISITE=='V0',]

# Merge data with available samples we have by patient ID
merged_data_V0=merge(data_noblanks, data_samples_V0[, c("fondacode", "GDF15 pg/ml")], by = "fondacode")
# sample id 0201805210 missing in database

table(merged_data_V0$arm, useNA = "ifany")

# Change age class from character to numeric
sapply(merged_data_V0["agedebutpremier_episode"], class)
merged_data_V0=transform(merged_data_V0, agedebutpremier_episode = as.numeric(agedebutpremier_episode))


# 2) VALUE RECODING

merged_data_V0 <- merged_data_V0 %>%
  
  mutate(
    # Gender: sex = feminin/masculin -> Female/Male
    
    sex = case_when(
      
      str_to_lower(str_trim(as.character(sex))) == "feminin" ~ "Female",
      
      str_to_lower(str_trim(as.character(sex))) == "masculin" ~ "Male",
      
      TRUE ~ NA_character_) %>% factor(levels = c("Female", "Male")),
    
    
    # BD diagnosis: arm -> BD type I / BD type II / BD not specified
    
    arm = case_when(
      
      str_detect(str_to_lower(str_trim(as.character(arm))), "type\\s*1|type\\s*i") ~ "BD type I",
      
      str_detect(str_to_lower(str_trim(as.character(arm))), "type\\s*2|type\\s*ii") ~ "BD type II",
      
      is.na(arm) | str_trim(as.character(arm)) == "Bipolaire non spécifié" ~ "BD not specified",
      
      TRUE ~ "BD not specified") %>% factor(levels = c("BD type I", "BD type II", "BD not specified")),
    
    
    # Smoking status
    
    suncf_cigarettes_lt = case_when(

      str_trim(as.character(suncf_cigarettes_lt)) == "Non fumeur" ~ "No",

      str_trim(as.character(suncf_cigarettes_lt)) == "Fumeur actuel" ~ "Yes",
      
      str_trim(as.character(suncf_cigarettes_lt)) == "Ex-fumeur" ~ "Remitted",

      TRUE ~ NA_character_) %>% factor(levels = c("No", "Yes", "Remitted")),


    # Substance disorder
    
    rad_tb_subst = case_when(

      str_trim(as.character(rad_tb_subst)) == "N" ~ "No",

      str_trim(as.character(rad_tb_subst)) == "Y" ~ "Yes",
      
      str_trim(as.character(rad_tb_subst)) == "U" ~ NA,

      TRUE ~ NA_character_) %>% factor(levels = c("No", "Yes")),


    # Alcohol status
    
    suoccur_alcool = ifelse(suoccur_alcool == 1, "Yes", "No"),
    
    # Cannabis status

    suoccur_cannabis = ifelse(suoccur_cannabis == 1, "Yes", "No"),
    
    # Treatment status
    
    Antidepressants_treat = ifelse(Antidepressant == 1, "Yes", "No"),
    
    Anxiolytics_treat = ifelse(anxio_hypno == 1, "Yes", "No"),
    
    Lithium_treat = ifelse(Lithium == 1, "Yes", "No"),
    
    Antipsychotics_treat = ifelse(FGA == 1 | APA == 1 | Antipsychotic == 1, "Yes", "No"),
    
    Thymoregulators_treat = ifelse(Valproate == 1 | thymoACAE == 1 | Thymoregulator == 1, "Yes", "No")
    )

# Clean the table and export

merged_data_V0 = merged_data_V0[-c(4:12)]
write.xlsx(merged_data_V0, "data/merged_data_samples_V0.xlsx")

# Verification for some variables 

table(merged_data_V0$sex, useNA = "ifany")

table(merged_data_V0$arm, useNA = "ifany")

table(merged_data_V0$rad_tb_subst, useNA = "ifany")


# 3) ARTICLE TABLE

# Socio demographic characteristics

tab_socio <- merged_data_V0 %>%
  
  select(sex, age, arm, agedebutpremier_episode, madrs_, ymrs_num) %>%
  
  tbl_summary(
    
    statistic = list(
      all_continuous() ~ "{mean} ± {sd}",
      all_categorical() ~ "{n} ({p}%)"),
    
    missing = "no",
    
    label = list(
      sex           ~ "Gender, n (%)",
      age           ~ "Age, mean ± SD",
      arm           ~ "BD diagnosis, n (%)",
      agedebutpremier_episode     ~ "Age at onset of BD, mean ± SD",
      madrs_         ~ "MADRS, mean ± SD",
      ymrs_num         ~ "YMRS, mean ± SD")
    ) ; tab_socio

# Comorbidities

tab_comorb <- merged_data_V0 %>%
  
  select(bmi, suncf_cigarettes_lt, rad_tb_subst, suoccur_alcool, suoccur_cannabis) %>%
  
  tbl_summary(
    
    statistic = list(
      all_continuous() ~ "{mean} ± {sd}",
      all_categorical() ~ "{n} ({p}%)"),
    
    missing = "ifany",
    missing_text = "Missing",
    
    label = list(
      bmi  ~ "BMI, mean ± SD",
      suncf_cigarettes_lt  ~ "Lifetime tobacco status (current smokers), n (%)",
      rad_tb_subst  ~ "Lifetime substance use disorder, n (%)",
      suoccur_alcool ~ "Lifetime alcohol use disorder, n (%)",
      suoccur_cannabis  ~ "Lifetime cannabis use disorder, n (%)"),
  
    # Only display "Yes" values
    value = list(
      suncf_cigarettes_lt ~ "Yes",
      rad_tb_subst     ~ "Yes",
      suoccur_alcool   ~ "Yes",
      suoccur_cannabis ~ "Yes"
    )) %>%
      
      # Remove the "Unknown" row for continuous variables only
      modify_table_body(
        ~.x %>%
          filter(!(row_type == "missing" & var_type == "continuous"))
      ) ; tab_comorb

# Current medication

tab_med <- merged_data_V0 %>%
  
  select(Antidepressants_treat, Anxiolytics_treat, Thymoregulators_treat, Antipsychotics_treat, Lithium_treat) %>%
  
  tbl_summary(
    
    statistic = everything() ~ "{n} ({p}%)",
    
    missing = "ifany",
    missing_text = "Missing",
    
    label = list(
      Antidepressants_treat  ~ "Antidepressants, n (%)",
      Anxiolytics_treat  ~ "Anxiolytics, n (%)",
      Thymoregulators_treat ~ "Mood stabilizers, n (%)",
      Antipsychotics_treat ~ "Antipsychotics, n (%)",
      Lithium_treat  ~ "Lithium, n (%)"),
    
    value = everything() ~ "Yes",
    
  ) ; tab_med

# Merge tables

tab_article <- tbl_stack(
  
  list(tab_socio, tab_comorb, tab_med),
  
  group_header = c("Socio demographic characteristics", "Comorbidities", "Current Medication")) %>%
  
  as_gt() %>%
  gt::tab_style(
    style = gt::cell_text(weight = "bold"),
    locations = gt::cells_row_groups(groups = everything())
  ) ; tab_article

# Export

tab_article %>%
  gt::tab_options(
    table.font.size = gt::px(9),          # police plus petite
    data_row.padding = gt::px(2),         # moins d'espace entre lignes
    column_labels.padding = gt::px(2),
    table.width = gt::pct(50)
  ) %>%
  gt::gtsave("outcome/desc_stats_gdf15.pdf")

