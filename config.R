library(stringr)

# Paths
# Edit these for your environment (or export PGSC_HOME in your shell).
# PGSC_HOME is this repo's root, where figures and Supp_tables are written.
PGSC_HOME <- Sys.getenv("PGSC_HOME", unset = ".")
# data_dir holds inputs not included in this repo (see README). It should contain:
# r2_out/      processed PGS R2 tables produced by the PGSC method repo
# sig_loci/    per-context significant GxC loci counts
# pop_counts/  per-context individual counts per population
# GENIE/       external GENIE results
data_dir  <- "/path/to/PGSC_data/"
# BioMe replication results, also user-supplied and not in this repo.
BioMe_dir <- paste0(data_dir, "BioMe_updated_results/")

# Derived paths, no need to edit.
home_dir        <- data_dir
fig_dir         <- paste0(PGSC_HOME, "/")
rdata_dir       <- paste0(PGSC_HOME, "/Rdata/")
BioMe_input_dir <- paste0(BioMe_dir, "analysis/")

# Phenotype lists
pheno_list <- c(
  "Alanine_aminotrans674206", "Alkaline_phosphatase674206",
  "Apolipoprotein_a674206", "Apolipoprotein_a_0",
  "Apolipoprotein_b674206", "Apolipoprotein_b_0",
  "Arm_fat-free_mass_left674178", "Arm_fat-free_mass_left_0",
  "Aspartate_aminotrans674206", "Basophill_count674178", "Bilirubin674206",
  "Bilirubin_0", "Birth_weight674178", "BMI674178", "Calcium674178",
  "Cholesterol674178", "Creatinine674178", "Creatinine_urine674206",
  "Cystatin_c674206", "DiastolicBP_auto674178", "Eosinophill_count674178",
  "FEV1_FVC_ratio674206", "Gamma_glutamyltransferase674178", "Glucose674178",
  "HbA1c674178", "HbA1c_0", "HDL674178", "HDL_0",
  "Heel_bone_mineral_density_Tscore674206", "Height674178", "Height_0",
  "Hip_circumference674178", "IGF-1674178", "LDL674178", "LDL_0",
  "Leukocyte_count674178", "Lipoprotein_a674206", "Lipoprotein_a_0",
  "Lymphocyte_count674206", "Mean_corpuscular_volume674178",
  "Monocyte_count674206", "Neutrophill_count674206", "Phosphate674206",
  "Platelet_count674178", "Platelet_volume674206", "Potassium_urine674206",
  "Protein674206", "Pulse_rate674178", "RBC674178",
  "Rheumatoid_factor674206", "Right_hand_grip_strength674178",
  "SHBG674178", "SHBG_0", "Sodium_urine674206", "SystolicBP_auto674178",
  "Testosterone674178", "Testosterone_0", "Triglycerides674178",
  "Urate674178", "Urea674178", "Vitamin_D674178",
  "Waist_circumference674178", "Whole_body_fat_mass674178", "WHRadjBMI_Zhu"
)

scale_list <- c(
  "Apolipoprotein_a674206", "Apolipoprotein_a_0",
  "Apolipoprotein_b674206", "Apolipoprotein_b_0",
  "Arm_fat-free_mass_left674178", "Arm_fat-free_mass_left_0",
  "Bilirubin674206", "Bilirubin_0", "HbA1c674178", "HbA1c_0",
  "HDL674178", "HDL_0", "Height674178", "Height_0",
  "LDL674178", "LDL_0", "Lipoprotein_a674206", "Lipoprotein_a_0",
  "SHBG674178", "SHBG_0", "Testosterone674178", "Testosterone_0"
)

pheno_drop <- c(
  "Apolipoprotein_a_0", "Apolipoprotein_b_0", "Arm_fat-free_mass_left_0",
  "Bilirubin_0", "HbA1c_0", "HDL_0", "LDL_0", "Lipoprotein_a_0",
  "Height_0", "SHBG_0", "Testosterone_0"
)

# Contexts
all_contexts <- c("sex", "age", "statins")
context_mapping <- list(
  sex      = list(clean_context = "Female vs. Male"),
  age      = list(clean_context = "Born in 1934-1951 vs. 1952-1971"),
  statins  = list(clean_context = "Never vs. Current statin usage")
)

context_list <- all_contexts
ctx_clean <- setNames(
  sapply(context_list, function(ctx) context_mapping[[ctx]]$clean_context),
  context_list
)
ctx_cols <- setNames(c('#156c4e', '#a796e8', '#004488'), ctx_clean)

# Populations
valid_pops <- c("white_euro", "afr", "asn")
Ancs <- c('European', 'African', 'Asian')
names(Ancs) <- valid_pops

# PGS type lists
type_list <- c("pgs", "ampPGS", "PGSC")
type_list_valid <- c("pgs", "ampPGS", "PGSC_v_pgs", "PGSC_v_ampPGS")

alt_pgs3 <- type_list_valid[1:3]   # pgs, ampPGS, PGSC_v_pgs  (portability)
alt_pgs2 <- type_list_valid[2:3]   # ampPGS, PGSC_v_pgs       (R2% change)

clean_meth3 <- c("pgs","ampPGS", "PGSC")
clean_meth2 <- c("ampPGS", "PGSC")
names(clean_meth3) <- alt_pgs3
names(clean_meth2) <- alt_pgs2

# Colors
base_cols3 <- c('#B2B2BF', '#6699CC', '#994455')
names(base_cols3) <- alt_pgs3
trans_col <- function(col) paste0(col, "99")
detail_col <- c("#DDDDDD","#9999A1","#2b2d42")
names(detail_col) <- c("grey_1","grey_2","grey_3")
# Misc constants
tau <- c(1e-10, 1e-8, 1e-6, 1e-4, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5)

col_drop_main  <- 4L
col_drop_valid <- 7L

pheno_patterns <- c(
  "IGF-1674178" = "IGF-1",
  "HbA1c674178" = "HbA1c"
)

# BioMe phenotype dictionary
pheno_dict <- c(
  alanine_aminotrans        = "Alanine_aminotrans674206",
  alkaline_phosphatase      = "Alkaline_phosphatase674206",
  apolipoprotein_a          = "Apolipoprotein_a674206",
  aspartate_aminotrans      = "Aspartate_aminotrans674206",
  basophill_count           = "Basophill_count674178",
  bilirubin                 = "Bilirubin674206",
  bmi                       = "BMI674178",
  calcium                   = "Calcium674178",
  cholesterol               = "Cholesterol674178",
  creatinine                = "Creatinine674178",
  creatinine_urine          = "Creatinine_urine674206",
  cystatin_c                = "Cystatin_c674206",
  diastolicbp_auto          = "DiastolicBP_auto674178",
  eosinophill_count         = "Eosinophill_count674178",
  fev1_fvc_ratio            = "FEV1_FVC_ratio674206",
  gamma_glutamyltransferase = "Gamma_glutamyltransferase674178",
  hba1c                     = "HbA1c674178",
  hdl                       = "HDL674178",
  height                    = "Height674178",
  ldl                       = "LDL674178",
  leukocyte_count           = "Leukocyte_count674178",
  lipoprotein_a             = "Lipoprotein_a674206",
  lymphocyte_count          = "Lymphocyte_count674206",
  mean_corpuscular_volume   = "Mean_corpuscular_volume674178",
  monocyte_count            = "Monocyte_count674206",
  neutrophill_count         = "Neutrophill_count674206",
  phosphate                 = "Phosphate674206",
  protein                   = "Protein674206",
  shbg                      = "SHBG674178",
  sodium_urine              = "Sodium_urine674206",
  systolicbp_auto           = "SystolicBP_auto674178",
  testosterone              = "Testosterone674178",
  triglycerides             = "Triglycerides674178",
  urate                     = "Urate674178",
  urea                      = "Urea674178",
  vitamin_d                 = "Vitamin_D674178"
)

# GENIE phenotype dictionary
match_list <- c("Alanine_aminotrans674206" = "Alanine aminotransferase","Alkaline_phosphatase674206" = "Alkaline phosphatase",
                "Apolipoprotein_a674206" = "Apolipoprotein A", "Apolipoprotein_b674206" = NA, "Arm_fat-free_mass_left674178" = NA,
                "Aspartate_aminotrans674206" = "Aspartate aminotransferase", "Basophill_count674178" = NA, "Bilirubin674206" = NA,
                "Birth_weight674178" = NA, "BMI674178" = "Body mass index", "Calcium674178" = NA, "Cholesterol674178" = "Cholesterol",
                "Creatinine674178" = "Creatinine", "Creatinine_urine674206" = "Creatinine in urine", "Cystatin_c674206" = "Cystatin-C",
                "DiastolicBP_auto674178" = "Diastolic blood pressure", "Eosinophill_count674178" = "Eosinophil count",
                "FEV1_FVC_ratio674206" = "FEV1-FVC ratio", "Gamma_glutamyltransferase674178" = "Gamma glutamyltransferase",
                "Glucose674178" = "Glucose", "HbA1c674178" = "Hemoglobin A1c", "HDL674178" = "HDL cholesterol",
                "Heel_bone_mineral_density_Tscore674206" = "BMD Heel T-score", "Height674178" = "Height", "Hip_circumference674178" = NA,
                "IGF-1674178" = "IGF-1", "LDL674178" = "LDL direct", "Leukocyte_count674178" = "White blood cell count",
                "Lipoprotein_a674206" = NA, "Lymphocyte_count674206" = "Lymphocyte count", "Mean_corpuscular_volume674178" = NA,
                "Monocyte_count674206" = "Monocyte count", "Neutrophill_count674206" = NA, "Phosphate674206" = "Phosphate",
                "Platelet_count674178" = "Platelet count", "Platelet_volume674206" = "Mean platelet volume",
                "Potassium_urine674206" = "Potassium in urine", "Protein674206" = NA, "Pulse_rate674178" = NA, "RBC674178" = "RBC count",
                "Rheumatoid_factor674206" = NA, "Right_hand_grip_strength674178" = NA, "SHBG674178" = "SHBG",
                "Sodium_urine674206" = "Sodium in urine", "SystolicBP_auto674178" = "Systolic blood pressure",
                "Testosterone674178" = "Testosterone", "Triglycerides674178" = "Triglycerides", "Urate674178" = "Urate",
                "Urea674178" = "Urea", "Vitamin_D674178" = "Vitamin D", "Waist_circumference674178" = NA,
                "Whole_body_fat_mass674178" = NA, "WHRadjBMI_Zhu" = "Waist-hip ratio")

# Functions
pheno_cleaner <- function(p) {
  p <- gsub("674178", "", p)
  p <- gsub("674206", "", p)
  p <- gsub("_", " ", p)
  p <- str_to_sentence(p)
  p <- gsub("Apolipoprotein a", "Apolipoprotein A", p)
  p <- gsub("Apolipoprotein b", "Apolipoprotein B", p)
  p <- gsub("Bmi", "BMI", p)
  p <- gsub("Cystatin c", "Cystatin C", p)
  p <- gsub("Diastolicbp auto", "Diastolic BP", p)
  p <- gsub("Fev1 fvc ratio", "FEV1 FVC ratio", p)
  p <- gsub("Gamma glutamyltransferase", "Gamma glutamyltrans", p)
  p <- gsub("Hba1c", "HbA1c", p)
  p <- gsub("Hdl", "HDL", p)
  p <- gsub("Heel bone mineral density tscore", "BMD Heel T-score", p)
  p <- gsub("Igf-1", "IGF-1", p)
  p <- gsub("Ldl", "LDL", p)
  p <- gsub("Rbc", "RBC", p)
  p <- gsub("Shbg", "SHBG", p)
  p <- gsub("Systolicbp auto", "Systolic BP", p)
  p <- gsub("Vitamin d", "Vitamin D", p)
  p <- gsub("Whradjbmi zhu", "WHR adj BMI", p)
  p
}

read_file_safely <- function(file_path) {
  if (file.exists(file_path)) {
    return(read.table(file_path, header = TRUE))
  } else {
    missing_info <- sub(".*processed_(.*)_PGS.*", "\\1", file_path)
    missing_files <<- append(missing_files, missing_info)
    return(NULL)
  }
}

scale_by_pgs <- function(mat, baseline) (mat / baseline) * 100

avg_ivw_fn <- function(SD_mat) {
  if (is.null(dim(SD_mat))) SD_mat <- matrix(SD_mat, ncol = 1)
  1 / sqrt(colSums(1 / SD_mat^2, na.rm = TRUE))
}

avg_se_fn <- function(SD_mat) {
  if (is.null(dim(SD_mat))) SD_mat <- matrix(SD_mat, ncol = 1)
  sqrt(apply(SD_mat^2, 2, mean, na.rm = TRUE)) / sqrt(colSums(!is.na(SD_mat)))
}

build_plot_matrices <- function(per_diff_r2, SD, CI25, CI75,
                                sort_col = "PGSC_v_pgs",
                                avg_r2s = NULL, avg_se = NULL) {
  ord         <- order(per_diff_r2[, sort_col], na.last = NA)
  per_diff_r2 <- per_diff_r2[ord, , drop = FALSE]
  SD          <- SD[ord,   , drop = FALSE]
  CI25        <- CI25[ord, , drop = FALSE]
  CI75        <- CI75[ord, , drop = FALSE]
  avg_r2s <- colMeans(per_diff_r2, na.rm = TRUE)
  avg_se  <- avg_se_fn(SD)
  list(
    data = rbind(avg_r2s,               rep(0, length(avg_r2s)), per_diff_r2),
    CI25 = rbind(avg_r2s - 1.96*avg_se, rep(0, length(avg_r2s)), CI25),
    CI75 = rbind(avg_r2s + 1.96*avg_se, rep(0, length(avg_r2s)), CI75)
  )
}

fmt_p <- function(p) {
  if (is.na(p))  return("")
  if (p < 0.001) return("p<0.001")
  if (p < 0.01)  return(sprintf("p=%.3f", p))
  sprintf("p=%.2f", p)
}

cor_label <- function(name, cor_res) {
  prefix <- if (nchar(name) > 0) paste0(name, " ") else ""
  paste0(prefix, "Cor: ", round(cor_res$estimate, 2), ", pv: ", round(cor_res$p.value, 2))
}

wald_z <- function(d, sd_a, sd_b)
  2 * pnorm(-abs(mean(d, na.rm=TRUE) / (sqrt(mean(sd_a^2 + sd_b^2, na.rm=TRUE)) / sqrt(sum(!is.na(d))))))

wald_z_ivw <- function(d, sd_a, sd_b) {
  var_i <- sd_a^2 + sd_b^2
  w <- 1 / var_i
  w[is.na(w) | is.na(d)] <- NA
  se_ivw <- 1 / sqrt(sum(w, na.rm = TRUE))
  2 * pnorm(-abs(mean(d, na.rm = TRUE) / se_ivw))
}

z_test_1s <- function(d)
  2 * pnorm(-abs(mean(d, na.rm=TRUE) / (sd(d, na.rm=TRUE) / sqrt(sum(!is.na(d))))))

r2diffplot <- function(data_matrix, CI25_matrix, CI75_matrix, filename,
                       meth, clean_meth, clean_phenos, clean_context, pop, context) {
  if(pop == "BioMe"){
    clean_pop <- "BioMe"
    h_val <- 12
  }else{
    clean_pop <- Ancs[pop]
    h_val <- 15
  }
  png(filename, width = 10, height = h_val, units = 'in', res = 300)
  par(mai = c(1, 2.95, 0.25, 0.25))

  plot(x = data_matrix[, meth[1]], y = seq_along(data_matrix[, meth[1]]),
       type = "p", pch = 16,
       col  = ifelse(data_matrix[, meth[1]] == 0, "white", base_cols[meth[1]]),
       xlab = paste0("R\u00b2 %Change in ", clean_pop, " population"), ylab = "",
       yaxt = "n", las = 1, cex.lab = 2, cex = 1.5, cex.axis = 1.5,
       xlim = c(min(data_matrix[, meth]) - 2, max(data_matrix[, meth]) + 2))

  axis_labels <- c("Avg across phenos", " ", clean_phenos[-c(1:2)])
  axis(2, at = seq_along(data_matrix[, meth[1]]), labels = axis_labels, las = 1, cex.axis = 1.3)

  for (m in meth) {
    points(data_matrix[, m], y = seq_along(data_matrix[, m]),
           type = "p", pch = 16, cex = 1.5,
           col  = ifelse(data_matrix[, m] == 0, "white", base_cols[m]))
    arrows(x0 = CI25_matrix[, m], y0 = seq_along(data_matrix[, m]),
           x1 = CI75_matrix[, m], y1 = seq_along(data_matrix[, m]),
           angle = 90, code = 3, lwd = 2, length = 0.05, col = base_cols[m])
  }

  abline(v = 0, col = "black")
  abline(h = 2, col = "black")
  legend('bottomright', bty = 'n', clean_context, cex = 1.5)
  if (context == "sex") {
    legend("right", legend = clean_meth, bty = 'n',
           col = base_cols[meth], pch = 16, cex = 1.5, lwd = 2)
  }
  dev.off()
}
