rm(list=ls())
source(file.path(Sys.getenv("PGSC_HOME", unset = "."), "config.R"))

# Part 1: Compile UKB output
output <- array(NA, dim = c(length(pheno_list), length(type_list), 9, length(all_contexts)),
                dimnames = list(pheno_list, type_list,
                                c('tau0','tau1','r2','w0','w0a','w0b','w1a','w1b','w1'),
                                all_contexts))
output_valid <- array(NA, dim = c(length(pheno_list), length(type_list_valid), 10, length(all_contexts), length(valid_pops)),
                      dimnames = list(pheno_list, type_list_valid,
                                      c('thresh','r2','w0','sd','se','pval','delta','w1','CI25','CI75'),
                                      all_contexts, valid_pops))
output_scaled <- array(NA, dim = c(length(pheno_list), length(type_list), 9, length(all_contexts)),
                       dimnames = list(pheno_list, type_list,
                                       c('tau0','tau1','r2','w0','w0a','w0b','w1a','w1b','w1'),
                                       all_contexts))
missing_files <- list()

for (context in all_contexts) {
  for (pheno in pheno_list) {
    matched       <- names(pheno_patterns)[sapply(names(pheno_patterns), grepl, pheno)]
    phenoNoDigits <- if (length(matched)) pheno_patterns[matched[1]] else if (grepl("_(0|2)$", pheno)) pheno else gsub("[0-9]+$", "", pheno)
    phenoLower   <- tolower(phenoNoDigits)

    pgsC                    <- read_file_safely(paste0(home_dir, "r2_out/processed_",         phenoLower, "_whitebrit_",  context, "_PGS.txt"))
    pgsC_scaled_whitebrit   <- read_file_safely(paste0(home_dir, "r2_out/scaled/processed_",  phenoLower, "_whitebrit_",  context, "_PGS_scaled.txt"))
    pgsC_valid_white_euro   <- read_file_safely(paste0(home_dir, "r2_out/valid/processed_",   phenoLower, "_white_euro_", context, "_PGS.txt"))
    pgsC_valid_afr          <- read_file_safely(paste0(home_dir, "r2_out/valid/processed_",   phenoLower, "_afr_",        context, "_PGS.txt"))
    pgsC_valid_asn          <- read_file_safely(paste0(home_dir, "r2_out/valid/processed_",   phenoLower, "_asn_",        context, "_PGS.txt"))

    if (is.null(pgsC)) next
    tau0 <- pgsC[pgsC[, "type"] == 'pgs', 'Threshold']
    if (length(tau0) != 1) tau0 <- tau0[1]

    for (pgs_type in type_list) {
      update_pgsC <- as.matrix(pgsC[pgsC[, "type"] == pgs_type, -col_drop_main])
      if (pgs_type %in% c("pgs", "ampPGS")) {
        output[pheno, pgs_type, -2, context] <- update_pgsC
      } else {
        output[pheno, pgs_type, , context] <- cbind(tau0, update_pgsC)
      }
    }

    if (!is.null(pgsC_scaled_whitebrit)) {
      tau0_scaled <- pgsC_scaled_whitebrit[pgsC_scaled_whitebrit[, "type"] == "pgs", "Threshold"]
      if (length(tau0_scaled) != 1) tau0_scaled <- tau0_scaled[1]
      for (pgs_type in type_list) {
        update_pgsC <- as.matrix(pgsC_scaled_whitebrit[pgsC_scaled_whitebrit[, "type"] == pgs_type, -col_drop_main])
        if (pgs_type %in% c("pgs", "ampPGS")) {
          output_scaled[pheno, pgs_type, -2, context] <- update_pgsC
        } else {
          output_scaled[pheno, pgs_type, , context] <- cbind(tau0_scaled, update_pgsC)
        }
      }
    }

    pgsC_valid_list <- setNames(list(pgsC_valid_white_euro, pgsC_valid_afr, pgsC_valid_asn), valid_pops)
    for (pgs_type_v in type_list_valid) {
      for (pop in names(pgsC_valid_list)) {
        pgs_data <- pgsC_valid_list[[pop]]
        if (is.null(pgs_data)) next
        valid_rows <- pgs_data[, "type"] == pgs_type_v
        if (any(valid_rows))
          output_valid[pheno, pgs_type_v, , context, pop] <- as.matrix(pgs_data[valid_rows, -col_drop_valid])
      }
    }
  }
}

if (length(missing_files) > 0) {
  missing_vec <- unlist(missing_files)
  cat(paste0("The following ", length(missing_files), " files were missing and skipped:\n"))
  print(head(missing_files))
  print(paste0("WB: ",  length(grep("whitebrit",  missing_vec, value=TRUE))))
  print(paste0("eur: ", length(grep("white_euro", missing_vec, value=TRUE))))
  print(paste0("afr: ", length(grep("afr",        missing_vec, value=TRUE))))
  print(paste0("asn: ", length(grep("asn",        missing_vec, value=TRUE))))
}

pheno_sig_drop <- as.character(na.omit(
  pheno_list[output_valid[, 'pgs', 'r2', context_list[3], 'white_euro'] < 0.01]
))
plot_phenos <- setdiff(pheno_list, c(pheno_drop, pheno_sig_drop))

save(output, output_valid, output_scaled,
     valid_pops, Ancs, pheno_list, pheno_drop, pheno_sig_drop, plot_phenos,
     all_contexts, context_list, context_mapping, type_list, type_list_valid,
     scale_list, tau, trans_col, pheno_cleaner,
     file = paste0(rdata_dir, "compiled_output.Rdata"))


# Part 2: Compute shared stats

# UKB per-pheno arrays
pgs_r2_arr <- array(NA, dimnames = list(plot_phenos, context_list, valid_pops),
                    dim = c(length(plot_phenos), length(context_list), length(valid_pops)))
per_diff_r2_arr <- array(NA, dimnames = list(plot_phenos, alt_pgs2, context_list, valid_pops),
                         dim = c(length(plot_phenos), length(alt_pgs2), length(context_list), length(valid_pops)))
SD_pct_arr   <- array(NA, dimnames = list(plot_phenos, alt_pgs2, context_list, valid_pops),
                      dim = c(length(plot_phenos), length(alt_pgs2), length(context_list), length(valid_pops)))
CI25_pct_arr <- array(NA, dimnames = list(plot_phenos, alt_pgs2, context_list, valid_pops),
                      dim = c(length(plot_phenos), length(alt_pgs2), length(context_list), length(valid_pops)))
CI75_pct_arr <- array(NA, dimnames = list(plot_phenos, alt_pgs2, context_list, valid_pops),
                      dim = c(length(plot_phenos), length(alt_pgs2), length(context_list), length(valid_pops)))

for (pop in valid_pops) {
  for (context in context_list) {
    base_r2 <- output_valid[plot_phenos, 'pgs', 'r2', context, pop]
    pgs_r2_arr[, context, pop]        <- base_r2
    per_diff_r2_arr[, , context, pop] <- scale_by_pgs(output_valid[plot_phenos, alt_pgs2, 'delta', context, pop], base_r2)
    SD_pct_arr[, , context, pop]      <- scale_by_pgs(output_valid[plot_phenos, alt_pgs2, 'sd',    context, pop], base_r2)
    CI25_pct_arr[, , context, pop]    <- scale_by_pgs(output_valid[plot_phenos, alt_pgs2, 'CI25',  context, pop], base_r2)
    CI75_pct_arr[, , context, pop]    <- scale_by_pgs(output_valid[plot_phenos, alt_pgs2, 'CI75',  context, pop], base_r2)
  }
}

# context_summary & pvals_vs0
context_summary <- array(NA, dimnames = list(context_list, valid_pops, alt_pgs2, c('r2','se')),
                         dim = c(length(context_list), length(valid_pops), length(alt_pgs2), 2))
pvals_vs0 <- array(NA, dimnames = list(context_list, valid_pops, alt_pgs2[1:2]),
                   dim = c(length(context_list), length(valid_pops), 2))

for (pop in valid_pops) {
  for (context in context_list) {
    keep    <- which(pgs_r2_arr[, context, pop] > 0.01)
    pdr     <- per_diff_r2_arr[keep, , context, pop]
    sdr     <- SD_pct_arr[keep, , context, pop]
    avg_r2s <- colMeans(pdr, na.rm=TRUE)
    avg_se  <- avg_se_fn(sdr)
    context_summary[context, pop, , ] <- cbind(avg_r2s, avg_se)
    pvals_vs0[context, pop, ] <- 2 * pnorm(-abs(avg_r2s[alt_pgs2[1:2]] / avg_se[alt_pgs2[1:2]]))
  }
}

# avg_R2diff & pval_R2diff
avg_R2diff  <- array(NA, dimnames = list(context_list, valid_pops[2:3], c("ampPGS","PGSC"), c('delta','se')),
                     dim = c(3, 2, 2, 2))
pval_R2diff <- array(NA, dimnames = list(context_list, valid_pops[2:3]), dim = c(3, 2))

for (anc in valid_pops[2:3]) {
  for (context in context_list) {
    keep          <- which(pgs_r2_arr[, context, anc] > 0.01)
    anc_amp_diff  <- per_diff_r2_arr[keep, 'ampPGS',     context, anc]
    anc_pgsc_diff <- per_diff_r2_arr[keep, 'PGSC_v_pgs', context, anc]
    anc_amp_sd    <- SD_pct_arr[keep,      'ampPGS',     context, anc]
    anc_pgsc_sd   <- SD_pct_arr[keep,      'PGSC_v_pgs', context, anc]
    avg_R2diff[context, anc, 'ampPGS', 'delta'] <- mean(anc_amp_diff,  na.rm=TRUE)
    avg_R2diff[context, anc, 'PGSC',   'delta'] <- mean(anc_pgsc_diff, na.rm=TRUE)
    avg_R2diff[context, anc, 'ampPGS', 'se']    <- avg_se_fn(anc_amp_sd)
    avg_R2diff[context, anc, 'PGSC',   'se']    <- avg_se_fn(anc_pgsc_sd)
    pval_R2diff[context, anc] <- wald_z(anc_pgsc_diff - anc_amp_diff, anc_pgsc_sd, anc_amp_sd)
  }
}

# avg_port & pval_port
avg_port  <- array(NA, dimnames = list(context_list, valid_pops[2:3], c("pgs","ampPGS","PGSC"), c("mean","se")),
                   dim = c(3, 2, 3, 2))
pval_port <- array(NA, dimnames = list(context_list, valid_pops[2:3], c("pgs_amp","pgs_pgsc","amp_pgsc")),
                   dim = c(3, 2, 3))

for (anc in valid_pops[2:3]) {
  for (context in context_list) {
    pheno_loc     <- plot_phenos[which(output_valid[plot_phenos, 'pgs', 'r2', context, anc] > 0.01)]
    anc_pgs_port  <- output_valid[pheno_loc, 'pgs',        'r2', context, anc] / output_valid[pheno_loc, 'pgs',        'r2', context, valid_pops[1]]
    anc_amp_port  <- output_valid[pheno_loc, 'ampPGS',     'r2', context, anc] / output_valid[pheno_loc, 'ampPGS',     'r2', context, valid_pops[1]]
    anc_pgsc_port <- output_valid[pheno_loc, 'PGSC_v_pgs', 'r2', context, anc] / output_valid[pheno_loc, 'PGSC_v_pgs', 'r2', context, valid_pops[1]]
    anc_pgs_sd    <- output_valid[pheno_loc, 'pgs',        'sd', context, anc]
    anc_amp_sd    <- output_valid[pheno_loc, 'ampPGS',     'sd', context, anc]
    anc_pgsc_sd   <- output_valid[pheno_loc, 'PGSC_v_pgs', 'sd', context, anc]
    avg_port[context, anc, 'pgs',    'mean'] <- mean(anc_pgs_port,  na.rm=TRUE)
    avg_port[context, anc, 'ampPGS', 'mean'] <- mean(anc_amp_port,  na.rm=TRUE)
    avg_port[context, anc, 'PGSC',   'mean'] <- mean(anc_pgsc_port, na.rm=TRUE)
    avg_port[context, anc, 'pgs',    'se']   <- avg_se_fn(anc_pgs_sd)
    avg_port[context, anc, 'ampPGS', 'se']   <- avg_se_fn(anc_amp_sd)
    avg_port[context, anc, 'PGSC',   'se']   <- avg_se_fn(anc_pgsc_sd)
    pval_port[context, anc, 'pgs_amp']  <- wald_z(anc_pgs_port  - anc_amp_port,  anc_pgs_sd, anc_amp_sd)
    pval_port[context, anc, 'pgs_pgsc'] <- wald_z(anc_pgs_port  - anc_pgsc_port, anc_pgs_sd, anc_pgsc_sd)
    pval_port[context, anc, 'amp_pgsc'] <- wald_z(anc_amp_port  - anc_pgsc_port, anc_amp_sd, anc_pgsc_sd)
  }
}

# port_gain_arr & labeled_df
all_phenos <- unique(unlist(lapply(valid_pops[2:3], function(anc)
  lapply(context_list, function(ctx)
    plot_phenos[which(output_valid[plot_phenos, 'pgs', 'r2', ctx, anc] > 0.01)]))))

port_gain_arr <- array(NA, dimnames = list(all_phenos, context_list, valid_pops[2:3]),
                       dim = c(length(all_phenos), length(context_list), 2))
for (anc in valid_pops[2:3]) {
  for (ctx in context_list) {
    pheno_loc <- intersect(all_phenos, plot_phenos[which(output_valid[plot_phenos, 'pgs', 'r2', ctx, anc] > 0.01)])
    pgs_port  <- output_valid[pheno_loc, 'pgs',        'r2', ctx, anc] / output_valid[pheno_loc, 'pgs',        'r2', ctx, valid_pops[1]]
    pgsc_port <- output_valid[pheno_loc, 'PGSC_v_pgs', 'r2', ctx, anc] / output_valid[pheno_loc, 'PGSC_v_pgs', 'r2', ctx, valid_pops[1]]
    port_gain_arr[pheno_loc, ctx, anc] <- pgsc_port - pgs_port
  }
}

labeled_rows <- list()
for (anc in valid_pops[2:3]) {
  for (context in context_list) {
    pheno_loc   <- plot_phenos[which(output_valid[plot_phenos, 'pgs', 'r2', context, anc] > 0.01)]
    we_d        <- scale_by_pgs(output_valid[pheno_loc, 'PGSC_v_pgs', 'delta', context, valid_pops[1]],
                                output_valid[pheno_loc, 'pgs',        'r2',    context, valid_pops[1]])
    pop_d       <- scale_by_pgs(output_valid[pheno_loc, 'PGSC_v_pgs', 'delta', context, anc],
                                output_valid[pheno_loc, 'pgs',        'r2',    context, anc])
    discord     <- abs(we_d - pop_d)
    top5_thresh <- sort(discord, decreasing=TRUE)[min(5, length(discord))]
    labeled_rows[[length(labeled_rows)+1]] <- data.frame(
      Phenotype  = pheno_loc,
      Population = unname(Ancs[anc]),
      Context    = unname(ctx_clean[context]),
      WE_delta   = we_d,
      Pop_delta  = pop_d,
      Label      = ifelse(discord >= top5_thresh, pheno_cleaner(pheno_loc), NA_character_),
      stringsAsFactors = FALSE
    )
  }
}
labeled_df <- do.call(rbind, labeled_rows)

# Part 3: BioMe statistics
pheno_count <- read.csv(paste0(BioMe_dir, "BioMe_WE_AJ_Phenotype_counts.csv"))

files <- list.files(BioMe_input_dir, pattern = "*analysis.txt", full.names = TRUE)
BioMe_data <- do.call(rbind, lapply(files, function(file) {
  df     <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
  df_sub <- subset(df, type %in% c("pgs", "ampPGS", "PGSCgen"))
  df_sub$trait <- sub("(.*)_analysis.txt", "\\1", basename(file))
  df_sub
}))
BioMe_data <- BioMe_data[BioMe_data$trait %in% names(pheno_dict), ]

biome_context    <- "sex"
biome_type_names <- setNames(c("pgs", "ampPGS", "PGSCgen"), alt_pgs3)

base_r2_biome <- BioMe_data[BioMe_data$type == "pgs", "R2"]
biome_data_summary <- array(NA,
                            dimnames = list(names(pheno_dict), alt_pgs3,
                                            c('R2','sd_raw','per_R2diff','sd_pct','CI25_pct','CI75_pct')),
                            dim = c(length(pheno_dict), length(alt_pgs3), 6))
biome_data_summary[, "pgs", "R2"]     <- base_r2_biome
biome_data_summary[, "pgs", "sd_raw"] <- BioMe_data[BioMe_data$type == "pgs", "sd"]

for (pgs_type in alt_pgs3[2:3]) {
  temp_data <- BioMe_data[BioMe_data$type == biome_type_names[pgs_type], ]
  biome_data_summary[, pgs_type, "R2"]         <- temp_data$R2
  biome_data_summary[, pgs_type, "sd_raw"]     <- temp_data$sd
  biome_data_summary[, pgs_type, "per_R2diff"] <- scale_by_pgs(temp_data$r2_delta, base_r2_biome)
  biome_data_summary[, pgs_type, "sd_pct"]     <- scale_by_pgs(temp_data$sd,       base_r2_biome)
  biome_data_summary[, pgs_type, "CI25_pct"]   <- scale_by_pgs(temp_data$CI_25,    base_r2_biome)
  biome_data_summary[, pgs_type, "CI75_pct"]   <- scale_by_pgs(temp_data$CI_75,    base_r2_biome)
}

biome_pgs_r2_all   <- biome_data_summary[, "pgs",        "R2"]
biome_per_diff_all <- biome_data_summary[, alt_pgs3[2:3], "per_R2diff"]
biome_sd_pct_all   <- biome_data_summary[, alt_pgs3[2:3], "sd_pct"]
biome_CI25_all     <- biome_data_summary[, alt_pgs3[2:3], "CI25_pct"]
biome_CI75_all     <- biome_data_summary[, alt_pgs3[2:3], "CI75_pct"]

biome_keep      <- which(biome_data_summary[, "pgs", "R2"] > 0.01)
pheno_keep_list <- na.exclude(pheno_count[as.numeric(pheno_count$no.median.obs) > 500, ])$UKB.pheno.match
crop_fn <- function(mat) mat[rownames(mat) %in% pheno_keep_list, , drop=FALSE]

biome_per_diff_crop <- crop_fn(biome_data_summary[biome_keep, alt_pgs3[2:3], "per_R2diff"])
biome_sd_pct_crop   <- crop_fn(biome_data_summary[biome_keep, alt_pgs3[2:3], "sd_pct"])
biome_CI25_crop     <- crop_fn(biome_data_summary[biome_keep, alt_pgs3[2:3], "CI25_pct"])
biome_CI75_crop     <- crop_fn(biome_data_summary[biome_keep, alt_pgs3[2:3], "CI75_pct"])
biome_r2_crop       <- crop_fn(biome_data_summary[biome_keep, alt_pgs3,      "R2"])
biome_sd_raw_crop   <- crop_fn(biome_data_summary[biome_keep, alt_pgs3,      "sd_raw"])

biome_traits     <- rownames(biome_per_diff_crop)
biome_traits_ukb <- pheno_dict[biome_traits]
keep_match       <- biome_traits_ukb %in% dimnames(output_valid)[[1]]
biome_traits     <- biome_traits[keep_match]
biome_traits_ukb <- biome_traits_ukb[keep_match]

amp_diff    <- biome_per_diff_crop[biome_traits, alt_pgs3[2]]
pgsc_diff   <- biome_per_diff_crop[biome_traits, alt_pgs3[3]]
amp_sd_pct  <- biome_sd_pct_crop[biome_traits,   alt_pgs3[2]]
pgsc_sd_pct <- biome_sd_pct_crop[biome_traits,   alt_pgs3[3]]

biome_r2diff_means <- c(mean(amp_diff,  na.rm=TRUE), mean(pgsc_diff, na.rm=TRUE))
biome_r2diff_ses   <- c(avg_se_fn(amp_sd_pct), avg_se_fn(pgsc_sd_pct))
names(biome_r2diff_means) <- names(biome_r2diff_ses) <- alt_pgs3[2:3]
biome_r2diff_pval <- wald_z(pgsc_diff - amp_diff, pgsc_sd_pct, amp_sd_pct)

eur_pgs_r2  <- output_valid[biome_traits_ukb, 'pgs',        'r2', biome_context, valid_pops[1]]
eur_amp_r2  <- output_valid[biome_traits_ukb, 'ampPGS',     'r2', biome_context, valid_pops[1]]
eur_pgsc_r2 <- output_valid[biome_traits_ukb, 'PGSC_v_pgs', 'r2', biome_context, valid_pops[1]]

biome_port_pgs  <- biome_r2_crop[biome_traits, alt_pgs3[1]] / eur_pgs_r2
biome_port_amp  <- biome_r2_crop[biome_traits, alt_pgs3[2]] / eur_amp_r2
biome_port_pgsc <- biome_r2_crop[biome_traits, alt_pgs3[3]] / eur_pgsc_r2

pgs_sd_raw  <- biome_sd_raw_crop[biome_traits, alt_pgs3[1]]
amp_sd_raw  <- biome_sd_raw_crop[biome_traits, alt_pgs3[2]]
pgsc_sd_raw <- biome_sd_raw_crop[biome_traits, alt_pgs3[3]]

biome_port_means <- c(mean(biome_port_pgs,  na.rm=TRUE),
                      mean(biome_port_amp,  na.rm=TRUE),
                      mean(biome_port_pgsc, na.rm=TRUE))
biome_port_ses   <- c(avg_se_fn(pgs_sd_raw), avg_se_fn(amp_sd_raw), avg_se_fn(pgsc_sd_raw))
names(biome_port_means) <- names(biome_port_ses) <- alt_pgs3

biome_port_pvals <- c(
  pgs_amp  = wald_z(biome_port_pgs - biome_port_amp,  pgs_sd_raw, amp_sd_raw),
  pgs_pgsc = wald_z(biome_port_pgs - biome_port_pgsc, pgs_sd_raw, pgsc_sd_raw),
  amp_pgsc = wald_z(biome_port_amp - biome_port_pgsc, amp_sd_raw, pgsc_sd_raw)
)

# Save summ_stats.Rdata
save(
  # shared meta
  plot_phenos, alt_pgs3, alt_pgs2, context_list, valid_pops, Ancs,
  ctx_clean, ctx_cols, base_cols3, clean_meth2,
  # helper functions
  scale_by_pgs, fmt_p, wald_z, z_test_1s,
  # UKB per-pheno arrays
  pgs_r2_arr, per_diff_r2_arr, SD_pct_arr, CI25_pct_arr, CI75_pct_arr,
  # context summaries
  context_summary, pvals_vs0,
  # avg R2% change
  avg_R2diff, pval_R2diff,
  # avg portability
  avg_port, pval_port,
  # per-pheno portability gain & scatter data
  port_gain_arr, labeled_df,
  # BioMe stats
  pheno_dict, biome_traits, biome_traits_ukb,
  biome_r2diff_means, biome_r2diff_ses, biome_r2diff_pval,
  biome_port_means, biome_port_ses, biome_port_pvals,
  # BioMe per-pheno arrays (all phenos)
  biome_pgs_r2_all, biome_per_diff_all, biome_sd_pct_all, biome_CI25_all, biome_CI75_all,
  # BioMe per-pheno arrays (cropped: R2>0.01 & pop count>500)
  biome_per_diff_crop, biome_sd_pct_crop, biome_CI25_crop, biome_CI75_crop,
  file = paste0(rdata_dir, "summ_stats.Rdata")
)

# Part 4: Supplementary CSV tables
setwd(paste0(fig_dir, "Supp_tables/summ_stats/"))
alt_pgs_tbl <- alt_pgs3

# Read per-phenotype misc stats (white_euro only): c_effect and pgs_gxc_cor
# c_effect    = standardized main effect of context C on scaled phenotype
# pgs_gxc_cor = Pearson correlation between scaled PGS and scaled PGxCS
misc_dir        <- paste0(home_dir, "r2_out/scaled/valid/")
c_effect_mat    <- matrix(NA_real_, nrow = length(pheno_list), ncol = length(context_list),
                           dimnames = list(pheno_list, context_list))
pgs_gxc_cor_mat <- matrix(NA_real_, nrow = length(pheno_list), ncol = length(context_list),
                            dimnames = list(pheno_list, context_list))
pgs_pv_mat      <- matrix(NA_real_, nrow = length(pheno_list), ncol = length(context_list),
                           dimnames = list(pheno_list, context_list))

for (ctx in context_list) {
  for (pheno in pheno_list) {
    matched       <- names(pheno_patterns)[sapply(names(pheno_patterns), grepl, pheno)]
    phenoNoDigits <- if (length(matched)) pheno_patterns[matched[1]] else if (grepl("_(0|2)$", pheno)) pheno else gsub("[0-9]+$", "", pheno)
    phenoLower    <- tolower(phenoNoDigits)
    misc_file <- paste0(misc_dir, "misc_", phenoLower, "_white_euro_", ctx, ".txt")
    if (!file.exists(misc_file)) next
    misc_df <- as.data.frame(t(read.table(misc_file, header = TRUE, stringsAsFactors = FALSE)))
    c_effect_mat[pheno, ctx]    <- as.numeric(misc_df$c_effect)
    pgs_gxc_cor_mat[pheno, ctx] <- as.numeric(misc_df$pgs_gxc_cor)
    pgs_pv_mat[pheno, ctx]      <- as.numeric(misc_df$pgs_pv)
  }
}

for (context in context_list) {
  for (pop in valid_pops) {
    output_loc  <- round(output_valid[pheno_list, , , context, pop], digits=4)
    R2          <- output_loc[, alt_pgs_tbl, "r2"]
    base_r2     <- R2[, "pgs"]
    per_diff_r2 <- scale_by_pgs(output_loc[, alt_pgs_tbl, "delta"], base_r2)
    CI25        <- scale_by_pgs(output_loc[, alt_pgs_tbl, "CI25"],  base_r2)
    CI75        <- scale_by_pgs(output_loc[, alt_pgs_tbl, "CI75"],  base_r2)
    PV          <- output_loc[, alt_pgs_tbl, "pval"]
    rhos <- round(output_scaled[, c('pgs','ampPGS','PGSC'), c('w0','w0a','w0b','w1','w1a','w1b'), context], digits=4)
    taus <- output_scaled[, c('pgs','ampPGS','PGSC'), c('tau0','tau1'), context]

    full_table <- data.frame(
      Phenotype    = pheno_cleaner(pheno_list),
      PGS_R2       = R2[,"pgs"],
      PGS_pv       = PV[,"pgs"],
      PGS_w0     = rhos[,'pgs','w0'],
      PGS_tau0     = taus[,'pgs','tau0'],
      ampPGS_R2    = R2[,"ampPGS"],
      ampPGS_pv    = PV[,"ampPGS"],
      ampPGS_R2pChange = per_diff_r2[,"ampPGS"],
      ampPGS_CI25  = CI25[,"ampPGS"],
      ampPGS_CI75  = CI75[,"ampPGS"],
      ampPGS_w0a = rhos[,'ampPGS','w0a'],
      ampPGS_w0b = rhos[,'ampPGS','w0b'],
      ampPGS_rho = rhos[,'ampPGS','w0b']/rhos[,'ampPGS','w0a'],
      ampPGS_tau0  = taus[,'ampPGS','tau0'],
      PGSC_R2      = R2[,"PGSC_v_pgs"],
      PGSC_pv      = PV[,"PGSC_v_pgs"],
      PGSC_R2pChange   = per_diff_r2[,"PGSC_v_pgs"],
      PGSC_CI25    = CI25[,"PGSC_v_pgs"],
      PGSC_CI75    = CI75[,"PGSC_v_pgs"],
      PGSC_w0    = rhos[,'PGSC','w0'],
      PGSC_w1a   = rhos[,'PGSC','w1a'],
      PGSC_w1b   = rhos[,'PGSC','w1b'],
      PGSC_rho1a = rhos[,'PGSC','w0']/rhos[,'PGSC','w1a'],
      PGSC_rho1b = rhos[,'PGSC','w0']/rhos[,'PGSC','w1b'],
      PGSC_tau0    = taus[,'PGSC','tau0'],
      PGSC_tau1    = taus[,'PGSC','tau1'],
      row.names = NULL, stringsAsFactors = FALSE
    )

    if (pop == "white_euro") {
      full_table$c_effect    <- c_effect_mat[pheno_list, context]
      full_table$pgs_gxc_cor <- pgs_gxc_cor_mat[pheno_list, context]
      full_table$PGS_pv_ols  <- pgs_pv_mat[pheno_list, context]
    }

    write.csv(full_table,
              file      = paste0("summ_stats_", context, "_", pop, ".csv"),
              row.names = FALSE)
  }
}
