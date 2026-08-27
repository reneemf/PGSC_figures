rm(list=ls())
source(file.path(Sys.getenv("PGSC_HOME", unset = "."), "config.R"))
load(paste0(rdata_dir, "compiled_output.Rdata"))

setwd(paste0(fig_dir, "ED_Fig2"))

genie <- read.csv(paste0(home_dir,'GENIE/ali_ajhg_results.csv'), head=T )

match_list <- match_list[!names(match_list) %in% c(pheno_drop, pheno_sig_drop)]
match_list_cleaned <- setNames(match_list, sapply(names(match_list), pheno_cleaner))
pheno_match_vals <- unname(match_list_cleaned)

# ── Panels (a)-(c): PGSC R2s vs opt rhos ────────────────────────────────────
pheno_list <- plot_phenos
alt_pgs    <- alt_pgs3              # pgs, ampPGS, PGSC_v_pgs
base_cols  <- base_cols3
names(base_cols) <- alt_pgs
rho_a_col  <- '#E88BB6'             # color for the "weighted towards A" points

r2 <- array(NA,dim=c(length(pheno_list), length(context_list)),
            dimnames = list(pheno_list, context_list) )
per_diff_r2 <- array(NA,dim=c(length(pheno_list), length(context_list), length(alt_pgs[2:3])),
                     dimnames = list(pheno_list, context_list, alt_pgs[2:3] ) )
ws <- array(NA,dim=c(length(pheno_list), length(context_list),length(alt_pgs),6),
            dimnames = list(pheno_list, context_list,alt_pgs,c("w0","w0a","w0b","w1","w1a","w1b")))

for(context in context_list){
  for(pheno in pheno_list){
    tryCatch({
      r2[pheno,context ] <- output_valid[pheno,'pgs','r2',context,'white_euro']
      per_diff_r2[pheno,context,alt_pgs[2:3]] <- scale_by_pgs(output_valid[pheno, alt_pgs[2:3], 'delta',context, 'white_euro'],
                                                              output_valid[pheno, 'pgs', 'r2',context, 'white_euro'])
      tau0 <- as.character(output_valid[pheno,'pgs', 'thresh',context,"white_euro"])
      tau1 <- as.character(output_valid[pheno,"PGSC_v_pgs", 'thresh',context,"white_euro"])
      ws[pheno,context,alt_pgs,"w0"] <- output_scaled[pheno, c('pgs','ampPGS','PGSC'),'w0', context]
      ws[pheno,context,"ampPGS",c('w0a','w0b')] <- output_scaled[pheno, 'ampPGS', c('w0a','w0b'), context]
      ws[pheno,context,"PGSC_v_pgs",c('w1a','w1b')] <- output_scaled[pheno, 'PGSC', c('w1a','w1b'), context]
    }, error=function(e) print( c(context, pheno ) ) )
  }
}

w1a_labels <- c("males", "born 1952-1971", "statin users")
w1b_labels <- c("females", "born 1934-1951", "statin non-users")
arrow_labels_a <- c(sex = "Weighted towards\nmale score",    age = "Weighted towards\nyounger score",    statins = "Weighted towards\nstatin user score")
arrow_labels_b <- c(sex = "Weighted towards\nfemale score",  age = "Weighted towards\nolder score",      statins = "Weighted towards\nnon-user score")
x_label_frac   <- c(sex = 0.57, age = 0.35, statins = 0.57)

# ── Panels (d)-(i): Rhos vs GxC heritability & heteroskedasticity (GENIE) ───
base_cols_genie <- c(base_cols3, rho_a_col)
names(base_cols_genie) <- c(alt_pgs3,"1b")

taus_genie <- array(NA,dim=c(length(pheno_list), length(context_list),2),
              dimnames = list(pheno_list, context_list, c('tau0','tau1')))
diff_r2_genie <- array(NA,dim=c(length(pheno_list), length(context_list), length(alt_pgs3)),
                     dimnames = list(pheno_list, context_list, alt_pgs3 ) )
per_diff_r2_genie <- array(NA,dim=c(length(pheno_list), length(context_list), length(alt_pgs3)),
                     dimnames = list(pheno_list, context_list, alt_pgs3 ) )
ws_genie <- array(NA,dim=c(length(pheno_list), length(context_list),length(alt_pgs3),5),
            dimnames = list(pheno_list, context_list,alt_pgs3,c("w0","w0a","w0b","w1a","w1b")))
merged_list <- list()

for(context in context_list){
  for(pheno in pheno_list){
      diff_r2_genie[pheno,context,alt_pgs3[2:3]] <- output_valid[pheno, alt_pgs3[2:3], 'delta',context, valid_pops[1]]
      per_diff_r2_genie[pheno,context,alt_pgs3[2:3]] <- scale_by_pgs(output_valid[pheno, alt_pgs3[2:3], 'delta',context, valid_pops[1]],
                                                              output_valid[pheno, 'pgs', 'r2',context, valid_pops[1]])
      ws_genie[pheno,context,alt_pgs3,"w0"] <- output_scaled[pheno, c('pgs','ampPGS','PGSC'),'w0', context]
      ws_genie[pheno,context,"ampPGS",c('w0a','w0b')] <- output_scaled[pheno, 'ampPGS',c('w0a','w0b'), context]
      ws_genie[pheno,context,"PGSC_v_pgs",c('w1a','w1b')] <- output_scaled[pheno, 'PGSC',c('w1a','w1b'), context]
      taus_genie[pheno,context,'tau1'] <- output_valid[pheno, alt_pgs3[3], 'thresh', context, valid_pops[1]]
  }

  # set up genie & link to R2 diffs/ws
  g_con <- if (context == "statins") "statin" else context
  genie_context <- genie[genie[,'env'] == g_con,]
  R2s <- as.data.frame(cbind(diff_r2_genie[pheno_list,context,alt_pgs3[2:3]],per_diff_r2_genie[pheno_list,context,alt_pgs3[2:3]]))
  colnames(R2s) <- c(paste0("diff_",alt_pgs3[2:3]), alt_pgs3[2:3])
  PGSC_ws <- as.data.frame(ws_genie[pheno_list,context,"PGSC_v_pgs",c("w0",'w1a','w1b')])
  R2s$clean_phenos <- pheno_cleaner(rownames(R2s))
  R2s$pheno_match <- pheno_match_vals
  PGSC_ws$pheno_match <- pheno_match_vals
  tau1_df <- data.frame(tau1 = taus_genie[pheno_list, context, 'tau1'], pheno_match = pheno_match_vals)
  merged_df <- na.omit(merge(
    genie_context[, c('pheno','h2g','h2gxe','h2nxe')],
    cbind(R2s, PGSC_ws[, !names(PGSC_ws) %in% "pheno_match"],
          tau1_df[, !names(tau1_df) %in% "pheno_match"]),
    by.x = "pheno", by.y = "pheno_match", all.x = TRUE
  ))
  merged_list[[context]] <- merged_df
}

# ── Combined 3x3 figure ──────────────────────────────────────────────────────
panel_labels <- paste0("(", letters[1:9], ")")

# panel plot-box size (fixed, unchanged from previous layout: 2.208in x 2.076in)
box_h        <- 2.076
mai_left     <- 0.7
mai_right    <- 0.3
mai_bottom   <- 0.6         # ~0.76cm gap between each panel row and its axis label (tightest that avoids tick/title overlap)
mai_top_titl <- 0.5        # top margin for every row (row 1: context title; rows 2-3: cor/p stat)

row_heights <- c(box_h + mai_top_titl + mai_bottom,
                 box_h + mai_top_titl + mai_bottom,
                 box_h + mai_top_titl + mai_bottom)

# small helpers shared by every panel, to avoid repeating the same calls 9x
panel_margins <- function(top) par(mai = c(mai_bottom, mai_left, top, mai_right))
panel_tag     <- function(label) mtext(label, side = 3, adj = -0.01, line = 0.15, cex = 0.7, font = 1)
cor_legend    <- function(x, y) mtext(cor_label("", cor.test(x, y, method = "pearson")), side = 3, line = 0.75, cex = 0.8)

# nudge the extreme-rho phenotype labels in panels (a)-(c) apart vertically so they don't overlap
resolve_label_overlap <- function(y_true, y_span, min_gap, y_lo, y_hi) {
  y_placed <- pmin(pmax(y_true, y_lo), y_hi)
  for (iter in 1:500) {
    moved <- FALSE
    y_target <- ifelse(y_true > 0, y_true + 0.03 * y_span, y_true)
    y_placed <- y_placed + 0.05 * (y_target - y_placed)
    for (i in 2:length(y_placed)) {
      if (y_placed[i-1] - y_placed[i] < min_gap) {
        mid <- (y_placed[i-1] + y_placed[i]) / 2
        y_placed[i-1] <- mid + min_gap/2;  y_placed[i] <- mid - min_gap/2
        moved <- TRUE
      }
    }
    for (i in (length(y_placed)-1):1) {
      if (y_placed[i] - y_placed[i+1] < min_gap) {
        mid <- (y_placed[i] + y_placed[i+1]) / 2
        y_placed[i] <- mid + min_gap/2;  y_placed[i+1] <- mid - min_gap/2
        moved <- TRUE
      }
    }
    y_placed <- pmin(pmax(y_placed, y_lo), y_hi)
    if (!moved) break
  }
  y_placed
}

# arrows + labels in the left margin of panels (a)-(c) showing which end of rho_c means what
draw_score_arrows <- function(usr, x_span, y_span, context) {
  x_arrow  <- usr[1] - 0.16 * x_span;  x_text <- usr[1] - 0.2 * x_span
  y_bottom <- usr[3] + 0.3  * y_span;  y_top  <- usr[4] - 0.3 * y_span
  arrows(x0 = x_arrow, y0 = y_bottom, x1 = x_arrow, y1 = y_top,
         code = 3, length = 0.08, lwd = 2, col = "black", xpd = TRUE)
  text(x_text, y_top    + 0.04 * y_span, arrow_labels_a[context], srt = 90, adj = 0, xpd = TRUE, cex = 0.9)
  text(x_text, y_bottom - 0.04 * y_span, arrow_labels_b[context], srt = 90, adj = 1, xpd = TRUE, cex = 0.9)
}

png("R2_vs_rhos_and_Rhos_combined.png", width = 9, height = sum(row_heights), units = 'in', res = 300)
layout(matrix(1:9, nrow = 3, byrow = TRUE), heights = row_heights)
# layout() (unlike mfrow) does not auto-shrink text for a 3x3 grid, so restore
# the cex mfrow would have applied to match the original panel proportions
par(mgp = c(3, 0.7, 0), cex = 0.66, lwd = 1.5)

# Row 1, panels (a)-(c): PGSC R2s vs opt rhos
for (c in seq_along(context_list)) {
  panel_margins(mai_top_titl)
  context <- context_list[c]
  y1     <- ws[,context,"PGSC_v_pgs","w1a"] / ws[,context,"PGSC_v_pgs","w0"]
  y2     <- ws[,context,"PGSC_v_pgs","w1b"] / ws[,context,"PGSC_v_pgs","w0"]
  x_vals <- per_diff_r2[,context,"PGSC_v_pgs"]

  y_range  <- range(c(y1, y2), na.rm = TRUE)
  ylim_use <- y_range + diff(y_range) * c(-0.05, 0.05)

  plot(x_vals, y1,
       xlim = range(x_vals, na.rm=TRUE) + c(0, 8), ylim = ylim_use,
       xlab="", ylab="", col = rho_a_col, pch = 16)
  segments(x0=x_vals, y0=y1, x1=x_vals, y1=y2, col="black")
  title(xlab = "R\u00b2 %Change", cex.lab = 1.5, line = 2.5)
  title(ylab = ifelse(context == "sex", expression(rho[c]), " "), cex.lab = 1.5, line = 3)
  title(main = context_mapping[[context]]$clean_context, cex.main = 1.5, font.main = 1, line = 1.5)
  points(x_vals, y1, col=rho_a_col, pch=16)
  points(x_vals, y2, col=base_cols["PGSC_v_pgs"], pch=16)
  abline(v=0); abline(h=0)
  panel_tag(panel_labels[c])
  legend("topleft", inset=c(ifelse(context=="sex", 0, 0.035),0.05),
         legend=c(str_to_sentence(w1a_labels[c]),str_to_sentence(w1b_labels[c])), bty = "n",
         col=c(rho_a_col, base_cols["PGSC_v_pgs"]), pch=16, cex=1)

  # --- labels for the most extreme rhos ---
  rho_extreme <- ifelse(abs(y1) >= abs(y2), y1, y2)
  label_idx   <- c(order(rho_extreme, decreasing=TRUE)[1:4],
                   order(rho_extreme, decreasing=FALSE)[1:4])
  usr    <- par("usr")
  x_span <- usr[2] - usr[1];  y_span <- usr[4] - usr[3]
  x_label <- usr[2] - x_label_frac[context] * x_span

  ord           <- order(rho_extreme[label_idx], decreasing=TRUE)
  y_true_sorted <- rho_extreme[label_idx][ord]
  labels_sorted <- pheno_cleaner(rownames(per_diff_r2)[label_idx])[ord]
  x_points      <- x_vals[label_idx][ord]

  y_placed <- resolve_label_overlap(y_true_sorted, y_span,
                                     min_gap = 0.09 * y_span,
                                     y_lo = usr[3] + 0.02 * y_span,
                                     y_hi = usr[4] - 0.02 * y_span)

  segments(x0=x_points, y0=y_true_sorted, x1=x_label, y1=y_placed, col=detail_col["grey_2"], lwd=0.8)
  text(x_label, y_placed, labels=labels_sorted, pos=4, cex=0.75)

  draw_score_arrows(usr, x_span, y_span, context)
}

# Rows 2-3, panels (d)-(i): PGSC rhos vs GxC heritability & heteroskedasticity (GENIE)
row_specs <- list(
  list(get_x = function(df) df$h2gxe,
       get_y = function(df) df$w1a/df$w0 - df$w1b/df$w0,
       xlab = "GxC heritability", ylab_prefix = "Diff", offset = 3),
  list(get_x = function(df) abs(df$h2nxe),
       get_y = function(df) abs(rowMeans(cbind(df$w1a/df$w0, df$w1b/df$w0))),
       xlab = "Heteroskedasticity", ylab_prefix = "Mean", offset = 6)
)

for (spec in row_specs) {
  for (i in seq_along(context_list)) {
    context   <- context_list[i]
    merged_df <- merged_list[[context]]
    x <- spec$get_x(merged_df)
    y <- spec$get_y(merged_df)

    panel_margins(mai_top_titl)
    plot(x, y, xlim = range(x), ylim = range(y), xlab = "", ylab = "",
         col = base_cols_genie["PGSC_v_pgs"], pch = 16)
    title(xlab = spec$xlab, cex.lab = 1.5, line = 2.5)
    title(ylab = bquote(.(paste0(spec$ylab_prefix, "(Gx", str_to_sentence(context), ") ")) * rho[c]),
          cex.lab = 1.5, line = 2.5)
    abline(h=0, col="black"); abline(v=0, col="black")
    panel_tag(panel_labels[spec$offset + i])
    cor_legend(x, y)
  }
}

dev.off()
