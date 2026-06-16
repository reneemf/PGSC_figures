rm(list=ls())
source(file.path(Sys.getenv("PGSC_HOME", unset = "."), "config.R"))
load(paste0(rdata_dir, "compiled_output.Rdata"))

setwd(paste0(fig_dir, "ED_Fig2"))

# for w1a w1b label axis as male v female or old v young
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

# PGSC R2s vs opt rhos
w1a_labels <- c("males", "born 1952-1971", "statin users")
w1b_labels <- c("females", "born 1934-1951", "statin non-users")
arrow_labels_a <- c(sex = "Weighted towards\nmale score",    age = "Weighted towards\nyounger score",    statins = "Weighted towards\nstatin user score")
arrow_labels_b <- c(sex = "Weighted towards\nfemale score",  age = "Weighted towards\nolder score",      statins = "Weighted towards\nnon-user score")
x_label_frac   <- c(sex = 0.57, age = 0.35, statins = 0.57)

png("R2_vs_rhos.png", width = 9, height = 3, units = 'in', res = 300)
par(mfrow = c(1,3), mar = c(4, 5, 2, 1), oma = c(0, 0, 0, 0))
panel_labels <- paste0("(", LETTERS[seq_along(context_list)], ")")
for (c in seq_along(context_list)) {
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
  title(xlab = "R\u00b2 %Change", cex.lab = 1.25, line = 2.5)
  title(ylab = ifelse(context == "sex", expression(rho[c]), " "), cex.lab = 1.25, line = 3)
  points(x_vals, y1, col=rho_a_col, pch=16)
  points(x_vals, y2, col=base_cols["PGSC_v_pgs"], pch=16)
  abline(v=0); abline(h=0)
  mtext(panel_labels[c], side=3, adj=0.008, line = -0.9, cex=0.6, font=1)
  legend("topleft", inset=c(ifelse(context=="sex", 0, 0.035),0.05),
         legend=c(str_to_sentence(w1a_labels[c]),str_to_sentence(w1b_labels[c])), bty = "n",
         col=c(rho_a_col, base_cols["PGSC_v_pgs"]), pch=16, cex=1)
  # labels for most extreme rhos
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

  min_gap  <- 0.09 * y_span
  y_lo     <- usr[3] + 0.02 * y_span
  y_hi     <- usr[4] - 0.02 * y_span
  y_placed <- pmin(pmax(y_true_sorted, y_lo), y_hi)

  for(iter in 1:500){
    moved    <- FALSE
    y_target_biased <- ifelse(y_true_sorted > 0, y_true_sorted + 0.03 * y_span, y_true_sorted)
    y_placed <- y_placed + 0.05 * (y_target_biased - y_placed)
    for(i in 2:length(y_placed)){
      if(y_placed[i-1] - y_placed[i] < min_gap){
        mid <- (y_placed[i-1] + y_placed[i]) / 2
        y_placed[i-1] <- mid + min_gap/2;  y_placed[i] <- mid - min_gap/2
        moved <- TRUE
      }
    }
    for(i in (length(y_placed)-1):1){
      if(y_placed[i] - y_placed[i+1] < min_gap){
        mid <- (y_placed[i] + y_placed[i+1]) / 2
        y_placed[i] <- mid + min_gap/2;  y_placed[i+1] <- mid - min_gap/2
        moved <- TRUE
      }
    }
    y_placed <- pmin(pmax(y_placed, y_lo), y_hi)
    if(!moved) break
  }

  segments(x0=x_points, y0=y_true_sorted, x1=x_label, y1=y_placed, col=detail_col["grey_2"], lwd=0.8)
  text(x_label, y_placed, labels=labels_sorted, pos=4, cex=0.75)

  # y axis arrows
  x_arrow <- usr[1] - 0.16 * x_span;  x_text <- usr[1] - 0.2 * x_span
  y_bottom <- usr[3] + 0.3 * y_span;  y_top   <- usr[4] - 0.3 * y_span
  arrows(x0=x_arrow, y0=y_bottom, x1=x_arrow, y1=y_top,
         code=3, length=0.08, lwd=2, col="black", xpd=TRUE)
  text(x_text, y_top    + 0.04*y_span, arrow_labels_a[context], srt=90, adj=0, xpd=TRUE, cex=0.9)
  text(x_text, y_bottom - 0.04*y_span, arrow_labels_b[context], srt=90, adj=1, xpd=TRUE, cex=0.9)
}

dev.off()







