library(ggplot2)
library(ggrepel)

# ── Shared parameters ─────────────────────────────────────────────────────────
FC_THRESH   <- 1.0      # |log2FC| cutoff
PADJ_THRESH <- 0.05     # adjusted p-value cutoff
TOP_N_LABEL <- 20       # number of top genes to label (by -log10 padj)

make_volcano <- function(df, title,
                         fc_thresh   = FC_THRESH,
                         padj_thresh = PADJ_THRESH,
                         top_n       = TOP_N_LABEL) {

  # Gene names from rownames
  df$gene <- rownames(df)

  # Guard against zeros / NAs before log
  df$p_val_adj[df$p_val_adj == 0] <- .Machine$double.xmin
  df$log10_padj <- -log10(df$p_val_adj)

  # Significance category
  df$sig <- "NS"
  df$sig[df$p_val_adj < padj_thresh & df$avg_log2FC >  fc_thresh] <- "Up"
  df$sig[df$p_val_adj < padj_thresh & df$avg_log2FC < -fc_thresh] <- "Down"
  df$sig <- factor(df$sig, levels = c("Up", "Down", "NS"))

  colours <- c("Up" = "#e74c3c", "Down" = "#3498db", "NS" = "grey70")

  # Top genes to label: highest -log10(padj) among significant
  sig_df    <- df[df$sig != "NS", ]
  label_df  <- sig_df[order(sig_df$log10_padj, decreasing = TRUE), ]
  label_df  <- head(label_df, top_n)

  # Count summary for subtitle
  n_up   <- sum(df$sig == "Up")
  n_down <- sum(df$sig == "Down")
  subtitle <- sprintf("Up: %d  |  Down: %d  (|log2FC| > %.1f, padj < %.2f)",
                      n_up, n_down, fc_thresh, padj_thresh)

  ggplot(df, aes(x = avg_log2FC, y = log10_padj, colour = sig)) +
    geom_point(size = 0.9, alpha = 0.65) +
    geom_vline(xintercept = c(-fc_thresh, fc_thresh),
               linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = -log10(padj_thresh),
               linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_text_repel(data      = label_df,
                    aes(label = gene),
                    size      = 2.8,
                    max.overlaps = 30,
                    segment.colour = "grey50",
                    segment.size   = 0.3,
                    show.legend    = FALSE) +
    scale_colour_manual(values = colours,
                        name   = NULL,
                        guide  = guide_legend(override.aes = list(size = 3,
                                                                   alpha = 1))) +
    labs(title    = title,
         subtitle = subtitle,
         x        = expression(log[2]~"Fold Change"),
         y        = expression(-log[10]~"(adj. p-value)")) +
    theme_classic(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 9, colour = "grey40"),
      legend.position = "top"
    )
}

# ── Plots ─────────────────────────────────────────────────────────────────────
p_immune <- make_volcano(degs_immune, "Immune — Differentially Expressed Genes")
p_endo   <- make_volcano(degs_endo,   "Endothelial — Differentially Expressed Genes")

ggsave("volcano_immune.png", p_immune, width = 7, height = 6, dpi = 200)
ggsave("volcano_endo.png",   p_endo,   width = 7, height = 6, dpi = 200)

message("Saved: volcano_immune.png  volcano_endo.png")
