##script for rachel
singlet_obj <- readRDS("integrated-singlet-obj.RDS")







DotPlot(cm_subset, features = c("Myh6", "Tnni3", "Myl3", "Myl2", "Myh7", "Tnni1", "Myl4", "Myl7" ))



# ============================================================
# Cardiomyocyte-only DotPlot, split by condition,
# genes grouped/colored as Embryonic vs Mature
# ============================================================

library(Seurat)
library(ggplot2)

# ---- 0. EDIT THESE ------------------------------------------------------
seu           <- seu                 # your Seurat object
celltype_col  <- "cell_type"         # meta.data column with annotations
cm_label      <- "Cardiomyocytes"    # value(s) identifying CMs
condition_col <- "condition"         # meta.data column with condition
cond_levels   <- NULL                # e.g. c("WT", "Mutant"); NULL = keep as-is
# -------------------------------------------------------------------------

# ---- 1. Subset to cardiomyocytes ----------------------------------------
Idents(seu) <- celltype_col
cm <- subset(seu, idents = cm_label)

DefaultAssay(cm) <- "RNA"            # plot normalized RNA, not SCT/integrated
# cm <- NormalizeData(cm)            # only if RNA data slot isn't populated

if (!is.null(cond_levels)) {
  cm@meta.data[[condition_col]] <- factor(cm@meta.data[[condition_col]],
                                          levels = cond_levels)
}
Idents(cm) <- condition_col

gene_groups <- list(
  Mature    = c("Myh6", "Tnni3", "Myl3", "Myl2"),
  Embryonic = c("Myh7", "Tnni1", "Myl4", "Myl7")
)

# drop anything missing from the object so DotPlot doesn't error
gene_groups <- lapply(gene_groups, function(g) g[g %in% rownames(cm_subset)])
stopifnot(lengths(gene_groups) > 0)

group_cols <- c(Mature = "#B2182B", Embryonic = "#2166AC")

# ---- 3. Build the dot plot ----------------------------------------------
# Pass a flat vector — no list, no faceting. (A named list would trip
# Seurat's internal facet_grid(facets = ~...), defunct in ggplot2 >= 3.5.)
feat_vec <- unlist(gene_groups, use.names = FALSE)
grp_map  <- setNames(rep(names(gene_groups), lengths(gene_groups)), feat_vec)

p <- DotPlot(
  cm_subset,
  features   = feat_vec,
  group.by   = "condition",
  cols       = c("lightgrey", "#08306B"),
  dot.scale  = 9,
  scale      = TRUE          # see note in step 5 if you have few conditions
) +
  labs(x = NULL, y = "Condition") +
  theme(
    axis.text.x  = element_blank(),   # labels move to the annotation strip
    axis.ticks.x = element_blank(),
    axis.text.y  = element_text(size = 11),
    panel.border = element_rect(colour = "grey70", fill = NA),
    plot.margin  = margin(5, 5, 0, 5)
  )

# ---- 4. Color bar annotating the gene groups ----------------------------
# Same discrete x scale and default expansion as the dot plot, so patchwork
# aligns the two panels column-for-column.
bar_df <- data.frame(
  gene  = factor(feat_vec, levels = feat_vec),
  group = factor(unname(grp_map[feat_vec]), levels = names(gene_groups))
)

p_bar <- ggplot(bar_df, aes(x = gene, y = 1, fill = group)) +
  geom_tile(width = 0.95, height = 1) +
  scale_fill_manual(values = group_cols, name = NULL) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x     = element_text(face = "italic", angle = 45, hjust = 1,
                                   size = 11, colour = "black"),
    axis.text.y     = element_blank(),
    panel.grid      = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(9, "pt"),
    plot.margin     = margin(0, 5, 5, 5)
  )

# ---- 4b. Stack them -----------------------------------------------------
library(patchwork)
final <- p / p_bar +
  plot_layout(heights = c(1, 0.10), guides = "keep")

print(final)

# ---- 5. Notes -----------------------------------------------------------
# * `scale = TRUE` z-scores each gene ACROSS THE GROUPS SHOWN. With only two
#   conditions that collapses to +/- one value per gene, which looks dramatic
#   but carries no magnitude information. With 2-3 conditions prefer:
#       DotPlot(cm, features = gene_groups, scale = FALSE) +
#         scale_color_viridis_c(name = "Mean expr\n(log-norm)")
# * Dot color = expression level, dot size = % of CMs expressing. The
#   embryonic/mature distinction is carried by the facets and label colors.
# * To also split by a second variable (e.g. sex or timepoint), use
#   group.by = "condition" and split.by = "sex" — note split.by switches the
#   color scale to one hue per split level.

# ---- 6. Save ------------------------------------------------------------
ggsave("cm_dotplot_embryonic_mature.pdf", final, width = 7, height = 3.4)
