library(clValid)
library(openxlsx)  # for saving Excel files

# --- Set working directory ---
setwd("C:/Users/ga53hil/Desktop/Granit_proteomics/13.02.26_Scientific_Data")

# --- Load and preprocess data ---
data <- read.table("RAW_P337_02B_proteinGroups.txt", 
header = TRUE, sep = "\t", quote = "", check.names = FALSE)

# --- Filter out unwanted rows ---
data_filtered <- data[
data$`Only identified by site` != "+" &
data$Reverse != "+" &
data$`Potential contaminant` != "+",
]

# --- Reset row names ---
row.names(data_filtered) <- NULL

# --- Save filtered data (TXT) ---
write.table(data_filtered,
"P337_02B_proteinGroups_filtered.txt",
sep = "\t", quote = FALSE, row.names = FALSE)

# --- Save filtered data (Excel) ---
write.xlsx(data_filtered,
"P337_02B_proteinGroups_filtered.xlsx",
overwrite = TRUE)


################################################################################
################################################################################
################################################################################

# --- Keep proteins with Peptides > 1 ---
# Coerce safely in case it's read as character
data_filtered$Peptides_num <- suppressWarnings(as.numeric(data_filtered$Peptides))

data_filtered2 <- subset(
data_filtered,
!is.na(Peptides_num) & Peptides_num > 1
)

# Drop helper column
data_filtered2$Peptides_num <- NULL

# --- Save result (TXT) ---
write.table(data_filtered2,
"P337_02B_proteinGroups_filtered_pepGT1.txt",
sep = "\t", quote = FALSE, row.names = FALSE)

# --- Save result (Excel) ---
library(openxlsx)
write.xlsx(data_filtered2,
"P337_02B_proteinGroups_filtered_pepGT1.xlsx",
overwrite = TRUE)


################################################################################
################################################################################
################################################################################

# --- Identify LFQ intensity columns (original names) ---
lfq_cols <- grep("^LFQ intensity\\s", colnames(data_filtered2), value = TRUE)

# --- Copy before transforming ---
log2_data_filtered2 <- data_filtered2

# --- Robust log2 transform of LFQ columns ---
log2_data_filtered2[lfq_cols] <- lapply(data_filtered2[lfq_cols], function(x) {
# 1) force to character, strip thousands separators, then to numeric
x_num <- suppressWarnings(as.numeric(gsub(",", "", as.character(x))))
# 2) treat zeros/negatives as missing (will be imputed later)
x_num[x_num <= 0] <- NA_real_
# 3) log2 and clean non-finite
y <- log2(x_num)
y[!is.finite(y)] <- NA_real_
y
})

# --- Save (TXT) ---
write.table(log2_data_filtered2,
"P337_02B_proteinGroups_filtered_pepGT1_log2LFQ.txt",
sep = "\t", quote = FALSE, row.names = FALSE)

# --- Save (Excel) ---
library(openxlsx)
write.xlsx(log2_data_filtered2,
"P337_02B_proteinGroups_filtered_pepGT1_log2LFQ.xlsx",
overwrite = TRUE)

################################################################################
################################################################################
################################################################################

################################################################################
# Add Gene_name via left join (Protein IDs ↔ UniProt Entry)
# Uses: log2_data_filtered2
################################################################################

# --- Load UniProt mapping (Entry -> Gene Names) ---
uni_file <- "C:/Users/ga53hil/Desktop/Granit_proteomics/13.02.26_Scientific_Data/FASTA_uniprotkb_taxonomy_id_9913_Bovine.tsv"
uniprot <- read.table(
uni_file,
header = TRUE, sep = "\t", quote = "",
check.names = FALSE, stringsAsFactors = FALSE,
comment.char = "", fill = TRUE
)

# Safety checks
if (!"Entry" %in% colnames(uniprot))
stop("Column 'Entry' not found in UniProt TSV.")
if (!"Gene Names" %in% colnames(uniprot))
stop("Column 'Gene Names' not found in UniProt TSV.")
if (!"Protein IDs" %in% colnames(log2_data_filtered2))
stop("Column 'Protein IDs' not found in log2_data_filtered2.")

# Keep only needed columns and drop empty gene-name rows
uni_map <- uniprot[, c("Entry", "Gene Names")]
uni_map$`Gene Names`[uni_map$`Gene Names` == ""] <- NA
uni_map <- uni_map[!is.na(uni_map$`Gene Names`), ]
uni_map <- uni_map[!duplicated(uni_map$Entry), ]

# Named lookup: Entry -> Gene Names
map <- setNames(uni_map$`Gene Names`, uni_map$Entry)

# Map "P12345;Q9XXXX" -> "GENE1;GENE2"
map_gene_names <- function(prot_ids) {
if (is.na(prot_ids) || prot_ids == "") return(NA_character_)
ids <- trimws(strsplit(prot_ids, ";", fixed = TRUE)[[1]])
genes <- unname(map[ids])
genes <- genes[!is.na(genes) & genes != ""]
if (!length(genes)) return(NA_character_)
paste(unique(genes), collapse = ";")
}

# --- Add Gene_name column ---
log2_data_filtered2$Gene_name <- vapply(
log2_data_filtered2$`Protein IDs`,
map_gene_names,
character(1)
)

# --- Save (TXT + Excel) ---
write.table(
log2_data_filtered2,
"P337_02B_proteinGroups_filtered_pepGT1_log2LFQ_gene_labeled.txt",
sep = "\t", quote = FALSE, row.names = FALSE
)

library(openxlsx)
write.xlsx(
log2_data_filtered2,
"P337_02B_proteinGroups_filtered_pepGT1_log2LFQ_gene_labeled.xlsx",
overwrite = TRUE
)


################################################################################
################################################################################
################################################################################

# Mapping from Roman numerals to numbers
roman_to_num <- c(
"I" = 1, "II" = 2, "III" = 3, "IV" = 4, "V" = 5,
"VI" = 6, "VII" = 7, "VIII" = 8, "IX" = 9, "X" = 10
)

# Loop through LFQ intensity columns and rename
new_colnames <- colnames(log2_data_filtered2)
for (i in seq_along(new_colnames)) {
if (grepl("^LFQ intensity", new_colnames[i])) {
# Extract the Roman numeral and replicate number
parts <- unlist(strsplit(gsub("^LFQ intensity ", "", new_colnames[i]), "-"))
roman <- trimws(parts[1])
rep <- trimws(parts[2])
# Map to new name format
new_colnames[i] <- paste0("T", roman_to_num[roman], "_", rep)
}
}

# Apply new column names
colnames(log2_data_filtered2) <- new_colnames

################################################################################
################################################################################
################################################################################

# --- Identify LFQ (T#_#) columns ---
lfq_cols <- grep("^T\\d+_\\d+$", colnames(log2_data_filtered2), value = TRUE)

# --- Clean NaN/Inf -> NA in LFQ columns ---
log2_data_filtered2[lfq_cols] <- lapply(log2_data_filtered2[lfq_cols], function(x) {
x[!is.finite(x)] <- NA_real_
x
})

# --- Build group list (T1..T10 detected from column names) ---
groups <- sort(unique(sub("^T(\\d+)_.*", "\\1", lfq_cols)))

# --- Keep rows where at least one group has >=70% valid values ---
passes_any_group <- Reduce(`|`, lapply(groups, function(g) {
cols_g <- grep(paste0("^T", g, "_\\d+$"), colnames(log2_data_filtered2), value = TRUE)
if (length(cols_g) == 0) return(rep(FALSE, nrow(log2_data_filtered2)))
thr <- ceiling(0.70 * length(cols_g))   # e.g., 0.70 * 8 = 5.6 -> 6
rowSums(!is.na(log2_data_filtered2[cols_g])) >= thr
}))

filtered_70 <- log2_data_filtered2[passes_any_group, , drop = FALSE]
row.names(filtered_70) <- NULL

# --- Save result (TXT) ---
write.table(filtered_70,
"P337_02B_proteinGroups_filtered_pepGT1_log2LFQ_valid70.txt",
sep = "\t", quote = FALSE, row.names = FALSE)

# --- Save result (Excel) ---
library(openxlsx)
write.xlsx(filtered_70,
"P337_02B_proteinGroups_filtered_pepGT1_log2LFQ_valid70.xlsx",
overwrite = TRUE)

################################################################################

################################################################################
############################## 6️⃣ CORRELATION QC (FULL) ########################
################################################################################
# Plot only. No removal here.

library(openxlsx)
library(pheatmap)
library(ggplot2)
library(grid)

out_dir <- file.path(getwd(), "Step6_Correlation_QC")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# LFQ columns
lfq_cols <- grep("^T\\d+_\\d+$", colnames(filtered_70), value = TRUE)
if (length(lfq_cols) < 2) stop("Not enough LFQ samples (T#_#) for correlation QC.")

# Numeric safeguard
df_group <- filtered_70[, lfq_cols, drop = FALSE]
df_group <- as.data.frame(lapply(df_group, function(x) suppressWarnings(as.numeric(x))))

# Correlation matrix
corr_mat <- cor(df_group, use = "pairwise.complete.obs", method = "pearson")

# Save matrix (exact values)
write.xlsx(corr_mat, file.path(out_dir, "filtered_70_correlation_matrix.xlsx"), overwrite = TRUE)

# Mean correlation per sample (exclude self-corr)
mean_corr <- apply(corr_mat, 1, function(x) mean(x[!is.na(x) & x < 0.999999], na.rm = TRUE))
corr_summary <- data.frame(
Sample = names(mean_corr),
Mean_Correlation = round(mean_corr, 3),
stringsAsFactors = FALSE
)
write.xlsx(corr_summary, file.path(out_dir, "filtered_70_correlation_summary.xlsx"), overwrite = TRUE)

# ---- Order samples by Time (T1..T10) then replicate number ----
get_time <- function(x) sub("_(.*)$", "", x)
get_rep  <- function(x) as.integer(sub("^T\\d+_", "", x))
time_levels <- paste0("T", 1:10)

ord <- order(factor(get_time(lfq_cols), levels = time_levels), get_rep(lfq_cols))
lfq_cols_ord <- lfq_cols[ord]
corr_ord <- corr_mat[lfq_cols_ord, lfq_cols_ord]

# ---- Palette (KEEP STABLE FOR ALL PLOTS) ----
time_pal <- setNames(
c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e",
"#e6ab02","#a6761d","#1f78b4","#ff7f00","#6a3d9a"),
time_levels
)

# ---- Annotation (Time bar) for heatmap ----
ann <- data.frame(Time = factor(get_time(lfq_cols_ord), levels = time_levels))
rownames(ann) <- lfq_cols_ord

ann_colors <- list(Time = time_pal)

# ---- Gaps between time points (after every 8 samples) ----
gaps <- seq(8, 8*9, by = 8)

heat_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(200)


################################################################################
# (A) OVERVIEW HEATMAP (80 samples) — LOWER TRIANGLE, custom margins (l/r/t/b)
################################################################################

# Keep LOWER triangle only
corr_tri <- corr_ord
corr_tri[upper.tri(corr_tri, diag = TRUE)] <- NA

# ---- hide the "Time" label safely (works across pheatmap versions) ----
# IMPORTANT: This renames ONLY for plot display; palette for (C) remains in time_pal
old_name <- colnames(ann)[1]
colnames(ann)[1] <- " "
names(ann_colors)[names(ann_colors) == old_name] <- " "

# ---- set padding (edit freely) ----
pad_l <- unit(0.5, "cm")
pad_r <- unit(2.0, "cm")
pad_t <- unit(0.5, "cm")
pad_b <- unit(0.5, "cm")

# Build heatmap grob
ph <- pheatmap(
corr_tri,
main = NA,
color = heat_palette,
border_color = NA,
cluster_rows = FALSE,
cluster_cols = FALSE,
show_rownames = FALSE,
show_colnames = FALSE,
annotation_row = ann,
annotation_colors = ann_colors,
gaps_row = gaps,
gaps_col = gaps,
na_col = "white",
fontsize = 32,
fontsize_title = 26,
fontsize_annotation = 32,
silent = TRUE
)

draw_with_padding <- function(g) {
pushViewport(viewport(
x = pad_l,
y = pad_b,
width  = unit(1, "npc") - pad_l - pad_r,
height = unit(1, "npc") - pad_t - pad_b,
just = c("left", "bottom")
))
grid.draw(g)
popViewport()
}

# ---- PDF ----
pdf(
file.path(out_dir, "A_overview_correlation_heatmap_80samples_LOWER_TRIANGLE.pdf"),
width = 12,
height = 11
)
grid.newpage()
draw_with_padding(ph$gtable)
dev.off()

# ---- PNG ----
png(
file.path(out_dir, "A_overview_correlation_heatmap_80samples_LOWER_TRIANGLE.png"),
width = 12,
height = 11,
units = "in",
res = 300
)
grid.newpage()
draw_with_padding(ph$gtable)
dev.off()

# (Optional) restore names if you want ann/ann_colors to stay "Time" downstream
colnames(ann)[1] <- old_name
names(ann_colors)[names(ann_colors) == " "] <- old_name


################################################################################
# (B) PER-TIMEPOINT HEATMAPS (8x8) — readable labels + numbers
################################################################################

for (tp in time_levels) {
cols_tp <- lfq_cols_ord[get_time(lfq_cols_ord) == tp]
if (length(cols_tp) < 2) next

cm_tp <- corr_mat[cols_tp, cols_tp]

pheatmap(
cm_tp,
main = paste0("Within-time correlation: ", tp, " (n=", length(cols_tp), ")"),
color = heat_palette,
border_color = NA,
cluster_rows = TRUE, cluster_cols = TRUE,
display_numbers = TRUE, number_format = "%.2f",
fontsize_row = 9, fontsize_col = 9, fontsize_number = 9,
filename = file.path(out_dir, paste0("B_withinTime_corr_", tp, ".png")),
width = 6, height = 5, dpi = 300
)

pdf(file.path(out_dir, paste0("B_withinTime_corr_", tp, ".pdf")), width = 6, height = 5)
pheatmap(
cm_tp,
main = paste0("Within-time correlation: ", tp, " (n=", length(cols_tp), ")"),
color = heat_palette,
border_color = NA,
cluster_rows = TRUE, cluster_cols = TRUE,
display_numbers = TRUE, number_format = "%.2f",
fontsize_row = 9, fontsize_col = 9, fontsize_number = 9
)
dev.off()
}


################################################################################
# (C) MEAN CORRELATION BARPLOT — ordered T1_1 .. T10_8 (no coord_flip)
################################################################################

# Factor Time
corr_summary$Time <- factor(get_time(corr_summary$Sample), levels = time_levels)

# Exact desired order: T1_1..T1_8, ..., T10_1..T10_8
sample_levels <- as.vector(unlist(lapply(time_levels, function(tp) paste0(tp, "_", 1:8))))
sample_levels <- sample_levels[sample_levels %in% corr_summary$Sample]
corr_summary$Sample <- factor(corr_summary$Sample, levels = sample_levels)

p_bar <- ggplot(corr_summary, aes(x = Sample, y = Mean_Correlation, fill = Time)) +
geom_col(width = 0.5) +

# Threshold line
geom_hline(yintercept = 0.6, linetype = 2) +

# IMPORTANT: use stable palette (time_pal), not ann_colors$Time (can be renamed in A)
scale_fill_manual(values = time_pal, guide = "none") +

scale_y_continuous(
limits = c(0, NA),
expand = c(0, 0)
) +

labs(
title = NULL,
x = NULL,
y = "Mean Pearson correlation"
) +

theme_classic(base_size = 50) +
theme(
plot.title   = element_text(size = 30, face = "bold"),
axis.title.x = element_text(size = 30),
axis.title.y = element_text(size = 30),
axis.text.x  = element_text(size = 20, angle = 90, vjust = 0.5, hjust = 1),
axis.text.y  = element_text(size = 20)
)

ggsave(
file.path(out_dir, "C_meanCorrelation_perSample_ordered_T1_1_to_T10_8.png"),
p_bar, width = 20, height = 12, dpi = 300, bg = "white"
)

ggsave(
file.path(out_dir, "C_meanCorrelation_perSample_ordered_T1_1_to_T10_8.pdf"),
p_bar, width = 20, height = 12
)


################################################################################
############################## 6️⃣ ADD-ON: CV CALCULATION #######################
################################################################################
# Protein-level CV% computed on LINEAR scale (2^(log2 LFQ)).
# Key change: compute TIMEPOINT MEANS first (replicates averaged), then CV across T1..T10.
# Also exports per-timepoint replicate CVs (within each T#).

library(openxlsx)

# Use the same out_dir from Step 6
if (!exists("out_dir")) out_dir <- file.path(getwd(), "Step6_Correlation_QC")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# LFQ columns (T#_#)
lfq_cols <- grep("^T\\d+_\\d+$", colnames(filtered_70), value = TRUE)
if (length(lfq_cols) < 2) stop("Not enough LFQ samples (T#_#) for CV calculation.")

# Ensure numeric
df_lfq <- as.data.frame(lapply(filtered_70[, lfq_cols, drop = FALSE], function(x) {
suppressWarnings(as.numeric(x))
}))

# Convert log2 -> linear intensities
df_lin <- 2^df_lfq

# Helpers
cv_pct <- function(x) {
x <- x[!is.na(x)]
if (length(x) < 2) return(NA_real_)
m <- mean(x)
if (!is.finite(m) || m <= 0) return(NA_real_)
100 * stats::sd(x) / m
}

get_time <- function(x) sub("_(.*)$", "", x)
time_levels <- paste0("T", 1:10)

# ---- (1) CV across TIMEPOINT MEANS (replicates averaged first) ----
tp_mean_mat <- sapply(time_levels, function(tp) {
cols_tp <- lfq_cols[get_time(lfq_cols) == tp]
if (length(cols_tp) == 0) return(rep(NA_real_, nrow(df_lin)))
rowMeans(df_lin[, cols_tp, drop = FALSE], na.rm = TRUE)
})

cv_across_time_means <- apply(tp_mean_mat, 1, cv_pct)

# ---- (2) Per-timepoint CV across replicates (within each T#) ----
cv_by_time <- sapply(time_levels, function(tp) {
cols_tp <- lfq_cols[get_time(lfq_cols) == tp]
if (length(cols_tp) < 2) return(rep(NA_real_, nrow(df_lin)))
apply(df_lin[, cols_tp, drop = FALSE], 1, cv_pct)
})

# ---- Build export table: keep original columns + CV columns ----
cv_table <- filtered_70
cv_table$CV_overall_pct <- round(cv_across_time_means, 2)

for (i in seq_along(time_levels)) {
cv_table[[paste0("CV_", time_levels[i], "_pct")]] <- round(cv_by_time[, i], 2)
}

# Optional: export timepoint means too (useful for debugging/plots)
tp_mean_df <- as.data.frame(tp_mean_mat)
colnames(tp_mean_df) <- paste0("Mean_", time_levels)
tp_mean_df <- cbind(filtered_70[, setdiff(colnames(filtered_70), lfq_cols), drop = FALSE], tp_mean_df)

# Small summary sheet
cv_summary <- data.frame(
Metric = c("Proteins (rows)",
"LFQ samples (cols)",
"Overall CV definition"),
Value  = c(nrow(filtered_70),
length(lfq_cols),
"CV across T1..T10 means (replicates averaged first)"),
stringsAsFactors = FALSE
)

# ---- Export ----
out_file <- file.path(out_dir, "filtered_70_CV.xlsx")
write.xlsx(
list(
CV_table = cv_table,
Timepoint_means = tp_mean_df,
CV_summary = cv_summary
),
out_file,
overwrite = TRUE
)

cat("\n✅ CV outputs saved to: ", out_file, "\n", sep = "")

################################################################################

################################################################################
####################### 6️⃣b REMOVE LOW-CORRELATION SAMPLES (NO RENAME) #########
################################################################################
# IMPORTANT:
# - filtered_70 stays unchanged
# - cor_remove_filter70 is the only updated proteinGroups table you use later
# - cor_remove_filter70_mat is the correlation matrix after removal
################################################################################

library(openxlsx)

# ------------------------- Settings you can change ----------------------------
corr_cutoff <- 0.60        # remove samples with mean correlation < cutoff
min_other_samples <- 3     # if too few comparisons, don't auto-remove (safety)
method <- "pearson"
# -----------------------------------------------------------------------------

# Output directory
if (!exists("out_dir")) out_dir <- file.path(getwd(), "Step6_Correlation_QC")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Identify LFQ columns (ONLY these can be removed)
lfq_cols <- grep("^T\\d+_\\d+$", colnames(filtered_70), value = TRUE)
if (length(lfq_cols) < 2) stop("Not enough LFQ samples (T#_#) for correlation filtering.")

# Numeric-only matrix for correlation (do NOT touch filtered_70)
df_lfq <- filtered_70[, lfq_cols, drop = FALSE]
df_lfq <- as.data.frame(lapply(df_lfq, function(x) suppressWarnings(as.numeric(x))))

# ------------------------- BEFORE correlation matrix --------------------------
corr_mat_before <- cor(df_lfq, use = "pairwise.complete.obs", method = method)
rownames(corr_mat_before) <- lfq_cols
colnames(corr_mat_before) <- lfq_cols

write.xlsx(corr_mat_before,
file.path(out_dir, "REMOVELOW_before_correlation_matrix.xlsx"),
overwrite = TRUE)

# Mean correlation per sample (exclude self-corr)
mean_corr <- apply(corr_mat_before, 1, function(x) {
x2 <- x[is.finite(x) & x < 0.999999]
if (length(x2) < min_other_samples) return(NA_real_)
mean(x2, na.rm = TRUE)
})

decision <- data.frame(
Sample = names(mean_corr),
Mean_Correlation = round(mean_corr, 4),
N_Compared = sapply(names(mean_corr), function(s) {
x <- corr_mat_before[s, ]
sum(is.finite(x) & x < 0.999999)
}),
Remove = ifelse(is.na(mean_corr), FALSE, mean_corr < corr_cutoff),
stringsAsFactors = FALSE
)

drop_samples <- decision$Sample[decision$Remove]

write.xlsx(decision,
file.path(out_dir, "REMOVELOW_decision_table.xlsx"),
overwrite = TRUE)

cat("\n---- Low-correlation removal summary ----\n")
cat("Cutoff:", corr_cutoff, "\n")
cat("LFQ samples evaluated:", nrow(decision), "\n")
cat("LFQ samples flagged for removal:", length(drop_samples), "\n")
if (length(drop_samples)) cat("Flagged:", paste(drop_samples, collapse = ", "), "\n")

# ------------------------- Create UPDATED TABLE: cor_remove_filter70 ----------
# Always start as a copy of filtered_70 (so same columns initially)
cor_remove_filter70 <- filtered_70

if (length(drop_samples) == 0) {

# No removal
cor_remove_filter70_mat <- corr_mat_before

cat("\n✅ No samples removed.\n")
cat("✅ cor_remove_filter70 is an unchanged copy of filtered_70.\n")

} else {

# Remove ONLY the flagged LFQ columns from the COPY
cor_remove_filter70 <- cor_remove_filter70[, setdiff(colnames(cor_remove_filter70), drop_samples), drop = FALSE]

# Remaining LFQ columns in the updated table
lfq_cols_after <- grep("^T\\d+_\\d+$", colnames(cor_remove_filter70), value = TRUE)

# Safety: must keep at least 2 LFQ columns
if (length(lfq_cols_after) < 2) stop("Removal would leave <2 LFQ samples. Aborting.")

# Recompute correlation AFTER removal
df_lfq2 <- cor_remove_filter70[, lfq_cols_after, drop = FALSE]
df_lfq2 <- as.data.frame(lapply(df_lfq2, function(x) suppressWarnings(as.numeric(x))))

cor_remove_filter70_mat <- cor(df_lfq2, use = "pairwise.complete.obs", method = method)
cor_remove_filter70_mat <- as.matrix(cor_remove_filter70_mat)
rownames(cor_remove_filter70_mat) <- lfq_cols_after
colnames(cor_remove_filter70_mat) <- lfq_cols_after

write.xlsx(
data.frame(Removed_Samples = drop_samples, stringsAsFactors = FALSE),
file.path(out_dir, "REMOVELOW_removed_samples.xlsx"),
overwrite = TRUE
)

cat("\n✅ Removed", length(drop_samples), "LFQ sample column(s) from cor_remove_filter70.\n")
cat("✅ filtered_70 remains unchanged.\n")
}

# ------------------------- Export outputs -------------------------------------
write.xlsx(cor_remove_filter70,
file.path(out_dir, "cor_remove_filter70_proteinGroups_table.xlsx"),
overwrite = TRUE)

write.xlsx(cor_remove_filter70_mat,
file.path(out_dir, "REMOVELOW_after_correlation_matrix.xlsx"),
overwrite = TRUE)

cat("\n✅ Outputs saved to:", out_dir, "\n")

################################################################################

#########################Total Peptides per Condition###########################

################################################################################
################################################################################
################################################################################

library(dplyr)
library(ggplot2)

# --- SAME group colors as before ---
group_colors <- c(
"T1"  = "#1b9e77",
"T2"  = "#d95f02",
"T3"  = "#7570b3",
"T4"  = "#e7298a",
"T5"  = "#66a61e",
"T6"  = "#e6ab02",
"T7"  = "#a6761d",
"T8"  = "#1f78b4",
"T9"  = "#ff7f00",
"T10" = "#6a3d9a"
)

groups <- list(
T1  = paste0("Peptides I-",     1:8),
T2  = paste0("Peptides II-",    1:8),
T3  = paste0("Peptides III-",   1:8),
T4  = paste0("Peptides IV-",    1:8),
T5  = paste0("Peptides V-",     1:8),
T6  = paste0("Peptides VI-",    1:8),
T7  = paste0("Peptides VII-",   1:8),
T8  = paste0("Peptides VIII-",  1:8),
T9  = paste0("Peptides IX-",    1:8),
T10 = paste0("Peptides X-",     1:8)
)

all_cols <- unlist(groups, use.names = FALSE)
cor_remove_filter70[all_cols] <- lapply(
cor_remove_filter70[all_cols],
function(x) suppressWarnings(as.numeric(x))
)

# Replicate sums
replicate_sums_mat <- sapply(groups, function(cols) {
colSums(cor_remove_filter70[, cols, drop = FALSE], na.rm = TRUE)
})

df <- data.frame(
Condition  = rep(names(groups), each = 8),
Replicate  = rep(1:8, times = length(groups)),
PeptideSum = as.vector(replicate_sums_mat)
)

df$Condition <- factor(df$Condition, levels = paste0("T", 1:10))

# Mean + SD
summary_df <- df %>%
group_by(Condition) %>%
summarise(
mean = mean(PeptideSum, na.rm = TRUE),
sd   = sd(PeptideSum,   na.rm = TRUE),
.groups = "drop"
)

y_max <- max(summary_df$mean + summary_df$sd, df$PeptideSum) * 1.10

# Plot (ONLY change: fill = Condition + scale_fill_manual)
plot_peptides <- ggplot(summary_df, aes(x = Condition, y = mean, fill = Condition)) +
geom_bar(stat = "identity", width = 0.85) +
geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
width = 0.3, linewidth = 1, color = "black") +
geom_point(
data = df, aes(x = Condition, y = PeptideSum),
position = position_jitter(width = 0.10, height = 0),
size = 3.5, color = "black", fill = "white",
shape = 21, stroke = 1.2
) +
scale_fill_manual(values = group_colors) +
labs(
x = "Time points",
y = "Number of peptides",
title = "Total Peptides per Condition"
) +
scale_y_continuous(expand = c(0, 0), limits = c(0, y_max)) +
theme_classic(base_size = 20) +
theme(
legend.position = "none",
axis.line        = element_line(color = "black", linewidth = 1.2),
axis.ticks       = element_line(color = "black", linewidth = 1),
axis.text.x      = element_text(size = 22, color = "black", face = "bold"),
axis.text.y      = element_text(size = 22, color = "black", face = "bold"),
axis.title.x     = element_text(size = 24, color = "black", face = "bold", margin = margin(t = 18)),
axis.title.y     = element_text(size = 24, color = "black", face = "bold", margin = margin(r = 18)),
plot.title       = element_text(size = 26, color = "black", face = "bold",
hjust = 0.5, margin = margin(b = 22)),
plot.background  = element_rect(fill = "white", color = NA),
panel.background = element_rect(fill = "white", color = NA)
)

print(plot_peptides)

# Export PNG
ggsave(
"peptide_totals_T1_T10_colored_opencircle.png",
plot_peptides, width = 10, height = 7, dpi = 300, bg = "white"
)

################################################################################


################################################################################
################################################################################
############Protein IDs per Sample with Cumulative and Shared Trends############
################################################################################
################################################################################

library(dplyr)
library(ggplot2)

# --- NEW: group colors (non-grey palette) ---
group_colors <- c(
"T1"  = "#1b9e77",
"T2"  = "#d95f02",
"T3"  = "#7570b3",
"T4"  = "#e7298a",
"T5"  = "#66a61e",
"T6"  = "#e6ab02",
"T7"  = "#a6761d",
"T8"  = "#1f78b4",
"T9"  = "#ff7f00",
"T10" = "#6a3d9a"
)

# 1. Define sample columns in T1_1 to T10_8 order
sample_cols <- c(
paste0("T1_", 1:8), paste0("T2_", 1:8), paste0("T3_", 1:8), paste0("T4_", 1:8),
paste0("T5_", 1:8), paste0("T6_", 1:8), paste0("T7_", 1:8), paste0("T8_", 1:8),
paste0("T9_", 1:8), paste0("T10_", 1:8)
)

# 2. Make protein presence matrix (TRUE = present, FALSE/NA = missing)
protein_matrix <- !is.na(as.matrix(cor_remove_filter70[sample_cols]))

# 3. Protein count per sample (number of TRUEs per column)
protein_counts_per_sample <- colSums(protein_matrix)
sample_names <- sample_cols

# 4. Calculate cumulative and shared trends
cumulative_protein <- sapply(1:length(sample_names), function(i) {
sum(rowSums(protein_matrix[, 1:i, drop = FALSE]) > 0)
})

shared_protein <- sapply(1:length(sample_names), function(i) {
sum(rowSums(protein_matrix[, 1:i, drop = FALSE]) == i)
})

# 5. Build dataframe for plotting
df <- data.frame(
Sample = factor(sample_names, levels = sample_names),
ProteinCount = protein_counts_per_sample,
Cumulative = cumulative_protein,
Shared = shared_protein
)

# --- NEW: add group (T1..T10) from sample names, for coloring bars ---
df$Group <- sub("_(.*)$", "", as.character(df$Sample))

# 6. Set y-axis upper limit for space above trend lines/bars
y_max <- max(df$ProteinCount, df$Cumulative, df$Shared, na.rm = TRUE) * 1.10

# 7. Set plot width automatically
n_samples <- length(sample_names)
plot_width <- max(12, n_samples * 0.25)

# 8. Plot with all sample names + group-colored bars
p <- ggplot(df, aes(x = Sample)) +
geom_bar(aes(y = ProteinCount, fill = Group), stat = "identity", width = 0.85) +
geom_line(aes(y = Cumulative, group = 1, color = "Cumulative"), linewidth = 1.2) +
geom_point(aes(y = Cumulative, color = "Cumulative"), size = 2.2) +
geom_line(aes(y = Shared, group = 1, color = "Shared"), linewidth = 1.2) +
geom_point(aes(y = Shared, color = "Shared"), size = 2.2) +
scale_fill_manual(values = group_colors) +
scale_color_manual(values = c("Cumulative" = "#008000",  # vivid green
"Shared" = "black")) +  # vivid orange
labs(
y = "Number of proteins",
x = NULL,
title = NULL
) +
scale_x_discrete(labels = sample_names) +
scale_y_continuous(expand = c(0, 0), limits = c(0, y_max)) +
theme_classic(base_size = 22) +
theme(
legend.title = element_blank(),
legend.position = "none",
axis.text.x = element_text(angle = 90, vjust = 0, hjust = 0, size = 14),
axis.text.y = element_text(size = 22, color = "black", face = "bold"),
axis.title.x = element_text(size = 24, color = "black", face = "bold", margin = margin(t = 18)),
axis.title.y = element_text(size = 24, color = "black", face = "bold", margin = margin(r = 18)),
plot.title   = element_text(size = 26, color = "black", face = "bold", hjust = 0.5, margin = margin(b = 22)),
plot.background  = element_rect(fill = "white", color = NA),
panel.background = element_rect(fill = "white", color = NA)
)

print(p)

# 9. Export PNG, white background, 300 dpi, auto width, height fixed
ggsave("Cumulative and Shared Trends_number_of_protein.png",
p, width = plot_width, height = 7, dpi = 300, bg = "white")

################################################################################
################################################################################
################################################################################


################################################################################
######## Total Proteins per Condition — WITHOUT and WITH Trends #################
################################################################################

library(dplyr)
library(ggplot2)

# Group colors
group_colors <- c(
"T1"="#1b9e77","T2"="#d95f02","T3"="#7570b3","T4"="#e7298a","T5"="#66a61e",
"T6"="#e6ab02","T7"="#a6761d","T8"="#1f78b4","T9"="#ff7f00","T10"="#6a3d9a"
)

# Samples
sample_cols <- unlist(lapply(1:10, \(t) paste0("T", t, "_", 1:8)))

# Numeric safeguard
cor_remove_filter70[sample_cols] <- lapply(cor_remove_filter70[sample_cols],
\(x) suppressWarnings(as.numeric(x)))

# Presence (LFQ > 0)
protein_matrix <- sapply(sample_cols, \(c) cor_remove_filter70[[c]] > 0)
protein_matrix[is.na(protein_matrix)] <- FALSE

# Per-sample counts
df <- data.frame(
Sample       = sample_cols,
Condition    = factor(sub("_(.*)$", "", sample_cols), levels = paste0("T", 1:10)),
ProteinCount = colSums(protein_matrix)
)

# Mean + SD per condition
sum_df <- df %>%
group_by(Condition) %>%
summarise(mean = mean(ProteinCount), sd = sd(ProteinCount), .groups = "drop")

# Common theme
th <- theme_classic(base_size = 20) +
theme(
legend.position = "none",
axis.line    = element_line(color = "black", linewidth = 1.2),
axis.ticks   = element_line(color = "black", linewidth = 1),
axis.text.x  = element_text(size = 22, color = "black", face = "bold"),
axis.text.y  = element_text(size = 22, color = "black", face = "bold"),
axis.title.x = element_text(size = 24, color = "black", face = "bold", margin = margin(t = 18)),
axis.title.y = element_text(size = 24, color = "black", face = "bold", margin = margin(r = 18)),
plot.title   = element_text(size = 26, color = "black", face = "bold", hjust = 0.5, margin = margin(b = 22))
)

# Base plot function
base_plot <- function(title_text, y_max) {
ggplot(sum_df, aes(Condition, mean, fill = Condition)) +
geom_bar(stat = "identity", width = 0.85) +
geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
width = 0.3, linewidth = 1, color = "black") +
geom_point(data = df, aes(Condition, ProteinCount),
position = position_jitter(width = 0.10),
size = 3.5, shape = 21, stroke = 1.2,
color = "black", fill = "white") +
scale_fill_manual(values = group_colors) +
labs(x = NULL, y = "Number of identified proteins", title = NULL) +
scale_y_continuous(expand = c(0, 0), limits = c(0, y_max)) +
th
}

# --- Plot 1: WITHOUT trends ---
y1 <- max(sum_df$mean + sum_df$sd, df$ProteinCount) * 1.10
p1 <- base_plot("Total Proteins per Condition", y1)
print(p1)
ggsave("protein_totals_T1_T10_colored.png", p1, width = 10, height = 7, dpi = 300, bg = "white")

# --- Plot 2: WITH cumulative + shared ---
conds <- paste0("T", 1:10)
cond_cols <- lapply(conds, \(t) paste0(t, "_", 1:8))

trend <- data.frame(
Condition = factor(conds, levels = conds),
Cumulative = sapply(1:10, \(i) { cols <- unlist(cond_cols[1:i]); sum(rowSums(protein_matrix[, cols, drop=FALSE]) > 0) }),
Shared     = sapply(1:10, \(i) { cols <- unlist(cond_cols[1:i]); sum(rowSums(protein_matrix[, cols, drop=FALSE]) == length(cols)) })
)

y2 <- max(y1, trend$Cumulative, trend$Shared) * 1.10
p2 <- base_plot("Total Proteins per Condition with Cumulative and Shared Trends", y2) +
geom_line(data = trend, aes(Condition, Cumulative, group = 1),
linewidth = 1.2, color = "#008000") +
geom_point(data = trend, aes(Condition, Cumulative),
size = 2.6, color = "#008000") +
geom_line(data = trend, aes(Condition, Shared, group = 1),
linewidth = 1.2, color = "black") +
geom_point(data = trend, aes(Condition, Shared),
size = 2.6, color = "black")

print(p2)
ggsave("protein_totals_T1_T10_colored_with_trends.png",
p2, width = 10, height = 7, dpi = 300, bg = "white")

################################################################################


################################################################################
# Perseus-style imputation (Total matrix) with fixed seed for reproducibility
# Normal( μ - 1.8·σ , (0.3·σ)^2 )
################################################################################

# 1) Select LFQ columns (T#_#)
expr_cols <- grep("^T\\d+_\\d+$", colnames(cor_remove_filter70), value = TRUE)
if (!length(expr_cols)) stop("No LFQ columns (T#_#) found for imputation.")

# 2) Build numeric matrix; treat non-finite as missing
X <- as.matrix(cor_remove_filter70[expr_cols])
storage.mode(X) <- "double"
X[!is.finite(X)] <- NA_real_

# 3) Compute μ and σ from ALL valid values across selected columns (Total matrix)
vals <- X[!is.na(X)]
if (length(vals) < 2) stop("Not enough valid values to estimate μ and σ for imputation.")
mu  <- mean(vals)
sig <- sd(vals)
if (!is.finite(sig) || sig == 0) stop("σ is not finite or zero; cannot perform Perseus-style imputation.")

# 4) Define imputation Gaussian (Width = 0.3, Down shift = 1.8)
mu_imp <- mu - 1.8 * sig
sd_imp <- 0.3 * sig

# 5) Impute all missing values with fixed seed
set.seed(1)  # change to any integer for different reproducible draws
miss <- which(is.na(X))
if (length(miss)) {
X[miss] <- rnorm(length(miss), mean = mu_imp, sd = sd_imp)
}

# 6) Replace and save
imputed_data <- cor_remove_filter70
imputed_data[expr_cols] <- X

# Save TXT
write.table(imputed_data,
"P337_02B_proteinGroups_filtered_pepGT1_log2LFQ_valid70_imputed_fixedseed.txt",
sep = "\t", quote = FALSE, row.names = FALSE)

# Save Excel
library(openxlsx)
write.xlsx(imputed_data,
"P337_02B_proteinGroups_filtered_pepGT1_log2LFQ_valid70_imputed_fixedseed.xlsx",
overwrite = TRUE)

################################################################################
################################################################################
################################################################################

################################################################################
# Quantile normalization after imputation (base R) + save as Excel
################################################################################

# 1) Select LFQ columns
expr_cols <- grep("^T\\d+_\\d+$", colnames(imputed_data), value = TRUE)

# 2) Quantile normalization function (base R)
qn <- function(m) {
r <- apply(m, 2, rank, ties.method = "min")
s <- apply(m, 2, sort)
m_avg <- rowMeans(s)
m[] <- m_avg[r]
m
}

# 3) Apply quantile normalization
quantile_norm_data <- imputed_data
quantile_norm_data[expr_cols] <- qn(as.matrix(imputed_data[expr_cols]))

# 4) Save as TXT
write.table(quantile_norm_data,
"P337_02B_proteinGroups_filtered_pepGT1_log2LFQ_valid70_imputed_fixedseed_quantilenorm.txt",
sep = "\t", quote = FALSE, row.names = FALSE)

# 5) Save as Excel
library(openxlsx)
write.xlsx(quantile_norm_data,
"P337_02B_proteinGroups_filtered_pepGT1_log2LFQ_valid70_imputed_fixedseed_quantilenorm.xlsx",
overwrite = TRUE)



################################################################################
################################################################################
################################################################################

################################################################################
# 4-PANEL sensitivity plots (ALL SAMPLES on x-axis, colored by Group):
# Show ALL protein distributions (dots) + THICK median ± SD overlay per sample
# Panels: RAW -> After 70% filter -> Imputed -> Quantile Normalized
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# --- Sample columns in order (80 samples) ---
sample_cols <- c(
paste0("T1_",  1:8), paste0("T2_",  1:8), paste0("T3_",  1:8), paste0("T4_",  1:8),
paste0("T5_",  1:8), paste0("T6_",  1:8), paste0("T7_",  1:8), paste0("T8_",  1:8),
paste0("T9_",  1:8), paste0("T10_", 1:8)
)

# --- Color palette for T1–T10 ---
group_colors <- setNames(
c(
"#1b9e77",  # teal green
"#d95f02",  # strong orange
"#7570b3",  # purple-blue
"#e7298a",  # magenta
"#66a61e",  # saturated green
"#e6ab02",  # golden yellow
"#a6761d",  # brown-orange (not grey)
"#1f78b4",  # deep blue
"#ff7f00",  # bright orange
"#6a3d9a"   # deep purple
),
paste0("T", 1:10)
)

# --- Helper: long-format with all protein values (robust; forces numeric) ---
make_long_all <- function(dat, label) {
cols_present <- intersect(sample_cols, colnames(dat))
if (length(cols_present) == 0) stop(paste0("No sample columns found in: ", label))

tmp <- as.data.frame(dat[, cols_present, drop = FALSE])
tmp[] <- lapply(tmp, function(x) suppressWarnings(as.numeric(x)))

out <- pivot_longer(
data = cbind(.row = seq_len(nrow(tmp)), tmp),
cols = - .row,
names_to = "Sample",
values_to = "LFQ"
) %>%
mutate(
Group   = factor(sub("_.*", "", Sample), levels = paste0("T", 1:10)),
Sample  = factor(Sample, levels = sample_cols),
Dataset = label
)
out
}

# --- Helper: compute per-sample MEDIAN and SD ---
make_median_sd <- function(df_long) {
df_long %>%
group_by(Sample, Group, Dataset) %>%
summarise(
median = median(LFQ, na.rm = TRUE),
sd     = sd(LFQ, na.rm = TRUE),
.groups = "drop"
)
}

# --- Plot function: dots + THICK median±SD overlay ---
plot_dist_medianSD <- function(df_long, df_sum, title_text, ylims) {
ggplot() +
# All protein values (background cloud)
geom_point(
data = df_long,
aes(x = Sample, y = LFQ, color = Group),
position = position_jitter(width = 0.20, height = 0),
size = 0.22, alpha = 0.05
) +
# SD error bars around median (THICKER)
geom_errorbar(
data = df_sum,
aes(x = Sample, ymin = median - sd, ymax = median + sd, color = Group),
width = 1,
linewidth = 2,
alpha = 1
) +
# Median points (BIGGER)
geom_point(
data = df_sum,
aes(x = Sample, y = median, color = Group),
size = 2
) +
scale_color_manual(values = group_colors) +
labs(
title = title_text,
x = NULL,
y = "log2 LFQ intensity"
) +
theme_classic(base_size = 20) +
theme(
axis.text.x = element_text(size = 20, angle = 90, hjust = 1, vjust = 1),
axis.text.y = element_text(size = 12, color = "black"),
axis.title.x = element_text(face = "plain"),
axis.title.y = element_text(face = "plain"),
plot.title  = element_text(face = "plain", hjust = 0.5, size = 32),
legend.position = "none"
) +
coord_cartesian(ylim = ylims)
}

# --- Build long data for your 4 stages (objects must exist from your pipeline) ---
L_raw <- make_long_all(log2_data_filtered2, "RAW (log2 LFQ)")
L_70  <- make_long_all(cor_remove_filter70, "After 70% filter")
L_imp <- make_long_all(imputed_data,        "After imputation")
L_qn  <- make_long_all(quantile_norm_data,  "After quantile normalization")

# OPTIONAL: If plotting becomes too slow/dense, uncomment this to sample dots per sample:
# set.seed(1)
# L_raw <- L_raw %>% group_by(Sample) %>% slice_sample(n = min(3000, n())) %>% ungroup()
# L_70  <- L_70  %>% group_by(Sample) %>% slice_sample(n = min(3000, n())) %>% ungroup()
# L_imp <- L_imp %>% group_by(Sample) %>% slice_sample(n = min(3000, n())) %>% ungroup()
# L_qn  <- L_qn  %>% group_by(Sample) %>% slice_sample(n = min(3000, n())) %>% ungroup()

# --- Summaries (MEDIAN ± SD) ---
S_raw <- make_median_sd(L_raw)
S_70  <- make_median_sd(L_70)
S_imp <- make_median_sd(L_imp)
S_qn  <- make_median_sd(L_qn)

# --- Common y limits across all panels (proteins + median±SD) ---
all_vals <- c(
L_raw$LFQ, L_70$LFQ, L_imp$LFQ, L_qn$LFQ,
S_raw$median - S_raw$sd, S_raw$median + S_raw$sd,
S_70$median  - S_70$sd,  S_70$median  + S_70$sd,
S_imp$median - S_imp$sd, S_imp$median + S_imp$sd,
S_qn$median  - S_qn$sd,  S_qn$median  + S_qn$sd
)
ylims <- range(all_vals, na.rm = TRUE)

# --- Build 4 plots ---
p1 <- plot_dist_medianSD(L_raw, S_raw, "all raw proteins", ylims)
p2 <- plot_dist_medianSD(L_70,  S_70,  "70% filter", ylims)
p3 <- plot_dist_medianSD(L_imp, S_imp, "Imputation", ylims)
p4 <- plot_dist_medianSD(L_qn,  S_qn,  "Quantile normalization", ylims)

# --- Combine into one 2x2 figure ---
p_all <- (p1 + p2) / (p3 + p4)
print(p_all)

# --- Save combined + individual ---
ggsave(
"Fig2_sensitivity_4panel_allProteins_medianSD_THICK.png",
p_all, width = 22, height = 12, dpi = 300, bg = "white"
)

ggsave("Fig2A_RAW_allProteins_medianSD_THICK.png",     p1, width = 20, height = 8, dpi = 300, bg = "white")
ggsave("Fig2B_valid70_allProteins_medianSD_THICK.png", p2, width = 20, height = 8, dpi = 300, bg = "white")
ggsave("Fig2C_imputed_allProteins_medianSD_THICK.png", p3, width = 20, height = 8, dpi = 300, bg = "white")
ggsave("Fig2D_QN_allProteins_medianSD_THICK.png",      p4, width = 20, height = 8, dpi = 300, bg = "white")


################################################################################
################################################################################
################################################################################

# ---- [Load Libraries] ----
library(mixOmics)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggrepel)

# ---- [Create Output Directory] ----
output_dir <- "sPLSDA_2D_Plots_T1_T10"
dir.create(output_dir, showWarnings = FALSE)

# ---- [Prepare Data: quantile_norm_data, columns T1_1 ... T10_8] ----
sample_cols <- c(
paste0("T1_", 1:8), paste0("T2_", 1:8), paste0("T3_", 1:8), paste0("T4_", 1:8),
paste0("T5_", 1:8), paste0("T6_", 1:8), paste0("T7_", 1:8), paste0("T8_", 1:8),
paste0("T9_", 1:8), paste0("T10_", 1:8)
)
X <- as.data.frame(t(quantile_norm_data[, sample_cols]))  # samples x proteins
sample_names <- rownames(X)
Y <- factor(sub("_.*", "", sample_names), levels = paste0("T", 1:10))

# ---- [Color palette for T1–T10] ----
group_colors <- setNames(
c(
"#1b9e77",  # teal green
"#d95f02",  # strong orange
"#7570b3",  # purple-blue
"#e7298a",  # magenta
"#66a61e",  # saturated green
"#e6ab02",  # golden yellow
"#a6761d",  # brown-orange (not grey)
"#1f78b4",  # deep blue
"#ff7f00",  # bright orange
"#6a3d9a"   # deep purple
),
paste0("T", 1:10)
)


custom_theme_splsda <- function(base_size = 14) {
theme_minimal(base_size = base_size) +
theme(
panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
plot.title = element_text(size = base_size + 4, face = "bold", hjust = 0.5),
axis.title = element_text(size = base_size + 2, face = "bold"),
axis.text = element_text(size = base_size),
legend.title = element_text(size = base_size, face = "bold"),
legend.text = element_text(size = base_size),
legend.position = "bottom"
)
}

# ---- [Tuning sPLS-DA] ----
set.seed(1)
optimal_ncomp <- 9
list_keepX <- c(25, 50, 100)

cat("🔍 Tuning sPLS-DA for T1-T10 groups...\n")
tune_result <- tune.splsda(X, Y, ncomp = optimal_ncomp,
validation = "Mfold", folds = 5,
dist = "centroids.dist", measure = "BER",
test.keepX = list_keepX, nrepeat = 10, progressBar = TRUE)
optimal_keepX <- tune_result$choice.keepX[1:optimal_ncomp]
print(optimal_keepX)

# ---- [keepX Barplot] ----
keepx_df <- data.frame(
Component = factor(paste0("Comp", seq_along(optimal_keepX))),
keepX = optimal_keepX
)
gg_keepx <- ggplot(keepx_df, aes(x = Component, y = keepX)) +
geom_bar(stat = "identity", fill = "#1f78b4", alpha = 0.85, width = 0.6) +
geom_text(aes(label = keepX), vjust = -0.5, size = 5) +
labs(title = "Optimal keepX per Component", x = "Component", y = "Variables Selected") +
theme_minimal(base_size = 14) +
theme(panel.border = element_rect(color = "black", fill = NA))
print(gg_keepx)
ggsave(file.path(output_dir, "keepX_per_component_T1_T10.png"),
plot = gg_keepx, dpi = 300, width = 6, height = 5, bg = "white")

# ---- [Fit Final Model] ----
splsda_model <- splsda(X, Y, ncomp = optimal_ncomp, keepX = optimal_keepX)

# ---- [Explained Variance] ----
expl_var <- apply(splsda_model$variates$X^2, 2, sum) / sum(splsda_model$X^2)

# ---- [Manual Ellipse Function] ----
desired_ellipse_level <- 0.95
compute_ellipse <- function(mean, cov, level = desired_ellipse_level, npoints = 100) {
angles <- seq(0, 2 * pi, length.out = npoints)
radius <- sqrt(qchisq(level, df = 2))
eig <- eigen(cov)
axes <- radius * t(eig$vectors %*% diag(sqrt(eig$values)))
ellipse <- t(axes %*% rbind(cos(angles), sin(angles))) + matrix(rep(mean, each = npoints), ncol = 2, byrow = FALSE)
df <- as.data.frame(ellipse)
colnames(df) <- c("comp1", "comp2")
return(df)
}

# ---- [Plot All Component Pairs] ----
for (i in 1:(optimal_ncomp - 1)) {
for (j in (i + 1):optimal_ncomp) {

plot_data <- data.frame(
comp1 = splsda_model$variates$X[, i],
comp2 = splsda_model$variates$X[, j],
Group = Y,
Sample = rownames(X)
)

group_centroids <- plot_data %>%
group_by(Group) %>%
summarise(comp1 = mean(comp1), comp2 = mean(comp2), count = n(), .groups = "drop")

ellipse_data <- plot_data %>%
group_by(Group) %>%
do({
group_data <- .[, c("comp1", "comp2")]   # <-- Fixed: base R select
ell <- compute_ellipse(colMeans(group_data), cov(group_data))
ell$Group <- unique(.$Group)
ell
}) %>% ungroup()

x_lab <- paste0("Component ", i, " (", round(expl_var[i] * 100, 1), "%)")
y_lab <- paste0("Component ", j, " (", round(expl_var[j] * 100, 1), "%)")

p <- ggplot(plot_data, aes(x = comp1, y = comp2, color = Group, fill = Group)) +
geom_point(size = 4, alpha = 0.9) +
geom_polygon(data = ellipse_data, aes(group = Group), alpha = 0.18, color = NA) +
geom_path(data = ellipse_data, aes(group = Group), linewidth = 1) +
geom_text_repel(data = group_centroids,
aes(label = paste0(Group, "\n(n=", count, ")")),
color = "black", size = 5, fontface = "bold",
max.overlaps = 100, box.padding = 0.6, point.padding = 0.6) +
scale_color_manual(values = group_colors) +
scale_fill_manual(values = group_colors) +
labs(
title = paste0("sPLS-DA: (Comp ", i, " vs ", j, ")"),
x = x_lab, y = y_lab
) +
custom_theme_splsda()

print(p)
ggsave(file.path(output_dir, paste0("sPLS-DA_T1-T10_Comp", i, "_vs_Comp", j, ".png")),
plot = p, dpi = 300, width = 10, height = 10, bg = "white")
}
}

cat("\n✅ All sPLS-DA 2D plots saved in folder:", output_dir, "\n")

# ---- [Summary Report] ----
cat("\n📊 sPLS-DA Summary Report\n")
cat("──────────────────────────────\n")
cat("✔ Number of Components Used:", optimal_ncomp, "\n")
cat("✔ keepX per Component:", paste(optimal_keepX, collapse = ", "), "\n")
cat("✔ Confidence Level for Ellipses:", desired_ellipse_level * 100, "%\n")
n_to_report <- min(2, length(expl_var))
cat("✔ Explained Variance (First", n_to_report, "components):",
paste0(round(expl_var[1:n_to_report] * 100, 1), collapse = "%, "), "%\n")
cat("✔ Output folder:", output_dir, "\n")

################################################################################
################################################################################
################################################################################

################################################################################
#                                PCA ANALYSIS                                   #
#        mixOmics-based PCA with Elbow Plot + Colored Ellipses + Borders        #
#                     All 2D Component Comparisons (T1–T10)                     #
################################################################################

# ---- [Load Libraries] ----
library(mixOmics)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(tidyr)

# ---- [Create Output Directory] ----
pca_output_dir <- "PCA_2D_Plots_T1_T10"
dir.create(pca_output_dir, showWarnings = FALSE)

# ---- [Prepare Data] ----
sample_cols <- c(
paste0("T1_", 1:8), paste0("T2_", 1:8), paste0("T3_", 1:8), paste0("T4_", 1:8),
paste0("T5_", 1:8), paste0("T6_", 1:8), paste0("T7_", 1:8), paste0("T8_", 1:8),
paste0("T9_", 1:8), paste0("T10_", 1:8)
)

X <- as.data.frame(t(quantile_norm_data[, sample_cols]))  # samples × proteins
Y <- factor(sub("_.*", "", rownames(X)), levels = paste0("T", 1:10))

# ---- [Color Palette for T1–T10] ----
group_colors <- c(
"T1"  = "#1b9e77",  # teal green
"T2"  = "#d95f02",  # strong orange
"T3"  = "#7570b3",  # purple-blue
"T4"  = "#e7298a",  # magenta
"T5"  = "#66a61e",  # saturated green
"T6"  = "#e6ab02",  # golden yellow
"T7"  = "#a6761d",  # brown-orange
"T8"  = "#1f78b4",  # deep blue
"T9"  = "#ff7f00",  # bright orange
"T10" = "#6a3d9a"   # deep purple
)

# ---- [Run PCA using mixOmics] ----
set.seed(1)
pca_res <- pca(X, ncomp = min(ncol(X), 15), center = TRUE, scale = TRUE)

# ---- [Extract Variance Explained] ----
var_explained <- pca_res$prop_expl_var$X * 100  # %

# ---- [Elbow Plot] ----
elbow_df <- data.frame(
PC = 1:length(var_explained),
Variance = var_explained,
Cumulative = cumsum(var_explained)
)

gg_elbow <- ggplot(elbow_df, aes(x = PC, y = Variance)) +
geom_line(color = "#1f78b4", linewidth = 1.2) +
geom_point(size = 3, color = "#1f78b4") +
geom_text(aes(label = sprintf("%.1f%%", Variance)), vjust = -0.6, size = 4) +
geom_line(aes(y = Cumulative / max(Cumulative) * max(Variance)),
color = "gray40", linetype = "dashed", linewidth = 0.9) +
scale_x_continuous(breaks = 1:length(var_explained)) +
labs(
title = "PCA Elbow Plot (mixOmics)",
x = "Principal Component",
y = "Variance Explained (%)"
) +
theme_minimal(base_size = 16) +
theme(
panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
plot.title = element_text(face = "bold", hjust = 0.5)
)

print(gg_elbow)
ggsave(file.path(pca_output_dir, "PCA_Elbow_Plot_T1_T10_mixOmics.png"),
gg_elbow, width = 8, height = 6, dpi = 300, bg = "white")

# ---- [Determine Optimal Components] ----
optimal_ncomp <- min(which(elbow_df$Cumulative >= 90))
if (is.infinite(optimal_ncomp)) optimal_ncomp <- 5
cat("✅ Optimal number of components based on cumulative variance:", optimal_ncomp, "\n")

# ---- [Extract Scores for Selected PCs] ----
pca_scores <- as.data.frame(pca_res$variates$X[, 1:optimal_ncomp])
pca_scores$Group <- Y
pca_scores$Sample <- rownames(X)

# ---- [Function: PCA Plot with Matching Fill + Border Colors + Colored Labels] ----
plot_pca_2d <- function(df, pcx, pcy, var_explained, group_colors, output_dir) {

pcx_name <- paste0("PC", pcx)
pcy_name <- paste0("PC", pcy)

centroids <- df %>%
group_by(Group) %>%
summarise(
!!pcx_name := mean(.data[[pcx_name]]),
!!pcy_name := mean(.data[[pcy_name]]),
count = n(),
.groups = "drop"
)

p <- ggplot(df, aes_string(x = pcx_name, y = pcy_name)) +
geom_point(
aes(color = Group), 
size = 3, 
alpha = 0.5,
shape = 16) +

stat_ellipse(
aes(color = Group, fill = Group, group = Group),
type = "norm", 
level = 0.95, 
geom = "polygon",
alpha = 0.15, 
linewidth = 1.2,
color = NA) +

# ✅ Labels now match group colors
geom_text_repel(
data = centroids,
aes_string(
x = pcx_name,
y = pcy_name,
label = "paste0(Group, '\\n(n=', count, ')')",
color = "Group"
),
size = 8,
fontface = "bold",
alpha = 1,   # 👈 group text label transparency
show.legend = FALSE,
max.overlaps = 100,
box.padding = 0.6,
point.padding = 0.6
) +

scale_color_manual(values = group_colors) +
scale_fill_manual(values = group_colors) +

labs(
title = NULL,
x = paste0(pcx_name, " (", round(var_explained[pcx], 1), "%)"),
y = paste0(pcy_name, " (", round(var_explained[pcy], 1), "%)")
) +

theme_minimal(base_size = 32) +
theme(
panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
plot.title = element_text(face = "bold", hjust = 0.5),
axis.title = element_text(face = "bold"),
legend.position = "none"
)

out_name <- paste0("PCA_T1-T10_", pcx_name, "_vs_", pcy_name, ".png")
ggsave(file.path(output_dir, out_name), p, width = 8, height = 8, dpi = 300, bg = "white")

return(p)
}

# ---- [Generate All PCA 2D Plots: All Component Comparisons] ----
for (i in 1:(optimal_ncomp - 1)) {
for (j in (i + 1):optimal_ncomp) {
plot_pca_2d(pca_scores, i, j, var_explained, group_colors, pca_output_dir)
}
}

cat("\n✅ All PCA 2D plots (mixOmics, All PCs, Filled + Colored Borders) saved in folder:", pca_output_dir, "\n")

# ---- [Summary Report] ----
cat("\n📊 PCA Summary Report (mixOmics)\n")
cat("──────────────────────────────────────\n")
cat("✔ Total Components Computed:", length(var_explained), "\n")
cat("✔ Optimal Components (≥90% cumulative variance):", optimal_ncomp, "\n")
cat("✔ Variance Explained (First 5 PCs):", paste0(round(var_explained[1:5], 1), collapse = "%, "), "%\n")
cat("✔ Output Folder:", pca_output_dir, "\n")

################################################################################
# End of PCA Section
################################################################################




################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

###############  Volcano plot with enrichment analysis  ########################

################################################################################
################################################################################
################################################################################

################################################################################
# Volcano plots (T1–T10 pairwise) using limma moderated t-statistic + BH FDR
# Outputs per comparison:
#   1) All_proteins_<Tearly_vs_Tlate>.txt
#   2) <Tlate>_up_<Tearly>_down_<comp>.txt
#   3) <Tearly>_up_<Tlate>_down_<comp>.txt
#   4) Volcano_<comp>.png
################################################################################

suppressPackageStartupMessages({
library(limma)
library(ggplot2)
library(dplyr)
})

# ---------------------------
# 0) Samples + colors
# ---------------------------
sample_cols <- c(
paste0("T1_", 1:8), paste0("T2_", 1:8), paste0("T3_", 1:8), paste0("T4_", 1:8),
paste0("T5_", 1:8), paste0("T6_", 1:8), paste0("T7_", 1:8), paste0("T8_", 1:8),
paste0("T9_", 1:8), paste0("T10_", 1:8)
)

group_colors <- setNames(
c(
"#1b9e77",  # T1
"#d95f02",  # T2
"#7570b3",  # T3
"#e7298a",  # T4
"#66a61e",  # T5
"#e6ab02",  # T6
"#a6761d",  # T7
"#1f78b4",  # T8
"#ff7f00",  # T9
"#6a3d9a"   # T10
),
paste0("T", 1:10)
)

# ---------------------------
# 1) Choose expression source
# ---------------------------
expr_source <- quantile_norm_data   # <-- or imputed_data (often better for testing)

missing_cols <- setdiff(sample_cols, colnames(expr_source))
if (length(missing_cols) > 0) {
stop("These sample columns are missing from expr_source: ",
paste(missing_cols, collapse = ", "))
}

# Build numeric matrix (proteins x samples)
E <- as.matrix(expr_source[, sample_cols, drop = FALSE])
storage.mode(E) <- "double"

# ---------------------------
# 2) Row names (Gene_name fallback Protein IDs)
# ---------------------------
if (all(c("Gene_name", "Protein IDs") %in% colnames(expr_source))) {

gene_ids_raw <- ifelse(
is.na(expr_source$Gene_name) | expr_source$Gene_name == "",
expr_source$`Protein IDs`,
expr_source$Gene_name
)

gene_ids <- vapply(gene_ids_raw, function(g) {
parts <- unlist(strsplit(as.character(g), "[ ;]+"))
parts <- unique(trimws(parts))
parts <- parts[nzchar(parts)]
paste(parts, collapse = ";")
}, character(1))

rownames(E) <- make.unique(gene_ids)

} else if (!is.null(rownames(expr_source)) && all(nzchar(rownames(expr_source)))) {

rownames(E) <- make.unique(rownames(expr_source))

} else {
rownames(E) <- make.unique(as.character(seq_len(nrow(E))))
}

# ---------------------------
# 3) Volcano parameters + output dir
# ---------------------------
output_dir <- getwd()
volcano_output_dir <- file.path(output_dir, "Volcano_Results_limma_BH")
dir.create(volcano_output_dir, showWarnings = FALSE)

pval_threshold <- 0.05
fc_threshold   <- log2(2)

cat("📊 Starting limma volcano analysis (moderated t + BH FDR) using expr_source:",
deparse(substitute(expr_source)), "\n")

# ---------------------------
# 4) Group -> columns mapping
# ---------------------------
group_names <- paste0("T", 1:10)
group_columns <- lapply(group_names, function(g) grep(paste0("^", g, "_"), colnames(E), value = TRUE))
names(group_columns) <- group_names

# ---------------------------
# 5) Function: run one comparison
# ---------------------------
run_one_comparison <- function(group_early, group_late) {

early_cols <- group_columns[[group_early]]
late_cols  <- group_columns[[group_late]]

if (!length(early_cols) || !length(late_cols)) return(invisible(NULL))

# Subset matrix to these samples
Esub <- E[, c(early_cols, late_cols), drop = FALSE]

# Design matrix for 2-group comparison
grp <- factor(
c(rep(group_early, length(early_cols)), rep(group_late, length(late_cols))),
levels = c(group_early, group_late)
)
design <- model.matrix(~0 + grp)
colnames(design) <- levels(grp)

# limma fit + moderated t
fit <- lmFit(Esub, design)
contr <- makeContrasts(contrasts = paste0(group_late, "-", group_early), levels = design)
fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2,  robust = TRUE)

# Extract results for all proteins
tab <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "none")

df <- data.frame(
Gene      = rownames(tab),
log2FC    = tab$logFC,
t_mod     = tab$t,
p_value   = tab$P.Value,
padj_BH   = tab$adj.P.Val,
AveExpr   = tab$AveExpr,
B         = tab$B,
negLog10P = -log10(pmax(tab$adj.P.Val, .Machine$double.xmin)),
stringsAsFactors = FALSE
)

# Labels for coloring
label_up   <- paste0(group_late, "_up_", group_early, "_down")   # higher in late
label_down <- paste0(group_early, "_up_", group_late, "_down")   # higher in early

df$group <- "Non-significant"
df$group[df$log2FC >  fc_threshold & df$padj_BH < pval_threshold] <- label_up
df$group[df$log2FC < -fc_threshold & df$padj_BH < pval_threshold] <- label_down

# Output folder for this comparison
comp <- paste0(group_early, "_vs_", group_late)
comp_dir <- file.path(volcano_output_dir, comp)
dir.create(comp_dir, showWarnings = FALSE)

# Save tables
write.table(df,
file = file.path(comp_dir, paste0("All_proteins_", comp, ".txt")),
sep = "\t", row.names = FALSE, quote = FALSE)

write.table(df[df$group == label_up, ],
file = file.path(comp_dir, paste0(group_late, "_up_", group_early, "_down_", comp, ".txt")),
sep = "\t", row.names = FALSE, quote = FALSE)

write.table(df[df$group == label_down, ],
file = file.path(comp_dir, paste0(group_early, "_up_", group_late, "_down_", comp, ".txt")),
sep = "\t", row.names = FALSE, quote = FALSE)

# Volcano plot
color_map <- c(
setNames(group_colors[group_late],  label_up),
setNames(group_colors[group_early], label_down),
"Non-significant" = "grey80"
)

p <- ggplot(df, aes(x = log2FC, y = negLog10P)) +
geom_point(aes(color = group), alpha = 0.65, size = 5) +
scale_color_manual(values = color_map, breaks = c(label_down, label_up, "Non-significant")) +
geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", linewidth = 1) +
geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", linewidth = 1) +
labs(
title = NULL,
x = "Log2 Fold Change",
y = "-Log10 BH-adjusted p-value",
color = NULL
) +
theme_classic(base_size = 30) +
theme(
plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
axis.title = element_text(size = 30, face = "bold"),
axis.text  = element_text(size = 30),
legend.text = element_text(size = 25),
legend.position = "bottom"
)

ggsave(file.path(comp_dir, paste0("Volcano_", comp, ".png")),
p, dpi = 300, width = 12, height = 12, bg = "white")

cat("✅ Volcano saved for", comp, "\n")
invisible(df)
}

# ---------------------------
# 6) Run all pairwise comparisons
# ---------------------------
for (i in 1:(length(group_names) - 1)) {
for (j in (i + 1):length(group_names)) {
run_one_comparison(group_names[i], group_names[j])
}
}

cat("🎉 limma volcano analysis complete. Results in:", normalizePath(volcano_output_dir), "\n")


################################################################################
################################################################################
################################################################################

################################################################################
# Volcano annotation script (Fully namespaced with :: to avoid select conflicts)
# - Reads:   All_proteins_<comp>.txt
# - Annotates ONLY: T#_up_T#_down_<comp>.txt files (raw up/down)
# - Writes:
#   1) <set_label>_annotated.txt
#   2) Merged_matrix_<set_label>.txt   (Whole volcano + annotation columns only)
################################################################################

suppressPackageStartupMessages({
library(org.Bt.eg.db)
library(AnnotationDbi)
library(dplyr)
library(tidyr)
library(ReactomePA)
library(biomaRt)
library(KEGGREST)
library(GO.db)
library(tools)
})

# --- Safety checks ---
if (!exists("volcano_output_dir")) stop("❌ 'volcano_output_dir' not found. Run volcano script first.")
if (!dir.exists(volcano_output_dir)) stop("❌ Folder does not exist: ", volcano_output_dir)

# --- Helper: split multi-gene entries ("A;B;C") and keep mapping back ---
split_genes <- function(gene_vec) {
gene_vec <- as.character(gene_vec)
gene_vec[is.na(gene_vec)] <- ""
gene_list <- strsplit(gene_vec, ";", fixed = TRUE)

data.frame(
Original    = rep(gene_vec, times = lengths(gene_list)),
Gene_Single = trimws(unlist(gene_list)),
stringsAsFactors = FALSE
) %>%
dplyr::filter(nzchar(Gene_Single))
}

# --- Cache KEGG pathway names once (optional but nice) ---
kegg_descriptions <- tryCatch({
desc <- utils::stack(KEGGREST::keggList("pathway", "bta"))
data.frame(
KEGG_Pathway      = sub("path:", "", desc$ind),
KEGG_Description  = desc$values,
stringsAsFactors  = FALSE
)
}, error = function(e) NULL)

# --- Cache KEGG ENTREZ -> pathway links once ---
kegg_link_all <- tryCatch(KEGGREST::keggLink("pathway", "bta"), error = function(e) NULL)

# --- Connect biomaRt ONCE (outside loops) ---
cow_mart <- tryCatch(
biomaRt::useMart("ensembl", dataset = "btaurus_gene_ensembl", host = "https://dec2021.archive.ensembl.org"),
error = function(e) NULL
)
human_mart <- tryCatch(
biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org"),
error = function(e) NULL
)

if (is.null(cow_mart) || is.null(human_mart)) {
message("⚠️ biomaRt connection failed (archive host). Reactome annotation will be skipped if marts are NULL.")
}

# --- List comparison folders ---
comparison_folders <- list.dirs(volcano_output_dir, recursive = FALSE, full.names = TRUE)
if (!length(comparison_folders)) stop("❌ No comparison folders found in: ", volcano_output_dir)

for (comp_folder in comparison_folders) {
comp_name <- basename(comp_folder)

# --- Read full volcano table (NEW naming) ---
all_file <- file.path(comp_folder, paste0("All_proteins_", comp_name, ".txt"))
if (!file.exists(all_file)) {
cat("❌ Missing All_proteins file for", comp_name, ":", all_file, "\n")
next
}

Whole_data <- utils::read.table(
all_file, header = TRUE, sep = "\t",
stringsAsFactors = FALSE, check.names = FALSE
)

# --- Find ONLY the raw up/down files produced by your volcano script ---
updown_files <- list.files(
comp_folder,
pattern = "^T\\d+_up_T\\d+_down_.*\\.txt$",
full.names = TRUE
)
updown_files <- updown_files[!grepl("_annotated\\.txt$", updown_files)]
updown_files <- updown_files[!grepl("^Merged_matrix_", basename(updown_files))]

if (!length(updown_files)) {
cat("❌ No raw up/down files found for", comp_name, "\n")
next
}

for (gene_file in updown_files) {
set_label <- tools::file_path_sans_ext(basename(gene_file))
cat("🔍 Annotating set:", set_label, "in", comp_name, "...\n")

gene_df <- utils::read.table(
gene_file, header = TRUE, sep = "\t", quote = "", fill = TRUE,
comment.char = "", stringsAsFactors = FALSE, check.names = FALSE
)

if (!"Gene" %in% colnames(gene_df)) {
cat("❌ 'Gene' column not found in", gene_file, "\n")
next
}

# --- Split multi-gene entries ---
gene_map <- split_genes(gene_df$Gene)
if (!nrow(gene_map)) {
cat("⚠️ No valid genes found in", gene_file, "\n")
next
}

# ==========================================================================
# 1) Cow SYMBOL -> Cow ENTREZ
# ==========================================================================
gene_entrez <- AnnotationDbi::select(
org.Bt.eg.db,
keys    = unique(gene_map$Gene_Single),
columns = c("ENTREZID", "SYMBOL"),
keytype = "SYMBOL"
) %>% as.data.frame()

gene_entrez <- gene_entrez %>%
dplyr::filter(!is.na(ENTREZID), !is.na(SYMBOL), nzchar(ENTREZID), nzchar(SYMBOL)) %>%
dplyr::distinct()

cow_entrez <- unique(gene_entrez$ENTREZID)

# ==========================================================================
# 2) GO annotation (cow ENTREZ)
# ==========================================================================
go_agg <- tryCatch({
if (!length(cow_entrez)) return(data.frame(SYMBOL = character(), GO_Pathway = character()))

go_anno <- AnnotationDbi::select(
org.Bt.eg.db,
keys    = cow_entrez,
columns = c("GO", "ONTOLOGY"),
keytype = "ENTREZID"
) %>% as.data.frame()

go_anno <- go_anno %>%
dplyr::filter(!is.na(GO), nzchar(GO))

go_terms <- AnnotationDbi::select(
GO.db,
keys    = unique(go_anno$GO),
columns = "TERM",
keytype = "GOID"
) %>% as.data.frame()

go_merged <- go_anno %>%
dplyr::left_join(go_terms, by = c("GO" = "GOID")) %>%
dplyr::left_join(gene_entrez, by = "ENTREZID") %>%
dplyr::mutate(GO_Pathway = paste(GO, ONTOLOGY, TERM, sep = " | ")) %>%
dplyr::select(SYMBOL, GO_Pathway) %>%
dplyr::filter(!is.na(SYMBOL), nzchar(SYMBOL))

go_merged %>%
dplyr::group_by(SYMBOL) %>%
dplyr::summarise(
GO_Pathway = paste(unique(stats::na.omit(GO_Pathway)), collapse = "; "),
.groups = "drop"
)
}, error = function(e) data.frame(SYMBOL = character(), GO_Pathway = character()))

# ==========================================================================
# 3) Reactome annotation via HUMAN ortholog ENTREZ (stable)
#    - getLDS Cow SYMBOL -> Human ENTREZ
#    - enrichPathway(readable=FALSE) returns ENTREZ IDs in geneID
# ==========================================================================
reactome_agg <- tryCatch({
if (is.null(cow_mart) || is.null(human_mart)) {
return(data.frame(SYMBOL = character(), Reactome_Pathway = character()))
}

cow_syms <- unique(gene_map$Gene_Single)
cow_syms <- cow_syms[nzchar(cow_syms)]
if (!length(cow_syms)) return(data.frame(SYMBOL = character(), Reactome_Pathway = character()))

chunks <- split(cow_syms, ceiling(seq_along(cow_syms) / 200))

orthologs <- do.call(rbind, lapply(chunks, function(chunk) {
tryCatch(
biomaRt::getLDS(
attributes  = "external_gene_name",
filters     = "external_gene_name",
values      = chunk,
mart        = cow_mart,
attributesL = c("entrezgene_id", "external_gene_name"),
martL       = human_mart
),
error = function(e) NULL
)
}))

if (is.null(orthologs) || !nrow(orthologs)) {
return(data.frame(SYMBOL = character(), Reactome_Pathway = character()))
}

colnames(orthologs) <- c("Cow_Gene", "Human_ENTREZID", "Human_Gene")
orthologs <- orthologs %>%
dplyr::filter(
!is.na(Cow_Gene), nzchar(Cow_Gene),
!is.na(Human_ENTREZID), nzchar(Human_ENTREZID)
) %>%
dplyr::distinct()

human_entrez <- unique(orthologs$Human_ENTREZID)
if (!length(human_entrez)) return(data.frame(SYMBOL = character(), Reactome_Pathway = character()))

enr <- ReactomePA::enrichPathway(
gene         = human_entrez,
organism     = "human",
readable     = FALSE,  # keep ENTREZ IDs in geneID
pvalueCutoff = 1,
qvalueCutoff = 1
)

df_enr <- as.data.frame(enr)
if (!nrow(df_enr)) return(data.frame(SYMBOL = character(), Reactome_Pathway = character()))

react <- df_enr %>%
dplyr::select(Reactome_Pathway = Description, geneID) %>%
tidyr::separate_rows(geneID, sep = "/") %>%
dplyr::rename(Human_ENTREZID = geneID) %>%
dplyr::left_join(orthologs, by = "Human_ENTREZID") %>%
dplyr::transmute(SYMBOL = Cow_Gene, Reactome_Pathway) %>%
dplyr::filter(!is.na(SYMBOL), nzchar(SYMBOL))

react %>%
dplyr::group_by(SYMBOL) %>%
dplyr::summarise(
Reactome_Pathway = paste(unique(stats::na.omit(Reactome_Pathway)), collapse = "; "),
.groups = "drop"
)
}, error = function(e) data.frame(SYMBOL = character(), Reactome_Pathway = character()))

# ==========================================================================
# 4) KEGG annotation (cow ENTREZ)
# ==========================================================================
kegg_agg <- tryCatch({
if (!length(cow_entrez) || is.null(kegg_link_all) || !length(kegg_link_all)) {
return(data.frame(SYMBOL = character(), kegg_results = character()))
}

formatted <- paste0("bta:", cow_entrez)
kmap <- kegg_link_all[names(kegg_link_all) %in% formatted]
if (!length(kmap)) return(data.frame(SYMBOL = character(), kegg_results = character()))

kegg_df <- data.frame(
ENTREZID     = sub("bta:", "", names(kmap)),
KEGG_Pathway = sub("path:", "", unname(kmap)),
stringsAsFactors = FALSE
)

kegg_df <- kegg_df %>%
dplyr::left_join(gene_entrez, by = c("ENTREZID" = "ENTREZID")) %>%
dplyr::filter(!is.na(SYMBOL), nzchar(SYMBOL))

if (!is.null(kegg_descriptions) && nrow(kegg_descriptions)) {
kegg_df <- kegg_df %>%
dplyr::left_join(kegg_descriptions, by = "KEGG_Pathway") %>%
dplyr::mutate(KEGG_Annotation = paste(KEGG_Description, KEGG_Pathway, sep = " | "))
} else {
kegg_df$KEGG_Annotation <- kegg_df$KEGG_Pathway
}

kegg_df %>%
dplyr::group_by(SYMBOL) %>%
dplyr::summarise(
kegg_results = paste(unique(stats::na.omit(KEGG_Annotation)), collapse = "; "),
.groups = "drop"
)
}, error = function(e) data.frame(SYMBOL = character(), kegg_results = character()))

# ==========================================================================
# 5) Merge annotations back to ORIGINAL multi-gene strings
# ==========================================================================
combined_annots <- gene_map %>%
dplyr::left_join(go_agg,       by = c("Gene_Single" = "SYMBOL")) %>%
dplyr::left_join(reactome_agg, by = c("Gene_Single" = "SYMBOL")) %>%
dplyr::left_join(kegg_agg,     by = c("Gene_Single" = "SYMBOL")) %>%
dplyr::group_by(Original) %>%
dplyr::summarise(
GO_Pathway       = ifelse(all(is.na(GO_Pathway)), NA, paste(unique(stats::na.omit(GO_Pathway)), collapse = "; ")),
Reactome_Pathway = ifelse(all(is.na(Reactome_Pathway)), NA, paste(unique(stats::na.omit(Reactome_Pathway)), collapse = "; ")),
kegg_results     = ifelse(all(is.na(kegg_results)), NA, paste(unique(stats::na.omit(kegg_results)), collapse = "; ")),
.groups = "drop"
)

annotated_set <- gene_df %>%
dplyr::left_join(combined_annots, by = c("Gene" = "Original"))

# --- Save annotated set (same rows as up/down file) ---
annotated_file <- file.path(comp_folder, paste0(set_label, "_annotated.txt"))
utils::write.table(annotated_set, annotated_file, sep = "\t", row.names = FALSE, quote = FALSE)

# ==========================================================================
# 6) Merge annotation columns into the FULL volcano table
#    IMPORTANT: merge only annotation columns to avoid duplicated stats columns
# ==========================================================================
ann_cols <- annotated_set %>%
dplyr::select(Gene, GO_Pathway, Reactome_Pathway, kegg_results) %>%
dplyr::mutate(Annotation_Label = "+") %>%
dplyr::distinct()

merged_matrix <- Whole_data %>%
dplyr::left_join(ann_cols, by = "Gene")

merged_file <- file.path(comp_folder, paste0("Merged_matrix_", set_label, ".txt"))
utils::write.table(merged_matrix, merged_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("✅ Annotation complete for", set_label, "in", comp_name, "\n")
}
}

cat("\n🎉 Annotation complete for all volcano comparison sets.\n")



################################################################################
################################################################################
################################################################################

# ==========================================================
# Fisher Exact Test Enrichment for Volcano Results (Cleaned for new naming)
# ==========================================================

library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)

# --- Working directory ---
setwd("C:/Users/ga53hil/Desktop/Granit_proteomics/13.02.26_Scientific_Data")

# --- Background file ---
bg_file <- "P337_02B_proteinGroups_filtered_pepGT1_log2LFQ_valid70_imputed_fixedseed_quantilenorm.txt"
background_data <- read.table(bg_file, header = TRUE, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
background_data$Gene_name[background_data$Gene_name == "" | is.na(background_data$Gene_name)] <- "Unknown"
background_data$Gene_name <- make.unique(background_data$Gene_name)
background_genes <- unique(background_data$Gene_name)

# --- Volcano results root ---
volcano_root <- file.path(getwd(), "Volcano_Results_limma_BH")

# --- Fisher test function ---
fisher_test_annotation <- function(df, source_col, sig_genes, background_genes, output_file) {
df <- df %>%
filter(!is.na(.data[[source_col]]) & .data[[source_col]] != "") %>%
mutate(Pathway = strsplit(as.character(.data[[source_col]]), ";")) %>%
unnest(Pathway) %>%
mutate(Pathway = str_trim(Pathway)) %>%
filter(!is.na(Pathway) & Pathway != "" & Pathway != "NA")

# --- Remove "NA | NA | NA" from GO ---
df <- df[!grepl("^NA\\s*\\|\\s*NA\\s*\\|\\s*NA$", df$Pathway, ignore.case = TRUE), ]

# --- Clean pathway display ---
df <- df %>%
mutate(Pathway_Name = case_when(
source_col == "kegg_results" ~ sub(" - Bos.*", "", Pathway),
grepl("\\|", Pathway) ~ sub(".*\\|\\s*", "", Pathway),
TRUE ~ Pathway
))

fisher_results <- data.frame()

for (i in unique(df$Pathway_Name)) {
term_genes <- unique(df$Gene[df$Pathway_Name == i & df$Gene %in% background_genes])
sig_in_pathway <- intersect(term_genes, sig_genes)

a <- length(sig_in_pathway)
b <- length(sig_genes) - a
c <- length(term_genes) - a
d <- length(background_genes) - a - b - c

mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
if (any(mat < 0) || any(!is.finite(mat))) next

test <- tryCatch(fisher.test(mat), error = function(e) NULL)
if (!is.null(test)) {
gene_ratio <- ifelse((a + b) > 0, a / (a + b), 0)
bg_ratio   <- ifelse((a + b + c + d) > 0, (a + c) / (a + b + c + d), 0)
enrichment_factor <- ifelse(bg_ratio > 0, gene_ratio / bg_ratio, NA)

fisher_results <- rbind(fisher_results, data.frame(
Pathway = i,
Sig_protein_volcano = a,
Sig_NotIn_Pathway   = b,
Nonsig_In_Pathway   = c,
Nonsig_NotIn_Pathway = d,
Genes = paste(sig_in_pathway, collapse = ";"),
p_value = test$p.value,
enrichment_factor = enrichment_factor
))
}
}

if (nrow(fisher_results) > 0) {
fisher_results$p_adj <- p.adjust(fisher_results$p_value, method = "BH")
fisher_results <- fisher_results[order(fisher_results$p_adj), ]
write.xlsx(fisher_results, output_file, rowNames = FALSE)
}
}

# --- Loop through Volcano comparison folders ---
comparisons <- list.dirs(volcano_root, recursive = FALSE, full.names = TRUE)

for (comp in comparisons) {
comp_name <- basename(comp)

# Find only new annotated up/down files (T#_up_T#_down_T#_vs_T#_annotated.txt)
annotated_files <- list.files(comp, 
pattern = "^T\\d+_up_T\\d+_down_T\\d+_vs_T\\d+_annotated\\.txt$", 
full.names = TRUE)

if (!length(annotated_files)) {
cat("❌ No annotated volcano sets for", comp_name, "\n")
next
}

for (annot_file in annotated_files) {
set_label <- tools::file_path_sans_ext(basename(annot_file))
cat("🔍 Running Fisher for:", set_label, "in", comp_name, "\n")

ann_df <- read.table(annot_file, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
if (!"Gene" %in% colnames(ann_df)) {
cat("❌ No Gene column in", annot_file, "\n")
next
}

sig_genes <- unique(ann_df$Gene)

# Output folder
out_dir <- file.path(comp, "Fisher")
dir.create(out_dir, showWarnings = FALSE)

# Run Fisher for GO, KEGG, Reactome
if ("GO_Pathway" %in% colnames(ann_df)) {
fisher_test_annotation(ann_df, "GO_Pathway", sig_genes, background_genes,
file.path(out_dir, paste0("Fisher_GO_", set_label, ".xlsx")))
}
if ("kegg_results" %in% colnames(ann_df)) {
fisher_test_annotation(ann_df, "kegg_results", sig_genes, background_genes,
file.path(out_dir, paste0("Fisher_KEGG_", set_label, ".xlsx")))
}
if ("Reactome_Pathway" %in% colnames(ann_df)) {
fisher_test_annotation(ann_df, "Reactome_Pathway", sig_genes, background_genes,
file.path(out_dir, paste0("Fisher_Reactome_", set_label, ".xlsx")))
}
}
}

cat("\n🎉 Fisher exact test complete for all Volcano results.\n")





################################################################################
################################################################################


################################################################################
################################################################################
################################################################################
# Volcano Enrichment Heatmaps (GO + KEGG + Reactome) [MERGED SOURCES]
# - Merges pathway names with same text, shows all sources in []
# - Filters: adj p < 0.01, protein count >= 5
# - Readable output with paging and keyword filter
# - AUTO x-axis gene label size + AUTO angle + AUTO plot width/height (per page)
#   so gene text stays readable across wide ranges of significant proteins.
################################################################################
################################################################################
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(openxlsx)
  library(viridisLite)
  library(scales)
})

# ---- User Parameters ----
volcano_root <- "Volcano_Results_limma_BH"
global_keywords <- c(
  # Metabolism
  "amino acid metabolism", "translation initiation", "rRNA processing",
  "mitochondrial translation", "lipid metabolism", "carbohydrate metabolism",
  "central carbon metabolism", "nucleotide metabolism", "glucose metabolism",
  "pyruvate metabolism", "Citrate cycle (TCA cycle)", "glycolysis", "gluconeogenesis",
  "beta oxidation", "fatty acid degradation", "steroid metabolism",
  "cholesterol metabolism", "glutathione metabolism", "methionine metabolism",
  "polyamine metabolism", "biosynthesis of amino acids",
  
  # Signaling
  "WNT signaling", "MAPK signaling", "interleukin signaling", "TNF signaling",
  "JAK-STAT signaling", "TGF-beta signaling", "interferon signaling",
  "NOTCH signaling", "EGFR signaling", "HIF-1 signaling", "TP53 signaling",
  "TCR signaling", "B cell receptor", "chemokine signaling",
  
  # Stress & Apoptosis
  "UPR", "autophagy", "apoptosis", "programmed cell death",
  "cellular senescence", "hypoxia response", "DNA repair", "DNA replication",
  "mitotic cell cycle", "cell cycle checkpoints",
  
  # Immune / Matrix / Transport
  "neutrophil degranulation", "extracellular matrix", "platelet activation",
  "lysosome", "endosome", "vesicle-mediated transport", "deubiquitination",
  "proteolysis", "biological oxidations", "organelle biogenesis",
  "ECM-receptor interaction",
  
  # Development / Morphogenesis
  "gastrulation", "somitogenesis", "chromatin organization", "histone modification",
  
  # Reproduction-specific
  "oocyte", "spermatogenesis", "meiosis", "gamete generation",
  "gonad development", "reproductive system", "reproductive development",
  "sexual reproduction", "fertilization", "androgen", "estrogen", "luteinizing hormone",
  "follicle-stimulating hormone", "hormone signaling"
)

P_ADJ_CUTOFF <- 0.01
MIN_PROTEINS <- 5
MAX_PATHWAYS_PER_PAGE <- 15
MAX_GENES_PER_PAGE <- 60
PATHWAY_WRAP_WIDTH <- 40

# Base figure size (used as minimums; actual size is auto per page)
A4_W <- 16
A4_H <- 10

BASE_FONTSIZE <- 14
AXIS_FONTSIZE_Y <- 20

LEGEND_TITLE_SIZE <- 16
LEGEND_TEXT_SIZE  <- 20
LEGEND_KEY_HEIGHT_PT <- 12
LEGEND_KEY_WIDTH_PT  <- 10

# ---- Helper functions ----
first_col <- function(df, cand) {
  c <- cand[cand %in% names(df)]
  if (length(c)) c[1] else NA_character_
}

count_genes <- function(x) {
  if (is.null(x) || is.na(x)) return(0L)
  toks <- unlist(stringr::str_split(as.character(x), "[;,/\\s]+"))
  toks <- stringr::str_trim(toks)
  toks <- toks[toks != ""]
  length(unique(toks))
}

# AUTO x-axis font size:
# - anchored at 30 genes => 30pt
# - scales down smoothly for more genes
# - never goes below 10pt (still readable)
auto_axis_fontsize_x <- function(n_genes, ref_genes = 30, ref_size = 30,
                                 min_size = 10, max_size = 30) {
  size <- ref_size * ref_genes / max(1, n_genes)
  size <- min(size, max_size)
  size <- max(size, min_size)
  round(size, 1)
}

# AUTO x-axis angle based on gene count (helps avoid overlap)
auto_angle_x <- function(n_genes) {
  if (n_genes > 50) 75 else if (n_genes > 35) 60 else 45
}

# AUTO plot width (inches) based on genes shown on the page
# This is the key to "gene text all can read across all range".
auto_plot_width_in <- function(n_genes,
                               base_w = A4_W,
                               min_w  = 16,
                               max_w  = 36,
                               inches_per_gene = 0.55) {
  w <- max(base_w, n_genes * inches_per_gene)
  w <- min(max(w, min_w), max_w)
  round(w, 1)
}

# Optional: AUTO plot height based on pathways on the page
auto_plot_height_in <- function(n_paths,
                                base_h = A4_H,
                                min_h  = 10,
                                max_h  = 18,
                                inches_per_path = 0.6) {
  h <- max(base_h, n_paths * inches_per_path)
  h <- min(max(h, min_h), max_h)
  round(h, 1)
}

# ---- Find all Fisher/ folders ----
all_dirs    <- list.dirs(volcano_root, recursive = TRUE, full.names = TRUE)
fisher_dirs <- all_dirs[grepl("[/\\\\]Fisher$", all_dirs)]

if (!length(fisher_dirs)) {
  stop("❌ No 'Fisher' folders found under: ", normalizePath(volcano_root))
}

# ---- Main loop: each Fisher directory (comparison) ----
for (fdir in fisher_dirs) {
  cat("\n📂 Processing:", fdir, "\n")
  
  files <- list.files(
    fdir,
    pattern = "^Fisher_(GO|KEGG|Reactome)_T\\d+_up_T\\d+_down_T\\d+_vs_T\\d+_annotated\\.xlsx$",
    full.names = TRUE
  )
  if (!length(files)) {
    cat("⚠️  No annotated Fisher enrichment files found in", fdir, "\n")
    next
  }
  
  # Parse file info for direction/comparison grouping
  info <- stringr::str_match(
    basename(files),
    "^Fisher_(GO|KEGG|Reactome)_(T\\d+_up_T\\d+_down_T\\d+_vs_T\\d+)_annotated\\.xlsx$"
  )
  colnames(info) <- c("full","Source","Direction")
  if (any(is.na(info))) next
  
  directions <- unique(info[, "Direction"])
  
  for (direction in directions) {
    cat("   🔎 Direction:", direction, "\n")
    
    f_direction <- files[info[, "Direction"] == direction]
    src         <- info[info[, "Direction"] == direction, "Source"]
    
    # Parse comp for titling (kept for filenames)
    direction_parts <- stringr::str_match(direction, "^(T\\d+)_up_(T\\d+)_down_(T\\d+_vs_T\\d+)$")
    comp <- direction_parts[4]
    
    pieces <- list()
    
    for (k in seq_along(f_direction)) {
      fpath       <- f_direction[k]
      source_type <- src[k]
      
      df0 <- tryCatch(openxlsx::read.xlsx(fpath), error = function(e) NULL)
      if (is.null(df0) || !"Pathway" %in% names(df0)) {
        cat("    ⚠️  Skipping (no Pathway):", basename(fpath), "\n")
        next
      }
      
      pcol <- first_col(df0, c("p_adj", "p.adj", "padj", "adj_p_value", "adj.P.Val", "p_adjust"))
      gcol <- first_col(df0, c("Genes", "Gene", "Gene_List", "GeneID", "Gene_Names", "Proteins", "Protein_List"))
      if (is.na(pcol) || is.na(gcol)) {
        cat("    ⚠️  Skipping (missing p-adj or gene column):", basename(fpath), "\n")
        next
      }
      
      df <- df0 %>%
        dplyr::mutate(
          .p         = suppressWarnings(as.numeric(.data[[pcol]])),
          .genes_raw = as.character(.data[[gcol]]),
          .gene_n    = vapply(.genes_raw, count_genes, integer(1))
        ) %>%
        dplyr::filter(is.finite(.p), .p < P_ADJ_CUTOFF, .gene_n >= MIN_PROTEINS) %>%
        dplyr::filter(
          stringr::str_detect(
            tolower(Pathway),
            paste(tolower(global_keywords), collapse = "|")
          )
        ) %>%
        dplyr::mutate(
          Pathway = stringr::str_trim(Pathway),
          Source  = source_type,
          log10_p = -log10(.p)
        ) %>%
        tidyr::separate_rows(dplyr::all_of(gcol), sep = "[;,/\\s]+") %>%
        dplyr::mutate(Gene = stringr::str_trim(.data[[gcol]])) %>%
        dplyr::filter(Gene != "")
      
      pieces[[length(pieces) + 1]] <- df
    }
    
    side_df <- dplyr::bind_rows(pieces)
    if (!nrow(side_df)) {
      cat("   ⚠️  No GO/KEGG/Reactome rows for", direction, "after filtering\n")
      next
    }
    
    # ---- Combine same-named pathways and collect sources ----
    source_levels <- c("KEGG", "GO", "Reactome")
    side_df <- side_df %>%
      dplyr::mutate(
        Pathway_Only = stringr::str_remove(Pathway, "\\s*\\[.*\\]$"),
        Source = factor(Source, levels = source_levels)
      ) %>%
      dplyr::group_by(Pathway_Only, Gene) %>%
      dplyr::summarise(
        log10_p = max(log10_p, na.rm = TRUE),
        SourceList = paste(intersect(source_levels, unique(as.character(Source))), collapse = ", "),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        Pathway = paste0(
          stringr::str_wrap(Pathway_Only, width = PATHWAY_WRAP_WIDTH),
          " [", SourceList, "]"
        )
      )
    
    # ---- Limit genes for readability: most frequent across pathways ----
    gene_order_by_freq <- side_df %>% dplyr::count(Gene, sort = TRUE) %>% dplyr::pull(Gene)
    keep_genes <- head(gene_order_by_freq, MAX_GENES_PER_PAGE)
    side_df <- side_df %>% dplyr::filter(Gene %in% keep_genes)
    
    mat <- side_df %>%
      dplyr::select(Pathway, Gene, log10_p) %>%
      tidyr::complete(Pathway, Gene = keep_genes)
    
    # ---- Pathway ordering: by number of proteins (highest first) ----
    ord <- mat %>%
      dplyr::filter(!is.na(log10_p)) %>%
      dplyr::group_by(Pathway) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(n)) %>%
      dplyr::pull(Pathway)
    
    mat <- mat %>%
      dplyr::mutate(
        Pathway = factor(Pathway, levels = ord),
        Gene    = factor(Gene, levels = keep_genes)
      )
    
    out_dir <- file.path(fdir, "Combined_Enrichment_Plots", direction)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    # ---- Page split for output ----
    pths_all <- levels(mat$Pathway)
    chunks <- split(pths_all, ceiling(seq_along(pths_all) / MAX_PATHWAYS_PER_PAGE))
    
    for (i in seq_along(chunks)) {
      sel <- chunks[[i]]
      
      d <- mat %>%
        dplyr::filter(Pathway %in% sel) %>%
        dplyr::mutate(Pathway = factor(Pathway, levels = rev(sel)))
      
      # ---- AUTO sizing per page ----
      n_genes_page <- length(unique(d$Gene))
      n_paths_page <- length(unique(d$Pathway))
      
      AXIS_FONTSIZE_X_AUTO <- auto_axis_fontsize_x(n_genes_page)
      ANGLE_X_AUTO <- auto_angle_x(n_genes_page)
      
      # Key improvement: auto width/height so genes are readable across ranges
      W_IN <- auto_plot_width_in(n_genes_page)
      H_IN <- auto_plot_height_in(n_paths_page)
      
      caption_txt <- paste0(
        "Filters: adj p < ", P_ADJ_CUTOFF,
        "; proteins \u2265 ", MIN_PROTEINS,
        "; genes shown capped at ", MAX_GENES_PER_PAGE,
        " by frequency."
      )
      
      p <- ggplot(d, aes(x = Gene, y = Pathway, fill = log10_p)) +
        geom_tile(color = "white", linewidth = 0.25, na.rm = FALSE) +
        scale_fill_gradientn(
          colours = c("green", "red"),
          na.value = "black",
          name = "-log10(p.adj)"
        ) +
        scale_x_discrete(position = "top", guide = guide_axis(check.overlap = TRUE)) +
        # Remove title + x label, keep gene tick labels
        labs(title = NULL, x = NULL, y = "Enriched Pathway", caption = caption_txt) +
        theme_minimal(base_size = BASE_FONTSIZE) +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          
          axis.text.x = element_text(
            angle = ANGLE_X_AUTO,
            hjust = 0,
            vjust = 0.1,
            size  = AXIS_FONTSIZE_X_AUTO,
            margin = margin(b = 2)
          ),
          axis.text.y = element_text(size = AXIS_FONTSIZE_Y),
          
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = AXIS_FONTSIZE_Y, margin = margin(r = 8)),
          
          plot.title = element_blank(),
          
          legend.position = "right",
          legend.title = element_text(size = LEGEND_TITLE_SIZE),
          legend.text  = element_text(size = LEGEND_TEXT_SIZE),
          legend.key.height = unit(LEGEND_KEY_HEIGHT_PT, "pt"),
          legend.key.width  = unit(LEGEND_KEY_WIDTH_PT,  "pt"),
          
          plot.caption = element_text(size = 10, hjust = 1),
          plot.margin = margin(t = 16, r = 16, b = 16, l = 16)
        ) +
        coord_cartesian(clip = "off", expand = FALSE)
      
      out_base <- file.path(out_dir, paste0("Enrichment_", direction, "_", comp, "_Page", i))
      
      ggsave(paste0(out_base, ".pdf"),
             plot = p, width = W_IN, height = H_IN, units = "in",
             device = cairo_pdf, bg = "white")
      
      ggsave(paste0(out_base, ".png"),
             plot = p, width = W_IN, height = H_IN, units = "in",
             dpi = 600, bg = "white", limitsize = FALSE)
      
      cat("   ✅ Saved:", out_base,
          "| genes:", n_genes_page,
          "| x-font:", AXIS_FONTSIZE_X_AUTO,
          "| angle:", ANGLE_X_AUTO,
          "| W:", W_IN, "H:", H_IN, "\n")
    }
  }
}

cat("\n🎉 All enrichment plots generated with merged pathway sources and readability-optimized auto sizing.\n")



################################################################################
################################################################################

# ----------------------------------------------------------
# Volcano Sankey: Top 5 Pathways PER Timepoint (Timepoint-specific)
# Keeps all your filters:
#  - keyword filter (global_keywords)
#  - dedup via Pathway_key (lower/trim/squish)
#  - FDR_CUTOFF and MIN_PROTEINS
# Key change vs your current script:
#  - selection is done PER timepoint, and we KEEP ONLY those pathway-timepoint pairs
#    (so no “shared pathway across many TPs unless it is top in those TPs too”)
# ----------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(openxlsx)
  library(RColorBrewer)
  library(networkD3)
  library(htmlwidgets)
  library(webshot2)
  library(jsonlite)
})

cat("🔎 Starting Volcano Sankey: Timepoint-specific Top 5 Per T1-T10\n")

volcano_root <- file.path(getwd(), "Volcano_Results_limma_BH")
if (!dir.exists(volcano_root)) stop("❌ 'Volcano_Results_limma_BH' folder not found at: ", volcano_root)

timepoints <- paste0("T", 1:10)
tp_colors  <- RColorBrewer::brewer.pal(10, "Paired")
names(tp_colors) <- timepoints

FDR_CUTOFF   <- 0.01
MIN_PROTEINS <- 5
TOP_N        <- 10   # <---- timepoint-specific top N

# --- Pathway normalization helpers (KEY FIX) ---
normalize_pathway <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("[\u00A0]", " ") %>%  # non-breaking spaces -> normal space
    stringr::str_squish() %>%                      # trim + collapse multiple spaces
    stringr::str_to_lower()                        # case-insensitive key
}
pretty_pathway <- function(key) {
  stringr::str_to_title(key)
}

global_keywords <- c(
  # Metabolism
  "amino acid metabolism", "translation initiation", "rRNA processing",
  "mitochondrial translation", "lipid metabolism", "carbohydrate metabolism",
  "central carbon metabolism", "nucleotide metabolism", "glucose metabolism",
  "pyruvate metabolism", "Citrate cycle (TCA cycle)", "glycolysis", "gluconeogenesis",
  "beta oxidation", "fatty acid degradation", "steroid metabolism",
  "cholesterol metabolism", "glutathione metabolism", "methionine metabolism",
  "polyamine metabolism", "biosynthesis of amino acids",
  
  # Signaling
  "WNT signaling", "MAPK signaling", "interleukin signaling", "TNF signaling",
  "JAK-STAT signaling", "TGF-beta signaling", "interferon signaling",
  "NOTCH signaling", "EGFR signaling", "HIF-1 signaling", "TP53 signaling",
  "TCR signaling", "B cell receptor", "chemokine signaling",
  
  # Stress & Apoptosis
  "UPR", "autophagy", "apoptosis", "programmed cell death",
  "cellular senescence", "hypoxia response", "DNA repair", "DNA replication",
  "mitotic cell cycle", "cell cycle checkpoints",
  
  # Immune / Matrix / Transport
  "neutrophil degranulation", "extracellular matrix", "platelet activation",
  "lysosome", "endosome", "vesicle-mediated transport", "deubiquitination",
  "proteolysis", "biological oxidations", "organelle biogenesis",
  "ECM-receptor interaction",
  
  # Development / Morphogenesis
  "gastrulation", "somitogenesis", "chromatin organization", "histone modification",
  
  # Reproduction-specific
  "oocyte", "spermatogenesis", "meiosis", "gamete generation",
  "gonad development", "reproductive system", "reproductive development",
  "sexual reproduction", "fertilization", "androgen", "estrogen", "luteinizing hormone",
  "follicle-stimulating hormone", "hormone signaling"
)

file_pattern <- "^Fisher_(GO|KEGG|Reactome)_(T\\d+)_up_(T\\d+)_down_(T\\d+_vs_T\\d+)_annotated\\.xlsx$"

# --- Helper: Read Fisher and tidy ---
read_and_tidy <- function(filepath) {
  bn <- basename(filepath)
  m  <- stringr::str_match(bn, file_pattern)
  if (any(is.na(m))) {
    cat("⚠️  Skipping unmatched filename: ", bn, "\n")
    return(NULL)
  }
  
  src      <- m[2]
  first_tp <- m[3]
  second_tp <- m[4]
  comp     <- m[5]
  side     <- if (startsWith(comp, first_tp)) "First" else "Second"
  up_in    <- first_tp
  
  df <- tryCatch(openxlsx::read.xlsx(filepath), error = function(e) NULL)
  if (is.null(df)) {
    cat("❌ Could not read file: ", bn, "\n")
    return(NULL)
  }
  
  gene_col <- if ("Genes" %in% names(df)) "Genes" else if ("Gene" %in% names(df)) "Gene" else NA_character_
  if (is.na(gene_col) || !"Pathway" %in% names(df)) {
    cat("⚠️  File lacks Pathway or Gene(s): ", bn, "\n")
    return(NULL)
  }
  
  if ("p_adj" %in% names(df)) df$p_adj <- suppressWarnings(as.numeric(df$p_adj)) else df$p_adj <- NA_real_
  
  df2 <- df %>%
    dplyr::mutate(
      Source = src, Up_in = up_in, Side = side, Comparison = comp,
      Pathway = stringr::str_squish(as.character(Pathway)),
      Pathway_key  = normalize_pathway(Pathway),
      Pathway_disp = pretty_pathway(Pathway_key)
    ) %>%
    dplyr::filter(!is.na(.data[[gene_col]]), .data[[gene_col]] != "") %>%
    dplyr::filter(
      !is.na(Pathway_key), Pathway_key != "",
      !(Pathway_key %in% c("na", "na | na | na"))
    )
  
  if (nrow(df2) == 0) {
    cat("⚠️  File has no usable rows after filter: ", bn, "\n")
    return(NULL)
  }
  
  # keyword filter: apply to normalized key to avoid missing due to case
  if (!is.null(global_keywords) && length(global_keywords) > 0) {
    pat <- paste0("(", paste0(global_keywords, collapse = "|"), ")")
    df2 <- df2 %>%
      dplyr::filter(stringr::str_detect(Pathway_key, stringr::str_to_lower(pat)))
    if (nrow(df2) == 0) {
      cat("⚠️  No keyword-matched pathways in: ", bn, "\n")
      return(NULL)
    }
  }
  
  df_long <- df2 %>%
    tidyr::separate_rows(.data[[gene_col]], sep = "[;,/\\s]+") %>%
    dplyr::mutate(Gene = stringr::str_trim(.data[[gene_col]])) %>%
    dplyr::filter(Gene != "") %>%
    dplyr::distinct(Pathway_key, Pathway_disp, Gene, Up_in, Side, Source, Comparison, p_adj)
  
  if (nrow(df_long) == 0) {
    cat("⚠️  No usable long-form rows: ", bn, "\n")
    return(NULL)
  }
  
  df_long <- df_long %>%
    dplyr::mutate(Timepoint = factor(Up_in, levels = timepoints))
  
  return(df_long)
}

# --- Load Fisher files ---
cat("🔎 Scanning all Fisher enrichment files...\n")
all_fisher_dirs <- list.dirs(volcano_root, recursive = TRUE, full.names = TRUE)
fdirs <- all_fisher_dirs[grepl(paste0(.Platform$file.sep, "Fisher$"), all_fisher_dirs)]
if (length(fdirs) == 0) stop("❌ No Fisher folders found under: ", volcano_root)
for (d in fdirs) cat("   Scanning:", normalizePath(d), "\n")

all_rows <- list()
for (d in fdirs) {
  files <- list.files(d, pattern = file_pattern, full.names = TRUE)
  if (length(files) == 0) next
  for (f in files) {
    r <- read_and_tidy(f)
    if (!is.null(r)) all_rows[[length(all_rows) + 1]] <- r
  }
}
if (length(all_rows) == 0) stop("❌ No enrichment rows collected. Check file names, contents or keyword filter.")
dat <- dplyr::bind_rows(all_rows)

# --- Deduplicate (USING Pathway_key) ---
dat <- dat %>%
  dplyr::group_by(Pathway_key, Pathway_disp, Gene, Timepoint, Side, Source, Comparison) %>%
  dplyr::summarise(
    p_adj = ifelse(all(is.na(p_adj)), NA_real_, min(p_adj, na.rm = TRUE)),
    .groups = "drop"
  )

# ----------------------------------------------------------
# TIMEPOINT-SPECIFIC SELECTION (Key change)
# 1) Compute pathway stats per timepoint
# 2) Apply MIN_PROTEINS and FDR_CUTOFF
# 3) Pick TOP_N per timepoint
# 4) Keep ONLY those Pathway_key + Timepoint combinations in the Sankey
# ----------------------------------------------------------

path_stats_tp <- dat %>%
  dplyr::group_by(Pathway_key, Timepoint) %>%
  dplyr::summarise(
    min_p = ifelse(all(is.na(p_adj)), NA_real_, min(p_adj, na.rm = TRUE)),
    gene_count = dplyr::n_distinct(Gene),
    Pathway_disp = dplyr::first(Pathway_disp),
    .groups = "drop"
  ) %>%
  dplyr::filter(gene_count >= MIN_PROTEINS, !is.na(min_p), min_p < FDR_CUTOFF) %>%
  dplyr::mutate(tp_num = as.numeric(stringr::str_replace(Timepoint, "T", "")))

top_tp <- path_stats_tp %>%
  dplyr::group_by(Timepoint) %>%
  dplyr::arrange(desc(gene_count), min_p) %>%
  dplyr::slice_head(n = TOP_N) %>%
  dplyr::ungroup()

tp_selected <- top_tp %>%
  dplyr::select(Timepoint, Pathway_key) %>%
  dplyr::distinct()

cat(sprintf("✅ Selected up to %d pathways per timepoint after filters (FDR<%.3g, proteins>=%d).\n",
            TOP_N, FDR_CUTOFF, MIN_PROTEINS))

# Keep only rows where pathway was selected for THAT timepoint
dat_sel <- dat %>%
  dplyr::inner_join(tp_selected, by = c("Timepoint", "Pathway_key"))

# Node ordering: by timepoint then strength
pathways_ordered_disp <- top_tp %>%
  dplyr::arrange(tp_num, desc(gene_count), min_p) %>%
  dplyr::pull(Pathway_disp) %>%
  unique()

# --- Sankey weighting (still your style), now inherently timepoint-specific ---
path_counts <- dat_sel %>%
  dplyr::group_by(Pathway_key) %>%
  dplyr::summarise(n_link = dplyr::n_distinct(Timepoint), .groups = "drop")

dat_w <- dat_sel %>%
  dplyr::left_join(path_counts, by = "Pathway_key") %>%
  dplyr::group_by(Pathway_key, Pathway_disp, Timepoint, Side, Comparison) %>%
  dplyr::summarise(
    weight = sum(1 / ifelse(is.na(n_link) | n_link == 0, 1, n_link)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(Timepoint = factor(as.character(Timepoint), levels = timepoints))

# --- Sankey plot function ---
plot_tp_specific_sankey <- function(df, label, outdir, pathways_disp_ordered) {
  if (is.null(df) || nrow(df) == 0) {
    cat("❌ No data for ", label, "\n")
    return(NULL)
  }
  
  pathways  <- pathways_disp_ordered
  tps_order <- timepoints
  
  nodes <- rbind(
    data.frame(name = pathways, type = "Pathway", stringsAsFactors = FALSE),
    data.frame(name = tps_order,  type = "Timepoint", stringsAsFactors = FALSE)
  )
  
  df <- df %>%
    dplyr::mutate(
      Pathway_disp = factor(Pathway_disp, levels = pathways),
      Timepoint    = factor(Timepoint, levels = tps_order)
    )
  
  links <- df %>%
    dplyr::mutate(
      source = match(as.character(Pathway_disp), nodes$name) - 1,
      target = match(as.character(Timepoint),    nodes$name) - 1,
      value  = weight,
      group  = as.character(Timepoint)
    ) %>%
    dplyr::select(source, target, value, group)
  
  # Keep missing timepoints visible (layout stabilizer)
  missing_tps <- setdiff(tps_order, unique(as.character(df$Timepoint)))
  if (length(missing_tps) > 0) {
    dummy_links <- data.frame(
      source = 0,
      target = match(missing_tps, nodes$name) - 1,
      value  = 0.0001,
      group  = missing_tps,
      stringsAsFactors = FALSE
    )
    links <- dplyr::bind_rows(links, dummy_links)
  }
  
  js_colors <- sprintf(
    'd3.scaleOrdinal().domain(%s).range(%s)',
    jsonlite::toJSON(timepoints),
    jsonlite::toJSON(unname(tp_colors))
  )
  
  sankey <- sankeyNetwork(
    Links = links,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value  = "value",
    NodeID = "name",
    fontSize = 36,
    nodeWidth = 36,
    sinksRight = TRUE,
    LinkGroup = "group",
    colourScale = js_colors,
    width = 1200, height = 800,
    nodePadding = 20,
    iterations = 0
  )
  
  outfile <- file.path(outdir, paste0("Sankey_", label, ".html"))
  htmlwidgets::saveWidget(sankey, file = outfile, selfcontained = TRUE, background = "#fff")
  
  # CSS to remove node rectangles and center text
  lines <- readLines(outfile, warn = FALSE)
  css_patch <- '
<style>
.node rect { fill: none !important; stroke: none !important; }
.node text {
  font-weight: bold;
  fill: #222;
  text-anchor: middle !important;
  x: 0 !important;
  transform: none !important;
  font-size: 32px !important;
}
body, html { background: #fff !important; }
</style>
'
idx <- grep("</head>", lines, fixed = TRUE)
if (length(idx) > 0) {
  lines <- append(lines, css_patch, after = idx[1] - 1)
  writeLines(lines, outfile)
}

cat("✅ Saved Sankey HTML: ", outfile, "\n")

outfile_png <- file.path(outdir, paste0("Sankey_", label, ".png"))
webshot2::webshot(
  url = outfile,
  file = outfile_png,
  vwidth = 1200,
  vheight = 750,
  delay = 3,
  zoom = 4
)
cat("✅ Also saved as PNG: ", outfile_png, "\n")

invisible(sankey)
}

# --- Output directory and plot ---
out_root <- file.path(volcano_root, "Sankey_Top5_TimepointSpecific_T1toT10_Manuscript")
if (!dir.exists(out_root)) dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

cat("🔎 Plotting Timepoint-specific Top pathways per Timepoint...\n")
plot_tp_specific_sankey(dat_w, paste0("Top", TOP_N, "_PerTimepoint_T1toT10"), out_root, pathways_ordered_disp)

cat("\n🎉 Sankey plot (Timepoint-specific, T1–T10, keyword-filtered, FDR<0.01 & proteins>=MIN_PROTEINS) saved in:\n",
    normalizePath(out_root), "\n")
cat("🎉 Done.\n")
