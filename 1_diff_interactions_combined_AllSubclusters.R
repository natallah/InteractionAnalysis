library(stats)
require(Seurat)

source("../0_format_interactionDB.R")
source("../1_calculate_genemeans.R")
source("../2_calculate_interactions.R")

kumar <- read.csv("../../../receptor_ligand_plus_strand/receptor_ligand_Ramilowski_modified.csv", header=TRUE)
dim(kumar)
#2567

#remove any duplicate entries
interaction.db <- kumar[!duplicated(kumar$Pair.Name),]
dim(interaction.db)
#2567
colnames(interaction.db)[2] <- "From"
colnames(interaction.db)[4] <- "To"
interaction.db <- interaction.db[interaction.db$Pair.Evidence == "literature supported",]

seurat.obj <- readRDS(file = "../../CC_corrected_combined_no766_may2020.RDS")

sample_names <- unique(seurat.obj@meta.data$Sample)

seurat.obj@meta.data %>% head(2)

DimPlot(seurat.obj, label = TRUE)

seurat.obj@meta.data %>% .[,c("Sample")] %>% table()

gsub(".*(small|large|strand).*", "\\1", seurat.obj[["Sample"]][[1]] ) %>% unique() %>% head()
# [1] "small"  "large"  "strand"
seurat.obj[["Sample_Type"]] <- gsub(".*(small|large|strand).*", "\\1", seurat.obj[["Sample"]][[1]] )

#label subclustering clusters
seurat.obj[["Subclusters"]] <- seurat.obj@meta.data$seurat_clusters
meta <- data.frame(as.character(seurat.obj@meta.data$Subclusters))
rownames(meta) <- rownames(seurat.obj@meta.data)
meta$as.character.seurat.obj.meta.data.Subclusters. <- as.character(meta$as.character.seurat.obj.meta.data.Subclusters.)
plasma <- readRDS(file="../../subclustering/plasma_cells_no766_noCD3CD8_RNA.RDS")
plasma.meta <- data.frame(as.character(plasma@meta.data$seurat_clusters))
rownames(plasma.meta)<- rownames(plasma@meta.data)
colnames(plasma.meta) <- "plasma_subclusters"
plasma.meta$plasma_subclusters <- paste ("plasma", plasma.meta$plasma_subclusters,
                                                       sep="_")
mast <- readRDS(file="../../subclustering/mast_cells_no766_noKRT_noCD3CD8_RNA.RDS")
mast.meta <- data.frame(as.character(mast@meta.data$seurat_clusters))
rownames(mast.meta)<- rownames(mast@meta.data)
colnames(mast.meta) <- "mast_subclusters"
mast.meta$mast_subclusters <- paste ("mast", mast.meta$mast_subclusters,
                                         sep="_")
Tcells <- readRDS(file="../../subclustering/T_cells_no766_All_clustered_no11_12_13_noCD19.RDS")
t.meta <- data.frame(as.character(Tcells@meta.data$seurat_clusters))
rownames(t.meta)<- row.names(Tcells@meta.data)
colnames(t.meta) <- "T_subclusters"
t.meta$T_subclusters <- paste ("T", t.meta$T_subclusters,
                                     sep="_")
test<- sub("\\_\\d+$","",row.names(Tcells@meta.data))

mac <- readRDS(file="../../subclustering/mac_cells_no766_nocluster3_7_8_9_11.RDS")
mac.meta <- data.frame(as.character(mac@meta.data$seurat_clusters))
rownames(mac.meta)<- rownames(mac@meta.data)
colnames(mac.meta) <- "mac_subclusters"
mac.meta$mac_subclusters <- paste ("mac", mac.meta$mac_subclusters,
                                   sep="_")
nk <- readRDS(file="../../subclustering/NK_cells_no766_noKT.RDS")
nk.meta <- data.frame(as.character(nk@meta.data$seurat_clusters))
rownames(nk.meta)<- rownames(nk@meta.data)
colnames(nk.meta) <- "NK_subclusters"
nk.meta$NK_subclusters <- paste ("NK", nk.meta$NK_subclusters,
                                   sep="_")
b <- readRDS(file="../../subclustering/B_cells_no766_noKT_no6_noCD3CD8_RNA.RDS")
b.meta <- data.frame(as.character(b@meta.data$seurat_clusters))
rownames(b.meta)<- rownames(b@meta.data)
colnames(b.meta) <- "B_subclusters"
b.meta$B_subclusters <- paste ("B", b.meta$B_subclusters,
                                 sep="_")

testmeta <- meta
testmatch <- match(rownames(plasma.meta),rownames(meta))
testmeta[testmatch,] <- plasma.meta[,1]
testmatch2 <- match(rownames(mast.meta),rownames(testmeta))
testmeta[testmatch2,] <- mast.meta[,1]
testmatch3 <- match(rownames(mac.meta),rownames(testmeta))
testmeta[testmatch3,] <- mac.meta[,1]
testmatch4 <- match(sub("\\_\\d+$","",row.names(Tcells@meta.data)),rownames(testmeta))
testmeta[testmatch4,] <- t.meta[,1]
testmatch5 <- match(rownames(nk.meta),rownames(testmeta))
testmeta[testmatch5,] <- nk.meta[,1]
testmatch6 <- match(rownames(b.meta),rownames(testmeta))
testmeta[testmatch6,] <- b.meta[,1]

#write.csv(testmeta,"metaData_Dec2020.csv")

seurat.obj[["Subclusters"]] <- testmeta
seurat.obj <- SetIdent(seurat.obj,value = "Subclusters")

genes_in_db <- unique(c(as.character(interaction.db$From),
                        as.character(interaction.db$To)))
genes_in_assay <- rownames(seurat.obj@assays$RNA)
#genes_in_assay <- rownames(seurat.obj@assays[[seurat.obj@active.assay]]@scale.data) 
genes_in_db2 <- genes_in_db[genes_in_db %in% genes_in_assay]
genes_in_db3 <- genes_in_db2


diff_between_samples <- read_tsv("diff_between_samples_padj0.05_allSubclusters.tsv") %>% as_tibble()

genes_in_db3 <- genes_in_db2[genes_in_db2 %in% diff_between_samples$gene]

diff_vln_list <- {
diff_between_samples %>%
  filter(gene %in% genes_in_db3) %>% 
  group_by(Comparisson) %>%
  top_n(6, avg_logFC) %>%
  select(gene) %>% 
  nest() %>% ungroup() %>% { 
     purrr::map2(
        .$data, gsub("(\\d+).*", "\\1", .$Comparisson),
        ~ VlnPlot(seurat.obj, .x$gene, group.by = "Sample_Type", pt.size = 0, ident = .y))
  }
}

split_seurat_obje <- Seurat::SplitObject(
        seurat.obj,
        split.by = "Sample")

names(split_seurat_obje)


genemeans_list <- lapply(
    split_seurat_obje,
    genemeans, 
    genes_use = genes_in_db3, 
    silent = TRUE,
    assay = "SCT")

head(genemeans_list[[1]])

cellcount_list <- lapply(
    split_seurat_obje,
    countcells, 
    genes_use = genes_in_db3, 
    silent = TRUE,
    assay = "SCT")

head(cellcount_list[[1]])
#   cluster gene cell_ratio
# 1       0  A2M 0.05766063
# 2       1  A2M 0.15183965

weighted_genemeans_list <- purrr::map2(genemeans_list, cellcount_list, ~ {
        ret <- mutate(full_join(.x, .y), weighted_mean_count = mean_count*cell_ratio)
        stopifnot(nrow(ret) == nrow(.x))
        return(ret)
    })


# USING ONLY GENE MEANS, NOT THE WEIGHTED

purrr::map_df(weighted_genemeans_list, ~ as_tibble(.x)[, "mean_count"], .id = "file") %>% as_tibble() %>%
  ggplot(aes(x = log1p(mean_count * 10000), colour = file, fill = file)) + geom_density(alpha = 0.01) 


interactions_list <- purrr::map(
    weighted_genemeans_list,
    function(x, interaction.db) {
      calculate_interactions(mutate(x, mean_count = mean_count),
                             interaction.db = interaction.db)
    },
    interaction.db = interaction.db
)

names(interactions_list) <- names(split_seurat_obje)

interactions_df <- purrr::map_df(interactions_list, ~ .x, .id = "Sample") %>% as_tibble()



require(tidyr)
nonzero <- interactions_df[interactions_df$interaction_score>0,]

# This section completes all interactions without a calculated score as 0
{
interactions_df <- complete(
  interactions_df, 
  Sample,
  nesting(Cluster.from, Cluster.to),
  nesting(Gene.from, Gene.to))

interactions_df[is.na(interactions_df$interaction_score),]

interactions_df$interaction_score[is.na(interactions_df$interaction_score)] <- 0

interactions_df$interaction_score  <- interactions_df$interaction_score
}

dim(interactions_df)
# [1] 43904     6

(length(unique(interactions_df$Cluster.from)) ^ 2) *
(nrow(unique(interactions_df[,c("Gene.from", "Gene.to")])))  *
(length(unique(interactions_df$Sample)))


# This function calculates the wilcox test for a chunk uf the former data frame
get_wilcox <- function(chunk, filter_regex = ".*(small|large).*") {
    filter(chunk, grepl(filter_regex, Sample)) %>%
    mutate(sample_type = gsub(filter_regex, "\\1",Sample)) %>%
    wilcox.test(interaction_score ~ sample_type, .)
}

get_t <- function(chunk, filter_regex = ".*(small|large).*") {
    filter(chunk, grepl(filter_regex, Sample)) %>%
    mutate(sample_type = gsub(filter_regex, "\\1",Sample)) %>%
    t.test(interaction_score ~ sample_type, .)
}

library(furrr)
plan(multiprocess)

nested_interaction_df <- group_by(
      interactions_df,
      Gene.from,
      Gene.to,
      Cluster.from,
      Cluster.to) %>% 
      nest()

table(purrr::map_dbl(nested_interaction_df$data, ~ nrow(.x)))
# 
#   17 
# 4961 

## Here I am filtering interactions that have only zeros in
## all samples 
all_zeros <- purrr::map_lgl(nested_interaction_df$data, ~ all(.x$interaction_score == 0)) 
table(all_zeros)
# all_zeros
# FALSE  TRUE 
# 1514252 2568873

nested_interaction_df <- nested_interaction_df[!all_zeros,]
nested_interaction_df
# # A tibble: 4,568 x 5
# # Groups:   Cluster.from, Cluster.to, Gene.from, Gene.to [4,568]
#    Cluster.from Cluster.to Gene.from Gene.to           data
#    <chr>        <chr>      <chr>     <chr>   <list<df[,2]>>
#  1 0            0          ANXA1     FPR1          [17 × 2]
#  2 0            0          ANXA1     FPR3          [17 × 2]
#  3 0            0          C1QA      CD93          [17 × 2]
#  4 0            0          C3        C3AR1         [17 × 2]
#  5 0            0          C3        ITGAX         [17 × 2]
#  6 0            0          CCL20     CCR6          [17 × 2]
#  7 0            0          CCL22     CCR4          [17 × 2]
#  8 0            0          CCL3      CCR4          [17 × 2]
#  9 0            0          CCL5      CCR4          [17 × 2]
# 10 0            0          CCL5      SDC1          [17 × 2]
# # … with 4,558 more rows

nested_interaction_df$data[[1]]
# # A tibble: 17 x 2
#    Sample            interaction_score
#    <chr>                         <dbl>
#  1 003_small                     0.828
#  2 004_small                     0.334
#  3 006_small                     0.348
#  4 007_small                     0.711
#  5 008_small                     0.546
#  6 009_small                     0.347
#  7 010_small                     0.857
#  8 012_large                     0.435
#  9 013_small                     1.41 
# 10 1144_small                    0.644
# 11 1157_large                    0.767
# 12 1195_large                    0.190
# 13 1196_small                    0.155
# 14 766_large                     0    
# 15 sample_strand_D17             0    
# 16 sample_strand_D27             0    
# 17 sample_strand_D35             0    


nested_interaction_df$data[[1]] %>% get_wilcox(".*(small|large).*")
# 
# 	Wilcoxon rank sum test
# 
# data:  interaction_score by sample_type
# W = 12, p-value = 0.3037
# alternative hypothesis: true location shift is not equal to 0
# 

nested_interaction_df$data[[1]] %>% get_wilcox(".*(small|strand).*")
# 
# 	Wilcoxon rank sum test with continuity correction
# 
# data:  interaction_score by sample_type

#  = 25.5, p-value = 0.07541
# alternative hypothesis: true location shift is not equal to 0
# 

combin_comparissons <- combn(c("strand","small", "large"), 2)
#      [,1]     [,2]     [,3]   
# [1,] "strand" "strand" "small"
# [2,] "small"  "large"  "large"

safe_w <- purrr::safely(get_wilcox)
safe_t <- purrr::safely(get_t) 
options(future.globals.maxSize = 200000*10^24)
for (i in seq_len(ncol(combin_comparissons))) {
    comp <- print(paste(combin_comparissons[, i], collapse = "|"))
# [1] "small|large"
    col_name <- make.names(comp)
# [1] "small.large"
    regex <- paste0(".*(", comp, ").*")
# [1] ".*(small|large).*"
    nested_interaction_df[[paste0("wilcox_test_", col_name)]] <- 
      furrr::future_map(
      nested_interaction_df$data,
      safe_w,
      .progress = TRUE, regex)


    nested_interaction_df[[paste0("t_test_", col_name)]] <- furrr::future_map(
      nested_interaction_df$data,
      safe_t,
      .progress = TRUE, regex)
}

colnames(nested_interaction_df)

for (i in grep(".*_test_.*", colnames(nested_interaction_df), value = TRUE)) {
  errors <- purrr::map(nested_interaction_df[[i]], ~ .x$error)
  head(errors)
  tot_errors <- sum(purrr::map_lgl(errors, ~ !is.null(.x)))
  message(tot_errors, " errors found in col ", i)
}

nested_interaction_df$t_test_strand.large[[1]][["result"]][["estimate"]]
#  mean in group large mean in group strand 
#            0.3479374            0.0000000 


comparisson_names <- grep(".*_test_.*", colnames(nested_interaction_df), value = TRUE) %>%
  gsub(".*_", "", .) %>% 
  unique()

for (comparisson in comparisson_names) {
  wilcox_name <- paste0("wilcox_test_", comparisson)
  t_name <- paste0("t_test_", comparisson)

  t_coef_names <- purrr::map(
    nested_interaction_df[[t_name]],
    ~ names(.x[["result"]][["estimate"]]))  %>% unlist() %>% unique()

  g <- {
    qplot(purrr::map_dbl(
        nested_interaction_df[[t_name]],
        ~ .x$result$estimate[[ t_coef_names[[1]] ]] - .x$result$estimate[[ t_coef_names[[2]] ]]),
      -log10(purrr::map_dbl(
          nested_interaction_df[[wilcox_name]],
          ~ .x$result$p.value)), geom = "jitter") + 
    geom_hline(yintercept = -log10(0.05), colour = "red") + 
    geom_vline(xintercept = 1.5, colour = "red") + 
    geom_vline(xintercept = -1.5, colour = "red") + 
    labs(x = paste(t_coef_names[[1]], t_coef_names[[2]], sep = " - "), y = "Wilcox Test: -log10(p.val)") +
    theme_bw()
  }
  ggsave(g, filename = paste0(comparisson, "_volcano_wilcox_test.png") , width = 6, height = 7)

  g <- {
    qplot(purrr::map_dbl(
        nested_interaction_df[[t_name]],
        ~ .x$result$estimate[[ t_coef_names[[1]] ]] - .x$result$estimate[[ t_coef_names[[2]] ]]),
      -log10(purrr::map_dbl(
          nested_interaction_df[[t_name]],
          ~ .x$result$p.value)), geom = "jitter") + 
    geom_hline(yintercept = -log10(0.05), colour = "red") + 
    geom_vline(xintercept = 0.2, colour = "red") + 
    geom_vline(xintercept = -0.2, colour = "red") + 
    labs(x = paste(t_coef_names[[1]], t_coef_names[[2]], sep = " - "), y = "T Test: -log10(p.val)") +
    theme_bw()
  }
  ggsave(g, filename = paste0(comparisson, "_volcano_T_test.png") , width = 6, height = 7)
}


results_df <- tibble(
  Gene.from = nested_interaction_df$Gene.from,
  Gene.to = nested_interaction_df$Gene.to,
  Cluster.from = nested_interaction_df$Cluster.from,
  Cluster.to = nested_interaction_df$Cluster.to)

for (comparisson in comparisson_names) {
  wilcox_name <- paste0("wilcox_test_", comparisson)
  t_name <- paste0("t_test_", comparisson)

  t_coef_names <- purrr::map(
    nested_interaction_df[[t_name]],
    ~ names(.x[["result"]][["estimate"]]))  %>% unlist() %>% unique()

  slim_coef_names <- gsub(".* ", "", t_coef_names)
  comp_name <- paste(slim_coef_names[[1]], slim_coef_names[[2]], sep = "_vs_")

  diff_estimate <- purrr::map_dbl(
    nested_interaction_df[[t_name]],
    ~ .x$result$estimate[[ t_coef_names[[1]] ]] - .x$result$estimate[[ t_coef_names[[2]] ]])

  wilcox_pval <- purrr::map_dbl(
    nested_interaction_df[[wilcox_name]],
    ~ .x$result$p.value)

  t_pval <- purrr::map_dbl(
    nested_interaction_df[[t_name]],
    ~ .x$result$p.value)

  results_df[[comp_name]] <- diff_estimate
  results_df[[paste0("w_pval_", comp_name)]] <- wilcox_pval
  results_df[[paste0("t_pval_", comp_name)]] <- t_pval

}

results_df

results_df <- {
  results_df %>%
    mutate(
        label = paste0(Gene.from, "->",
                       Gene.to, "::",
                       Cluster.from, "->", 
                       Cluster.to)) 
   }


for (comparisson in comparisson_names) {
  wilcox_name <- paste0("wilcox_test_", comparisson)
  t_name <- paste0("t_test_", comparisson)

  t_coef_names <- purrr::map(
    nested_interaction_df[[t_name]],
    ~ names(.x[["result"]][["estimate"]]))  %>% unlist() %>% unique()

  slim_coef_names <- gsub(".* ", "", t_coef_names)
# [1] "small"  "strand"
  comp_name <- paste(slim_coef_names[[1]], slim_coef_names[[2]], sep = "_vs_")
# [1] "small_vs_strand"
  signif_results_df <- filter(results_df, abs(results_df[[comp_name]]) > 0.2, results_df[[paste0("w_pval_", comp_name)]] < 0.05)

  labeldata <- unique(rbind(
     top_n(signif_results_df %>% ungroup(), 10, !! sym(comp_name)),
     top_n(signif_results_df %>% ungroup(), -10, !! sym(comp_name))))

  stopifnot(nrow(labeldata) < 30)

  g <- {
    results_df %>%
    ggplot(aes_string(
        y = paste0("-log10(", paste0("w_pval_", comp_name), ")"),
        x = comp_name, label = "label")) +
      geom_point() +
      ggrepel::geom_label_repel(data = labeldata)
  }

  ggsave(g, filename = paste0(comparisson, "_labelled_volcano_wilcox_test.png") , width = 6, height = 7)
}
results_df
# # A tibble: 4,568 x 14
#    Gene.from Gene.to Cluster.from Cluster.to small_vs_strand t_pval_small_vs… w_pval_small_vs… large_vs_strand t_pval_large_vs… w_pval_large_vs…
#    <chr>     <chr>   <chr>        <chr>                <dbl>            <dbl>            <dbl>           <dbl>            <dbl>            <dbl>
#  1 ANXA1     FPR1    0            0                 0.00618           0.0137         0.000445         0.00348            0.119            0.126 
#  2 ANXA1     FPR3    0            0                 0.00156           0.0455         0.00263          0.00287            0.119            0.309 

#adjust p-values for multiple testing

results_df$w_padj_small_vs_strand <- p.adjust(results_df$w_pval_small_vs_strand,
                                                method= "fdr",n=length(results_df$w_pval_small_vs_strand))
results_df$w_padj_large_vs_strand <- p.adjust(results_df$w_pval_large_vs_strand,
                                                method= "fdr",n=length(results_df$w_pval_large_vs_strand))
results_df$w_padj_large_vs_small <- p.adjust(results_df$w_pval_large_vs_small,
                                                method= "fdr",n=length(results_df$w_pval_large_vs_small))

readr::write_csv(results_df, "all_results_df_dec2020.csv")

signif_results_df <- results_df %>% filter_at(vars(matches("w_padj")), any_vars(. < 0.10 ))
# # A tibble: 1,788 x 14
#    Gene.from Gene.to Cluster.from Cluster.to small_vs_strand t_pval_small_vs… w_pval_small_vs… large_vs_strand t_pval_large_vs… w_pval_large_vs…
#    <chr>     <chr>   <chr>        <chr>                <dbl>            <dbl>            <dbl>           <dbl>            <dbl>            <dbl>
#  1 ANXA1     FPR1    0            0                 0.00618           0.0137         0.000445         0.00348            0.119            0.126 
#  2 ANXA1     FPR3    0            0                 0.00156           0.0455         0.00263          0.00287            0.119            0.309 
#  3 C1QA      CD93    0            0                 0.000264          0.0137         0.00420          0.000955           0.0436           0.121 
#  4 C3        C3AR1   0            0                 0.000854          0.0137         0.00852          0.000937           0.119            0.0946
#  5 C3        ITGAX   0            0                 0.000223          0.0137         0.0185           0.000349           0.119            0.0978
#  6 CCL20     CCR6    0            0                 0.0224            0.00699        0.000543         0.0102             0.285            0.169 
#  7 CCL22     CCR4    0            0                 0.000629          0.0137         0.0171           0.000720           0.119            0.124 
#  8 CCL3      CCR4    0            0                 0.0289            0.0141         0.0000324        0.0327             0.0497           0.0154
#  9 CCL5      CCR4    0            0                 0.201             0.0141         0.0000156        0.218              0.0497           0.0143
# 10 CD28      CD86    0            0                 0.000960          0.0137         0.00112          0.00102            0.119            0.0893
# # … with 1,778 more rows, and 4 more variables: large_vs_small <dbl>, t_pval_large_vs_small <dbl>, w_pval_large_vs_small <dbl>, label <chr>





non_auto_signif_results_df <- filter(signif_results_df, Cluster.from != Cluster.to)
readr::write_csv(non_auto_signif_results_df, "non_auto_signif_results_df_dec2020.csv")

#T-Cells to macs or vice versa (mac_0 through mac_9 and T_0 through T_9)

mac_T_non_auto_signif_results_df <- non_auto_signif_results_df %>% filter_at(vars(matches("Cluster")), all_vars(str_detect(.,pattern= "mac_|T_" )))
readr::write_csv(mac_T_non_auto_signif_results_df, "mac_T_non_auto_signif_results_df_dec2020.csv")

mac_to_X_non_auto_signif_results_df <- non_auto_signif_results_df %>% filter_at(vars(matches("Cluster.From")), all_vars(str_detect(.,pattern= "mac_" )))
mac_to_mast_non_auto_signif_results_df <- mac_to_X_non_auto_signif_results_df %>% filter_at(vars(matches("Cluster.To")), all_vars(str_detect(.,pattern= "mast_" )))
readr::write_csv(mac_to_mast_non_auto_signif_results_df, "mac_to_mast_non_auto_signif_results_df_dec2020.csv")

mast_to_X_non_auto_signif_results_df <- non_auto_signif_results_df %>% filter_at(vars(matches("Cluster.From")), all_vars(str_detect(.,pattern= "mast_" )))
mast_to_mac_non_auto_signif_results_df <- mast_to_X_non_auto_signif_results_df %>% filter_at(vars(matches("Cluster.To")), all_vars(str_detect(.,pattern= "mac_" )))
readr::write_csv(mast_to_mac_non_auto_signif_results_df, "mast_to_mac_non_auto_signif_results_df_dec2020.csv")


g2 <- {
  mac_T_non_auto_signif_results_df %>%
    ggplot(aes_string(
      y = paste0("-log10(", paste0("w_pval_", comp_name), ")"),
      x = comp_name, label = "label")) +
    geom_point() +
    ggrepel::geom_label_repel(data = labeldata)
}
png(file="Subcluster_combined_dec2020.png", width=15, height=8, res=300, units='in')
DimPlot(seurat.obj, reduction="umap", group.by="Subclusters")
dev.off()
ggsave(g, filename = paste0(comparisson, "_labelled_volcano_wilcox_test_mac_Tcells_dec2020.png") , width = 6, height = 7)
# # A tibble: 1,590 x 14
#    Gene.from Gene.to Cluster.from Cluster.to small_vs_strand t_pval_small_vs… w_pval_small_vs… large_vs_strand t_pval_large_vs… w_pval_large_vs…
#    <chr>     <chr>   <chr>        <chr>                <dbl>            <dbl>            <dbl>           <dbl>            <dbl>            <dbl>
#  1 ANXA1     FPR1    0            1                0.732              0.00699       0.000268         0.426               0.0571          0.0281 
#  2 ANXA1     FPR3    0            1                0.265              0.0490        0.0239           0.291               0.114           0.0621 
#  3 C3        C3AR1   0            1                0.00449            0.0137        0.00425          0.00559             0.119           0.138  
#  4 C3        ITGAX   0            1                0.00838            0.0137        0.00833          0.00869             0.119           0.0809 
#  5 CCL20     CCR6    0            1                0.00465            0.0137        0.00102          0.00425             0.119           0.165  
#  6 CCL22     CCR4    0            1                0.0000239          0.0137        0.0123           0.0000366           0.119           0.145  
#  7 CD28      CD86    0            1                0.130              0.00699       0.00000117       0.102               0.0571          0.0614 
#  8 CD40LG    CD40    0            1                0.0617             0.00699       0.0000142        0.0462              0.285           0.144  
#  9 CSF1      CSF1R   0            1                0.0877             0.00699       0.00106          0.162               0.0571          0.0331 
# 10 CTLA4     CD86    0            1                0.108              0.00699       0.00000120       0.0748              0.0571          0.00529
# # … with 1,580 more rows, and 4 more variables: large_vs_small <dbl>, t_pval_large_vs_small <dbl>, w_pval_large_vs_small <dbl>, label <chr>

# FPlot_top_int <- function(Gene.from, Gene.to, Cluster.from, Cluster.to, label) {
#   g1 <- FeaturePlot(seurat.obj, c(Gene.from), 
#               split.by = "Sample Type",
#               label = TRUE,
#               order = TRUE,
#               col = viridis::viridis(7, direction =-1), 
#               cells = Cells(seurat.obj)[Idents(seurat.obj) %in% c(Cluster.from)])
# 
#   g2 <- FeaturePlot(seurat.obj, c(Gene.to), 
#               split.by = "Sample Type",
#               label = TRUE,
#               order = TRUE,
#               col = viridis::viridis(7, direction =-1), 
#               cells = Cells(seurat.obj)[Idents(seurat.obj) %in% c(Cluster.to)])
#   return(cowplot::plot_grid(g1, g2) + ggtitle(label))
# }
# 
# VPlot_top_int <- function(Gene.from, Gene.to, Cluster.from, Cluster.to, label) {
#   print(c(Gene.from, Gene.to, Cluster.from, Cluster.to, label))
#   idents <- print(as.character(c(Cluster.from, Cluster.to)))
#   g1 <- VlnPlot(seurat.obj, c(Gene.from), 
#               split.by = "Sample Type",
#               log = TRUE,
#               pt.size = 0,
#               idents = idents[[1]])
# 
#   g2 <- VlnPlot(seurat.obj, c(Gene.to), 
#               split.by = "Sample Type",
#               log = TRUE,
#               pt.size = 0,
#               idents = idents[[2]])
#               
#   return(cowplot::plot_grid(g1, g2) + ggtitle(label))
# 
# }
# 
# non_auto_signif_results_df %>% arrange(-abs(LGvsSM)) %>% top_n(10, (LGvsSM))
# # # A tibble: 10 x 8
# #    Gene.from Gene.to Cluster.from Cluster.to  t_pval  w_pval LGvsSM label              
# #    <chr>     <chr>   <fct>        <fct>        <dbl>   <dbl>  <dbl> <chr>              
# #  1 ICAM1     ITGAX   5            4          0.108   0.0240   1.47  ICAM1->ITGAX::5->4 
# #  2 ANXA1     FPR1    3            5          0.00767 0.00799 -0.869 ANXA1->FPR1::3->5  
# #  3 IL1A      IL1R2   9            2          0.0318  0.0240  -0.751 IL1A->IL1R2::9->2  
# #  4 ANXA1     FPR1    9            5          0.0147  0.0140  -0.688 ANXA1->FPR1::9->5  
# #  5 ANXA1     FPR1    1            5          0.0109  0.0240  -0.677 ANXA1->FPR1::1->5  
# #  6 IL1RN     IL1R1   1            12         0.0235  0.0282  -0.655 IL1RN->IL1R1::1->12
# #  7 IL1A      IL1R2   0            2          0.00909 0.00200 -0.619 IL1A->IL1R2::0->2  
# #  8 IL1RN     IL1R2   0            2          0.0312  0.00799 -0.608 IL1RN->IL1R2::0->2 
# #  9 IL1B      IL1R2   13           6          0.00742 0.0337  -0.604 IL1B->IL1R2::13->6 
# # 10 IL1A      IL1R1   3            12         0.0220  0.0282  -0.520 IL1A->IL1R1::3->12 
# 
# gglist <- non_auto_signif_results_df %>% arrange(-abs(LGvsSM)) %>% top_n(10, abs(LGvsSM)) %>%
#   select(-LGvsSM, -w_pval, -t_pval) %>%
#   mutate(Cluster.from = as.character(Cluster.from), Cluster.to = as.character(Cluster.to)) %>%
#   pmap(FPlot_top_int)
# 
# gglist[[1]]
# gglist[[2]]
# gglist
# 
# 
# gglist <- non_auto_signif_results_df %>% arrange(-abs(LGvsSM)) %>% top_n(10, (LGvsSM)) %>%
#   select(-LGvsSM, -w_pval, -t_pval) %>%
#   mutate(Cluster.from = as.character(Cluster.from), Cluster.to = as.character(Cluster.to)) %>%
#   pmap(VPlot_top_int)
# 
# gglist[[1]]
# gglist[[2]]
# gglist[[2]]
# 
# non_auto_signif_results_df %>% top_n(10, -(w_pval)) 
# # # A tibble: 14 x 8
# #    Gene.from Gene.to  Cluster.from Cluster.to    t_pval  w_pval  LGvsSM label              
# #    <chr>     <chr>    <fct>        <fct>          <dbl>   <dbl>   <dbl> <chr>              
# #  1 IL1A      IL1R2    0            2          0.00909   0.00200  -0.619 IL1A->IL1R2::0->2  
# #  2 IL1A      IL1R1    0            7          0.000245  0.00200  -7.13  IL1A->IL1R1::0->7  
# #  3 ANXA1     FPR1     0            9          0.0000819 0.00200 -45.8   ANXA1->FPR1::0->9  
# #  4 ANXA1     FPR2     0            9          0.0000774 0.00200 -37.1   ANXA1->FPR2::0->9  
# #  5 IL1A      IL1R1    0            9          0.0000976 0.00200 -11.4   IL1A->IL1R1::0->9  
# #  6 IL1A      IL1R2    0            9          0.0000317 0.00200 -15.4   IL1A->IL1R2::0->9  
# #  7 IL1B      IL1R2    0            9          0.00240   0.00200 -79.6   IL1B->IL1R2::0->9  
# #  8 TNF       TNFRSF1B 1            9          0.000125  0.00200 -46.0   TNF->TNFRSF1B::1->9
# #  9 ANXA1     FPR1     12           5          0.00152   0.00200  -1.02  ANXA1->FPR1::12->5 
# # 10 IL1B      IL1R1    0            9          0.000642  0.00400 -61.5   IL1B->IL1R1::0->9  
# # 11 IL1RN     IL1R2    0            9          0.00158   0.00400 -13.0   IL1RN->IL1R2::0->9 
# # 12 ANXA1     FPR2     3            9          0.00263   0.00400 -69.4   ANXA1->FPR2::3->9  
# # 13 ANXA1     FPR2     5            9          0.00242   0.00400 -93.8   ANXA1->FPR2::5->9  
# # 14 IL1A      IL1R1    9            0          0.00464   0.00400 -11.7   IL1A->IL1R1::9->0  
# 
# 
# gglist <- non_auto_signif_results_df %>% top_n(10, -(w_pval)) %>%
#   select(-LGvsSM, -w_pval, -t_pval) %>%
#   mutate(Cluster.from = as.character(Cluster.from), Cluster.to = as.character(Cluster.to)) %>%
#   pmap(FPlot_top_int)
# 
# gglist[[1]]
# gglist[[2]]
# gglist[[3]]
# gglist
# 
# gglist <- non_auto_signif_results_df %>% top_n(10, -(w_pval)) %>%
#   select(-LGvsSM, -w_pval, -t_pval) %>%
#   mutate(Cluster.from = as.character(Cluster.from), Cluster.to = as.character(Cluster.to)) %>%
#   pmap(VPlot_top_int)
# 
# gglist[[3]]
# gglist

################# ===================
# adj.pvals <- purrr::map_dbl(nested_interaction_df$t_tests, ~ .x$result$p.value) %>% p.adjust(method = "BH")
# table(adj.pvals < 0.05)
# # 
# #  FALSE   TRUE 
# # 125669     10 
# 
# nested_interaction_df[adj.pvals < 0.05,]
# 
# g <- {
#   labeldata <- results_df[adj.pvals < 0.05, ]
#   results_df %>%
#   ggplot(aes(y = -log10(t_pval), x = BPHvsNORM, label = label)) +
#     geom_point() +
#     ggrepel::geom_label_repel(data = labeldata)
# }
# ggsave(g, filename = "weighted_labelled_volcano_BH_ttest.png", width = 6, height = 7)
# 
# adj.pvals <- purrr::map_dbl(nested_interaction_df$wilcox_tests, ~ .x$result$p.value) %>% p.adjust(method = "BH")
# table(adj.pvals < 0.05)
# # 
# #  FALSE 
# # 125679 
# 
# adj.pvals %>% table()
# # .
# # 0.158522802810297 0.232736974246821 0.323964041045546 0.364266178173463 
# #             28359             15179                 1               151 
# # 0.379661171525447 0.381194446198859 0.392733303864871 0.490857471573108 
# #               569             12402                17                86 
# #  0.49662276250847 0.513727398408545 0.537485257276963 0.538426829895721 
# #             26524                41                 1               560 
# # 0.546832304756905 0.580336783388689 0.630756457967253 0.654109558454458 
# #               740               448                27             31249 
# # 0.684515107647935 0.704602487286336 0.734787348121943 0.772943216024556 
# #                 2                80              2655                41 
# # 0.789323909787501 0.811543369967584 0.902517296743695 0.908953770317918 
# #              1875               672                 2               104 
# #                 1 
# #              3894 
# ################# ===================
# 
# purrr::map_dbl(nested_interaction_df$wilcox_tests, ~ .x$result$p.value) %>% head()
# # [1] 0.03571429 0.30169958 0.03571429 0.30169958 0.32911399 0.03247332
# 
# purrr::map_dbl(nested_interaction_df$wilcox_tests, ~ .x$result$p.value) %>% table()
# # .
# #  0.016804604239087 0.0261791970636472 0.0324733203260696 
# #               3103                292              18301 
# # 0.0357142857142857 0.0357700822324908  0.056663942878564 
# #               3595               3071                 11 
# # 0.0667980684751383 0.0714285714285714 0.0719240820130502 
# #                 17                487                605 
# #  0.078982579263783 0.0806252626513428  0.112230924681785 
# #               4108              10959                  2 
# #  0.126630457947617  0.133701125808964  0.142857142857143 
# #                153                574                451 
# #  0.171857339062799  0.177112629766732  0.221699993764876 
# #              12944                 17                107 
# #  0.230243350577407               0.25  0.266236376195817 
# #                570                376                581 
# #  0.293825847189637  0.301699582478348  0.305507086861254 
# #                 16              11376                 24 
# #  0.329113985978608  0.340616892098009  0.359396770820519 
# #              15196                 39                583 
# #  0.368227133821696  0.392857142857143  0.427124088792106 
# #                772                483                 18 
# #  0.453355504597787  0.494524667883871  0.525234049633751 
# #                704                  4                 57 
# #  0.541192898741999  0.548595465799934  0.571428571428571 
# #                333               3075               2016 
# #  0.605576616335346  0.633737059218191  0.652782845261911 
# #              27599                  7                122 
# #  0.696270340114023  0.732678261375707  0.750751184646232 
# #               2794                 23                137 
# #  0.759981527157728  0.785714285714286  0.873807137111765 
# #               2106                647                  1 
# #  0.880791022511061                  1 
# #                 95               4583 
# 
# purrr::map_dbl(nested_interaction_df$wilcox_tests, ~ .x$result$p.value) %>% { table(. < 0.05) }
# # 
# #  FALSE   TRUE 
# # 104772  28362 
# 
# head(adj.pvals)
# # [1] 0.1679083 0.5040060 0.1679083 0.5040060 0.5040060 0.1679083
# 
# 
# 
interactions_df

#separate interactions based on sample type
late_interactions_df <- filter(interactions_df, Sample == "012_large"|Sample == "1157_large"|Sample=="1195_large"|Sample=="766_large")
early_interactions_df <- filter(interactions_df,Sample=="003_small"|Sample=="004_small"|Sample=="006_small"|Sample=="007_small"|Sample=="008_small"|Sample=="009_small"|Sample=="010_small"|Sample=="013_small"|Sample=="1144_small"|Sample=="1196_small")
normal_interactions_df <- filter(interactions_df,Sample=="sample_strand_D17"|Sample=="sample_strand_D27"|Sample=="sample_strand_D35")

median_normal_interactions_df <- group_by(
  normal_interactions_df, 
  Gene.from, Gene.to, 
  Cluster.from, Cluster.to) %>%
  summarise(median_int_score = median(interaction_score))

median_small_interactions_df <- group_by(
  early_interactions_df, 
  Gene.from, Gene.to, 
  Cluster.from, Cluster.to) %>%
  summarise(median_int_score = median(interaction_score))

median_large_interactions_df <- group_by(
  late_interactions_df, 
  Gene.from, Gene.to, 
  Cluster.from, Cluster.to) %>%
  summarise(median_int_score = median(interaction_score))
median_normal_interactions_df <- filter(median_normal_interactions_df, median_int_score > 0)
median_small_interactions_df <- filter(median_small_interactions_df, median_int_score > 0)
median_large_interactions_df <- filter(median_large_interactions_df, median_int_score > 0)

readr::write_csv(median_normal_interactions_df, "non_0_median_interactions_normal_allSubtypes_dec2020.csv")
readr::write_csv(median_small_interactions_df, "non_0_median_interactions_small_allSubtypes_dec2020.csv")
readr::write_csv(median_large_interactions_df, "non_0_median_interactions_large_allSubtypes_dec2020.csv")
