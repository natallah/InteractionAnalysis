
require(Seurat)
require(purrr)
require(tidyverse)


source("../0_format_interactionDB.R")
source("../1_calculate_genemeans.R")
source("../2_calculate_interactions.R")

# seurat.obj <- readRDS(file = "../step1_cellRangerCounts/cellranger_seurat/combinedOct14/combined_seurat.RDS")
# seurat.obj <- SCTransform(seurat.obj, vars.to.regress = "percent.mt", verbose = FALSE)
# saveRDS(file = "sctransformed_seurat.RDS", seurat.obj)

dir(pattern = "res")
# [1] "res_1.1_imm_clustered_combined_seurat.RDS"

seurat.obj <- readRDS(file = "../../CC_corrected_combined_no766_may2020.RDS")
seurat.obj
seurat.obj@meta.data %>% head()
#####test
get_condition_de(condition_col = "Sample_Type", ident_col = "seurat_clusters", seurat_obj = seurat.obj)

condition_combs <- combn(unique(seurat.obj[["Sample_Type"]][["Sample_Type"]]), 2)
condition_combs <- as.data.frame(condition_combs)
params <- expand.grid(
  condition.2 = seq_along(colnames(condition_combs)),
  idents = unique(as.character(seurat.obj[["seurat_clusters"]][["seurat_clusters"]]))
)
params <- as.data.frame(params)
params$condition.1 <- unlist(condition_combs[1, params$condition.2])
params$condition.2 <- unlist(condition_combs[2, params$condition.2])

print(params)

seurat.obj[["TMP_IDENTS"]] <- paste0(
  seurat.obj[["Sample_Type"]][["Sample_Type"]],
  "_", 
  Idents(seurat.obj))

Idents(seurat.obj) <- "TMP_IDENTS"


tmp_fun <- function(condition.1, condition.2, idents) {
  message(paste0("Running for: ", condition.1, 
                 " vs ", condition.2,
                 " cluster ", idents))
  
  diff_tbl <- FindMarkers(
    seurat.obj,
    ident.1 = paste0(condition.1, "_", idents),
    ident.2 = paste0(condition.2, "_", idents),
    verbose = FALSE)
  
  diff_tbl$gene <- rownames(diff_tbl)
  diff_tbl <- as_tibble(diff_tbl)
  
  return(diff_tbl)
}
  
####end test

get_condition_de <- function(condition_col, ident_col, seurat_obj) {
  condition_combs <- combn(unique(seurat_obj[[condition_col]][[condition_col]]), 2)
  condition_combs <- as.data.frame(condition_combs)


  params <- expand.grid(
    condition.2 = seq_along(colnames(condition_combs)),
    idents = unique(as.character(seurat_obj[[ident_col]][[ident_col]]))
  )

  params <- as.data.frame(params)
  params$condition.1 <- unlist(condition_combs[1, params$condition.2])
  params$condition.2 <- unlist(condition_combs[2, params$condition.2])

  print(params)

  seurat_obj[["TMP_IDENTS"]] <- paste0(
      seurat_obj[[condition_col]][[condition_col]],
      "_", 
      Idents(seurat_obj))

  Idents(seurat_obj) <- "TMP_IDENTS"

  tmp_fun <- function(condition.1, condition.2, idents) {
    message(paste0("Running for: ", condition.1, 
                   " vs ", condition.2,
                   " cluster ", idents))

    diff_tbl <- FindMarkers(
      seurat.obj,
      ident.1 = paste0(condition.1, "_", idents),
      ident.2 = paste0(condition.2, "_", idents),
      verbose = FALSE)

    diff_tbl$gene <- rownames(diff_tbl)
    diff_tbl <- as_tibble(diff_tbl)

    return(diff_tbl)
  }

  safe_fun <- safely(tmp_fun)

  outs <- purrr::pmap(params[, c("condition.1", "condition.2", "idents")], safe_fun)

  errors <- purrr::map_lgl(outs, ~ is.null(.x$error))

  if ( ! all(errors)) {
    message("Errors were found, returning list")
    names(outs) <- paste0(params$idents, "_", params$condition.1, "vs",params$condition.2)
    return(outs)
  } else {
    names(outs) <- paste0(params$idents, "_", params$condition.1, "vs",params$condition.2)
    outs <- purrr::map_df(outs, .id = "Comparisson", ~ .x$result)
    return(outs)
  }

}

gsub(".*(small|large|strand).*", "\\1", seurat.obj[["Sample"]][["Sample"]]) %>% unique() %>% head()
seurat.obj[["Sample_Type"]] <- gsub(".*(small|large|strand).*", "\\1", seurat.obj[["Sample"]][["Sample"]])

table(seurat.obj@meta.data[, c("seurat_clusters", "Sample_Type")])
#                 Sample_Type

table(seurat.obj@meta.data[, c("seurat_clusters", "Sample")])


foo <- get_condition_de(condition_col = "Sample_Type", ident_col = "seurat_clusters", seurat_obj = seurat.obj)

foo_lst <- foo
errors <- purrr::map_lgl(foo, ~ is.null(.x$error))
names(foo[!errors] )
# [1] "10_smallvsstrand" "10_largevsstrand" "6_smallvsstrand"  "6_largevsstrand" 
foo[!errors] 
# [[1]]
# [[1]]$result
# NULL
# 
# [[1]]$error
# <simpleError in FindMarkers.default(object = data.use, slot = data.slot, counts = counts,     cells.1 = ident.1, cells.2 = ident.2, features = features,     reduction = reduction, logfc.threshold = logfc.threshold,     test.use = test.use, min.pct = min.pct, min.diff.pct = min.diff.pct,     verbose = verbose, only.pos = only.pos, max.cells.per.ident = max.cells.per.ident,     random.seed = random.seed, latent.vars = latent.vars, min.cells.feature = min.cells.feature,     min.cells.group = min.cells.group, pseudocount.use = pseudocount.use,     ...): Cell group 2 has fewer than 3 cells>
# 
# 
# [[2]]
# [[2]]$result
# NULL
# 
# [[2]]$error
# <simpleError in FindMarkers.default(object = data.use, slot = data.slot, counts = counts,     cells.1 = ident.1, cells.2 = ident.2, features = features,     reduction = reduction, logfc.threshold = logfc.threshold,     test.use = test.use, min.pct = min.pct, min.diff.pct = min.diff.pct,     verbose = verbose, only.pos = only.pos, max.cells.per.ident = max.cells.per.ident,     random.seed = random.seed, latent.vars = latent.vars, min.cells.feature = min.cells.feature,     min.cells.group = min.cells.group, pseudocount.use = pseudocount.use,     ...): Cell group 2 has fewer than 3 cells>
# 
# 
# [[3]]
# [[3]]$result
# NULL
# 
# [[3]]$error
# <simpleError in WhichCells.Seurat(object = object, idents = ident.2): Cannot find the following identities in the object: strand_6>
# 
# 
# [[4]]
# [[4]]$result
# NULL
# 
# [[4]]$error
# <simpleError in WhichCells.Seurat(object = object, idents = ident.2): Cannot find the following identities in the object: strand_6>
# 
# 

foo <- purrr::map_df(foo[errors], ~ .x$result, .id = "Comparisson") %>% as_tibble()
foo
# # A tibble: 19,791 x 7
#    Comparisson        p_val avg_logFC pct.1 pct.2 p_val_adj gene    
#    <chr>              <dbl>     <dbl> <dbl> <dbl>     <dbl> <chr>   
#  1 1_smallvslarge 0.           -0.439 0.047 0.304 0.        IGLC3   
#  2 1_smallvslarge 0.           -0.497 0.046 0.334 0.        IGHG1   
#  3 1_smallvslarge 0.           -0.500 0.124 0.389 0.        IGHA1   
#  4 1_smallvslarge 0.           -0.798 0.269 0.617 0.        IGKC    
#  5 1_smallvslarge 0.           -0.879 0.246 0.628 0.        HLA-DRB5
#  6 1_smallvslarge 0.           -1.21  0.071 0.372 0.        IGLC2   
#  7 1_smallvslarge 5.01e-288     0.488 0.83  0.63  1.17e-283 ETS2    
#  8 1_smallvslarge 1.41e-284     0.283 0.999 0.997 3.30e-280 EIF1    
#  9 1_smallvslarge 1.96e-279     0.546 0.82  0.578 4.59e-275 NFKB1   
# 10 1_smallvslarge 7.20e-279     0.628 0.94  0.812 1.68e-274 SOD2    
# # … with 19,781 more rows


signif_changes <- foo %>% filter(p_val_adj < 0.05) %>% mutate(cluster = gsub("(\\d+).*", "\\1", Comparisson))
# # A tibble: 8,227 x 8
#    Comparisson        p_val avg_logFC pct.1 pct.2 p_val_adj gene     cluster
#    <chr>              <dbl>     <dbl> <dbl> <dbl>     <dbl> <chr>    <chr>  
#  1 1_smallvslarge 0.           -0.439 0.047 0.304 0.        IGLC3    1      
#  2 1_smallvslarge 0.           -0.497 0.046 0.334 0.        IGHG1    1      
#  3 1_smallvslarge 0.           -0.500 0.124 0.389 0.        IGHA1    1      
#  4 1_smallvslarge 0.           -0.798 0.269 0.617 0.        IGKC     1      
#  5 1_smallvslarge 0.           -0.879 0.246 0.628 0.        HLA-DRB5 1      
#  6 1_smallvslarge 0.           -1.21  0.071 0.372 0.        IGLC2    1      
#  7 1_smallvslarge 5.01e-288     0.488 0.83  0.63  1.17e-283 ETS2     1      
#  8 1_smallvslarge 1.41e-284     0.283 0.999 0.997 3.30e-280 EIF1     1      
#  9 1_smallvslarge 1.96e-279     0.546 0.82  0.578 4.59e-275 NFKB1    1      
# 10 1_smallvslarge 7.20e-279     0.628 0.94  0.812 1.68e-274 SOD2     1      
# # … with 8,217 more rows

write_tsv(signif_changes, "diff_between_samples_padj0.05_allSubclusters.tsv")

top_changes <- signif_changes %>% group_by(Comparisson) %>% top_n(10, abs(avg_logFC))
# # A tibble: 277 x 8
# # Groups:   Comparisson [28]
#    Comparisson        p_val avg_logFC pct.1 pct.2 p_val_adj gene     cluster
#    <chr>              <dbl>     <dbl> <dbl> <dbl>     <dbl> <chr>    <chr>  
#  1 1_smallvslarge 0.           -0.798 0.269 0.617 0.        IGKC     1      
#  2 1_smallvslarge 0.           -0.879 0.246 0.628 0.        HLA-DRB5 1      
# 10 1_smallvslarge 9.22e- 36     0.737 0.241 0.161 2.16e- 31 HAMP     1      
# # … with 267 more rows
