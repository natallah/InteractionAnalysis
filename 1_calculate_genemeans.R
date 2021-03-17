
# > find_genes.seurat(Seurat::pbmc_small, "CD3E", assay = "RNA")
# [1] "CD3E"
find_genes.seurat <- function(seurat.obj, genes, database, assay = "rna") {
  stopifnot(assay %in% names(seurat.obj@assays))
  names_in_assay <- rownames(seurat.obj@assays[[assay]]@data)

  genes_found <- genes[genes %in% names_in_assay]
  return(genes_found)
}



# > genemeans(Seurat::pbmc_small, c("CD3E", "CD4"))
#   cluster gene mean_count
# 1       0 CD3E  2.2786085
# 2       1 CD3E  0.1417045
# 3       2 CD3E  1.7344702
#
# > genemeans(Seurat::pbmc_small, c("CD3E", "CD4", "SomeFaeGene"))

genemeans <- function(seurat.obj, genes_use, silent = FALSE, assay = 'rna') {
  if (silent) {
      FetchDataLocal <- function(...) {
          suppressMessages({
              suppressWarnings({
                  Seurat::FetchData(...)
              })
          })
      }
  } else {
      FetchDataLocal <- Seurat::FetchData
  }

  Seurat::DefaultAssay(seurat.obj) <- assay

  seurat_df <- FetchDataLocal(seurat.obj, vars = c(genes_use, "ident"))
  seurat_df <- dplyr::group_by(seurat_df, ident)
  seurat_df <- dplyr::summarise_each(seurat_df, mean)
  seurat_df <- dplyr::rename(seurat_df, cluster = ident)
  seurat_df <- reshape2::melt(seurat_df, variable.name = "gene", value.name = "mean_count")
  return(seurat_df)
}


# > countcells(Seurat::pbmc_small, c("CD3E", "S100A8"), assay = "RNA")
#   cluster   gene cell_ratio
# 1       0   CD3E     0.4500
# 2       1   CD3E     0.3125
# 3       2   CD3E     0.2375
# 4       0 S100A8     0.4500
# 5       1 S100A8     0.3125
# 6       2 S100A8     0.2375
countcells <- function(seurat.obj, genes_use, silent = FALSE, assay = 'rna') {
  if (silent) {
      FetchDataLocal <- function(...) {
          suppressMessages({
              suppressWarnings({
                  Seurat::FetchData(...)
              })
          })
      }
  } else {
      FetchDataLocal <- Seurat::FetchData
  }

  Seurat::DefaultAssay(seurat.obj) <- assay

  seurat_df <- FetchDataLocal(seurat.obj, vars = c(genes_use, "ident"))
  num_cells <- nrow(seurat_df)
  seurat_df <- dplyr::group_by(seurat_df, ident)
  seurat_df <- dplyr::summarise_each(seurat_df, length)
  seurat_df <- dplyr::rename(seurat_df, cluster = ident)
  seurat_df <- reshape2::melt(seurat_df, variable.name = "gene", value.name = "cell_ratio")
  seurat_df <- dplyr::mutate(seurat_df, cell_ratio = cell_ratio/num_cells) 
  return(seurat_df)
}
