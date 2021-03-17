
#' Calculate the interaction weights
#' 
#' @param interaction.db a data.frame with columns To and From
#' @param genemeans a data.frame with columns cluster, gene and mean_count
#'
#' > calculate_interactions(head(interaction.db), genemeans_list[[1]])
#'    Gene.from Gene.to Cluster.from Cluster.to interaction_score
#'  1       A2M    LRP1            0          0      0.0006330224
#'  2       A2M    LRP1            0          1      0.0000000000
#'  3       A2M    LRP1            0          2      0.0004072501
#'  4       A2M    LRP1            0          3      0.0136198115
#'  5       A2M    LRP1            0          4      0.0342077661
#'  6       A2M    LRP1            0          5      0.0000000000
calculate_interactions <- function(interaction.db, genemeans) {
  slim_interaction_db <- dplyr::select(interaction.db, From, To)

  tmp_from <- left_join(interaction.db, genemeans, by = c("From" = "gene"))
  tmp_to <- left_join(interaction.db, genemeans, by = c("To" = "gene"))

  tmp_join <- full_join(tmp_from, tmp_to,
            by = c("From" = "From", "To" = "To"), 
            suffix = c(".from", ".to"))

  tmp_join$interaction_score <- tmp_join$mean_count.from * tmp_join$mean_count.to
  tmp_join <- tmp_join[!is.na(tmp_join$interaction_score),]

  return(dplyr::select(
      tmp_join,
      Gene.from = From, Gene.to = To,
      Cluster.from = cluster.from, Cluster.to = cluster.to,
      interaction_score))
}


# Calculates null interactions for a receptor, randomizing the ligands
calculate_null_interactions_receptor <- function(interaction.db, genemeans, receptor, nperm = 100) {
  interaction.db <- interaction.db[interaction.db$From == receptor]

  interaction.db$From <- sample(interaction.db$From)
  calculate_interactions(interaction.db, genemeans)
}



calculate_null_interactions_receptor_all <- function(interaction.db, genemeans, nperm = 10) {
  main <- function() {
      interaction.db$From <- sample(interaction.db$From)
      return(calculate_interactions(interaction.db, genemeans))
  }

  purrr::map_df(
      seq_len(nperm),
       ~ { print(.x) ; main() }
  )
}

calculate_interaction_signif <- function(interaction_score_tbl, null_interaction_score) {

  

}
