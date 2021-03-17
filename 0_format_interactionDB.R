#load packages
require('tidyverse')
require('dplyr')

#Kumar et al Cell Rep. 2018 Nov 6;25(6):1458-1468.e4. doi: 10.1016/j.celrep.2018.10.047.
kumar <- read.csv("../../../receptor_ligand_plus_strand/receptor_ligand_Ramilowski_modified.csv", header=TRUE)
dim(kumar)
#2567

#remove any duplicate entries
interaction.db <- kumar[!duplicated(kumar$Pair.Name),]
dim(interaction.db)
#2567
colnames(interaction.db)[2] <- "From"
colnames(interaction.db)[4] <- "To"

#include only literature supported
interaction.db <- interaction.db[interaction.db$Pair.Evidence == "literature supported",]
save.image(file="interaction_combinedCells.RData")
