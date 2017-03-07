library(Paschold2012)
library(dplyr)

str(Paschold2012)
counts = Paschold2012 %>%
  mutate(genotype_replicate = paste(genotype,replicate,sep="_")) %>%
  select(GeneID, genotype_replicate, total) %>%
  tidyr::spread(genotype_replicate, total)

zero_rows <- counts$GeneID[which(sapply(1:nrow(counts), function(row) sum(counts[row,-1]) == 0))]
y <- filter(counts, !(GeneID %in% zero_rows)) %>%
  select(-GeneID) %>%
  as.matrix()

saveRDS(y, "filtered_counts.rds")