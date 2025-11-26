################################################################################
## Script Name:        Co-phylo with linked tips
## Purpose:            Plot multiple phylogenies side-by-side (identical taxa)
##                    and link tips together
## Author:             James Baxter
## Date Created:       2025-06-24
################################################################################

############################### SYSTEM OPTIONS #################################
options(
  scipen = 6,     # Avoid scientific notation
  digits = 7      # Set precision for numerical display
)
memory.limit(30000000)

############################### DEPENDENCIES ###################################
# Load required libraries
library(tidyverse)
library(magrittr)


################################### DATA #######################################
# Read and inspect data
s1 <- read.tree("./2025Jun24/global_analysis/USUV_2025Jun24_alldata_aligned_formatted_noFLI_NFLG_subsampled.fasta.contree")
s2 <- read.tree("./2025Jun24/global_analysis/USUV_2025Jun24_alldata_aligned_formatted_noFLI_ENV_subsampled.fasta.contree")
s3 <- read.tree("./2025Jun24/global_analysis/USUV_2025Jun24_alldata_aligned_formatted_noFLI_NS5_subsampled.fasta.contree")


################################### MAIN #######################################
# Main analysis or transformation steps
tree1 <- fortify(s1)
tree2 <- fortify(s2)
tree3 <- fortify(s3)

tree2$x <- tree2$x + max(tree1$x) +.1
tree3$x <- tree3$x + max(tree2$x) +.1

jingmen_s <- bind_rows(tree1, tree2, tree3) %>%
  filter(!is.na(label) & isTip == TRUE) %>%
  mutate(of_interest = grepl('^PQ1445', label))

ggtree(tree1) +
  geom_tree(data = tree2) +
  geom_tree(data = tree3) +
  #geom_tiplab(data = tree1) +
  #geom_tiplab(data = tree2, hjust = 0) +
  #geom_tiplab(data = tree3, hjust = 0) +
  geom_line(
    aes(x, y, group = label, colour = of_interest),
    data = jingmen_s,
    alpha = 0.5
  ) + 
  annotate("text", x = c(0, 0.2, 0.35), y = 520, label = c("NFLG", "ENV", 'NS5')) +
  scale_y_continuous(limits = c(0, 530))
################################### OUTPUT #####################################
# Save output files, plots, or results

#################################### END #######################################
################################################################################
