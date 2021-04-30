library("Analysis5204")
library("ggplot2")

# Create a Results directory in the top directory
dir.create("Results/")

olink_dist(plasma_npx, "Olink CARDIOVASCULAR III", "Results/test.png")
