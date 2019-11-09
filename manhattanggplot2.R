#!/bin/bash
# questo tutorial l'ho preso da qui http://www.r-bloggers.com/manhattan-plot-with-simulated-data-using-ggplot2/
exec R --vanilla -q --slave -e "source(file=pipe(\"tail -n +4 $0\"))" --args $@

#debug: exec R --vanilla --verbose -e "source(file=pipe(\"tail -n +4 $0\"))" --args $@
### The above line starts R and then reads in this script, starting at line 4
### (taken from the R Wiki at http://rwiki.sciviews.org/doku.php?id=tips:scriptingr).

# simulated_marker_effects.r
#
# Use ggplot2 to produce Manhattan plots of simulated marker effects.
#
# Author:       John B. Cole (john.cole@ars.usda.gov)
# Changes:      07/20/2011      Original program

library("ggplot2")
library("RColorBrewer")

# Create 30 chromosomes with 1,000 SNP per chromosome. This assumes that all chromosomes are the same
# length and have the same number of SNP.
chromosomes <- data.frame(chrome = rep(seq(1:30),each=1000), seq_chrome = rep(seq(1:1000), each=30),
                          snp_id = seq(1:30000))

# Simulate random SNP locations on each chromosome, and add an offset for each chromosome
# so that the SNP location is relative to the beginning of the genome, not an individual
# chromosome.
print("Adding snp_locs to chromosomes dataframe")
for ( c in 1:30 ) {
  if ( c == 1 ) { snp_loc <- sort(as.integer(runif(1000, 1, 100000001))) }
  else {
    new_snp_loc <- sort(as.integer(runif(1000, 1, 100000001))) + ( (c-1) * 100000000 )
    snp_loc <- c(snp_loc, new_snp_loc)
  }
}
# Sort the SNP locations so that the SNP on each chromosome are sequential.
sort(snp_loc)
chromosomes$snp_loc <- snp_loc
print("Structure of the chromosomes dataframe")
str(chromosomes)

# Simulate an additive genetic effect for each SNP from a t-distribution.
print("Simulating marker effects")
solutions <- data.frame(snp_id = seq(1:30000), effect = abs(rt(30000,df=100)))

# Merge the SNP effects into the dataframe with the other marker information.
print("Merging solutions with chromosomes")
yld <- merge(solutions, chromosomes, by="snp_id")
yld$chrome <- as.factor(yld$chrome)
yld[with(yld, order(chrome, snp_id)),]

#agg <- data.frame(seq_chrome=yld$seq_chrome, chrome=yld$chrome)
print("Finding chromosome starts and ends")
agg <- data.frame(snp_loc=yld$snp_loc, chrome=yld$chrome)

# Several things happen in this section of code:
agg$chrome <- as.character(agg$chrome)
# 1. Find the location of the first and last SNP on each chromosome
chrome_starts <- aggregate(agg, by = list(agg$chrome), FUN = min, na.rm = T)
chrome_ends <- aggregate(agg, by = list(agg$chrome), FUN = max, na.rm = T)
chrome_starts$chrome <- as.numeric(chrome_starts$chrome)
chrome_starts <- chrome_starts[order(chrome_starts$chrome),]
chrome_ends$chrome <- as.numeric(chrome_ends$chrome)
chrome_ends <- chrome_ends[order(chrome_ends$chrome),]
chrome_starts <- chrome_starts[c(2,3)]
chrome_ends <- chrome_ends[c(2,3)]
names(chrome_starts) <- c("chrome_start","chrome")
names(chrome_ends) <- c("chrome_end","chrome")
# 2. Merge the starts and ends into s single dataframe. These data are then used to construct
#    a list of the locations on the x-axis of each tick-mark. This is done by calculating the
#    midpoint between the first and last marker on each chromosome.
chrome_start_end <- merge(chrome_starts, chrome_ends, by="chrome")
chrome_start_end$length <- chrome_start_end$chrome_end - chrome_start_end$chrome_start
chrome_start_end$tick <- chrome_start_end$chrome_start + floor( chrome_start_end$length / 2  )
print("Chromosome starts and ends")
print(chrome_start_end)
# This is the list of ticks that will be merged into the data frame with the SNP information.
ticks <- chrome_start_end[c("chrome", "tick")]

# Merge ticks back into dataset
print("Merging SNP effects with chromosome starts and ends")
yld <- merge(yld, ticks, by="chrome")

# Change the chromosome 30 label to "X"
yld$label <- as.character(yld$chrome)
yld$label[yld$chrome == 30] <- "X"

# Construct our color palette by taking the first 7 colors from the Brewer palette.
colours <- rep(c(brewer.pal(n = 7, name = "Set1")),5)

# Build the actual plot
print ("Generating the plot")
manhattan <- ggplot(data=yld, aesthetics=aes(snp_loc, effect, color=chrome)) +
  geom_point() +
  scale_colour_manual(values = colours) +
  theme(legend.position = "none") +
  scale_x_continuous(name = "Chromosome", breaks = unique(yld$tick), labels = unique(yld$label)) +
  ylab("Marker Effects") +
  #opts(title = "Distribution of simulated marker effects") #deprecated
  ggtitle("Distribution of simulated marker effects") #nuovo title
ggsave("simulated_manhattan_plot.png", width = 10, height = 5, dpi = 300)