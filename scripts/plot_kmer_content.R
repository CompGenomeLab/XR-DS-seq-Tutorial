#!/usr/bin/env Rscript

#### Libraries ####
suppressPackageStartupMessages({
  library(ggplot2)
  library(argparser)
  library(tidyr)
  library(dplyr)
  library(logger)
})

set.seed(1)

######## Arguments ##########
p <- arg_parser("Produce nucleotide/dinucleotide content plots")
p <- add_argument(p, "-i", help = "input table")
p <- add_argument(p, "-o", help = "output image (png)")
p <- add_argument(p, "-l", help = "log file")
p <- add_argument(p, "-k", help = "kmer: '1' (nucleotide) or '2' (dinucleotide)")
p <- add_argument(p, "-s", help = "sample name in caption")
p <- add_argument(p, "-f", help = "filter (only for k=2). Ex: 'CC,CT,TC,TT'", nargs = 1, default = "")
argv <- parse_args(p)

# log file
log_appender(appender_file(argv$l))

#### Read & reshape ####
log_info("Reading file: {argv$i}")
nuc_table <- read.table(argv$i, header = TRUE, check.names = FALSE)

log_info("Renaming columns...")
x_order <- character()
if (argv$k == "1") {
  for (i in 2:ncol(nuc_table)) {
    colnames(nuc_table)[i] <- paste0(i - 1)
    x_order <- c(x_order, paste0(i - 1))
  }
  mytitle <- "Nucleotide Content of Oligomers"
} else if (argv$k == "2") {
  for (i in 2:ncol(nuc_table)) {
    colnames(nuc_table)[i] <- paste0(i - 1, "-", i)
    x_order <- c(x_order, paste0(i - 1, "-", i))
  }
  mytitle <- "Dinucleotide Content of Oligomers"
} else {
  stop("k must be '1' or '2'")
}

dt_organized <- nuc_table %>%
  gather(positions, counts, 2:ncol(nuc_table))
colnames(dt_organized) <- c("nucleotides", "positions", "counts")

# Per-position frequency (each position sums to 100%)
log_info("Computing per-position frequencies...")
dt_organized <- dt_organized %>%
  group_by(positions) %>%
  mutate(freq = 100 * counts / sum(counts)) %>%
  ungroup()

# Optional filtering for k=2 AFTER frequency calculation
if (argv$k == "2" && nzchar(argv$f)) {
  filt <- unlist(strsplit(argv$f, ","))
  log_info("Filtering dinucleotides (post-frequency): {paste(filt, collapse=',')}")
  dt_organized <- dt_organized %>% filter(nucleotides %in% filt)
}

# Order positions on x-axis
dt_organized$positions <- factor(dt_organized$positions, levels = x_order)

# For mononucleotides, use specified order T, C, G, A
if (argv$k == "1") {
  dt_organized$nucleotides <- factor(dt_organized$nucleotides,
                                     levels = c("A","G","C","T"))
}

#### Plot ####
log_info("Plotting...")
p <- ggplot(dt_organized, aes(x = positions, y = freq, fill = nucleotides)) +
  geom_col() +
  ylim(0, 100) +
  ylab("Frequency (%)") +
  xlab(ifelse(argv$k == "1", "Position in oligomers", "Dinucleotide position")) +
  labs(title = mytitle, subtitle = NULL, caption = argv$s) +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x  = element_text(size = 12, vjust = 0.6, angle = 65),
    axis.text.y  = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 12),
    plot.caption = element_text(size = 12, face = "italic"),
    plot.title   = element_text(size = 14, face = "bold")
  )

ggsave(argv$o, p, width = 10, height = 8, dpi = 150)
log_info("Saved to {argv$o}")


