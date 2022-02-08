# Introduction ------------------------------------------------------------
# Making a metaplot of read coverage accross all samples. Following the 
# blog post located here: 
# https://groverj3.github.io/articles/2019-06-28_making-better-metaplots-with-ggplot-part-2.html

library(tidyverse)

# Importing the data ------------------------------------------------------
# Created with deepTools bamCoverage, computeMatrix, and plotProfile.

read_deeptools_table <- function(file) {
  n <- max(count.fields(file, sep = '\t'), na.rm = TRUE)
  x <- readLines(file)
  .splitvar <- function(x, sep, n) {
    var <- unlist(strsplit(x, split = sep))
    length(var) <- n
    return(var)
  }
  x <- do.call(cbind, lapply(x, .splitvar, sep = '\t', n = n))
  x <- apply(x, 1, paste, collapse = '\t')
  plot_table <- na.omit(read.csv(text = x, sep = '\t')[-1,])  # Remove first row with "gene" label
  return(plot_table)
}

table_input <- read_deeptools_table(file.path(getwd(), "data", "metaplot.tab"))

table_input <- gather(table_input, 'sample', 'score', -bin.labels, -bins)

# Adding factor for us vs genewiz
table_input$usvgenewiz <- NA
table_input$usvgenewiz[grep("R1", table_input$sample)] <- "Us"
table_input$usvgenewiz[grep("R2", table_input$sample)] <- "Genewiz"


# Plotting ----------------------------------------------------------------
plot_colors <- c("#009292", "#d1214d")
metaplot <- ggplot(table_input, aes(x = bins, y = as.numeric(score), color = sample)) +
  geom_vline(xintercept = c(100, 235), 
             linetype = 'dotted',
             size = 1) +
  geom_line(size = 1) +
  scale_color_manual(values = c(rep(plot_colors, 5))) +
  scale_x_continuous(breaks = table_input$bins,
                     labels = table_input$bin.labels) +
  labs(title = "Library coverage across all transcripts",
       x = "Relative transcript position",
       y = "Normalized coverage") +
  theme_bw() +
  theme(axis.title = element_text(size = 26, face = 'bold'),
        axis.text = element_text(size = 22, face = 'bold', color = 'black'),
        axis.text.x = element_text(size = 20, face = 'bold', color = 'black'),
        plot.title = element_text(size = 28, face = 'bold', margin = margin(0, 0, 10, 0)),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
        panel.border = element_blank(),
        axis.line = element_line(size = 1, color = 'black'),
        axis.ticks = element_line(size = 1, color = 'black'),
        axis.ticks.length = unit(8, 'pt'),
        axis.ticks.x = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'cm'),
        panel.grid = element_blank(),
        legend.position = 'none')

ggsave(filename = './plots/metaplot_together.png',
       plot = metaplot,
       device = 'png',
       width = 12,
       height = 8,
       dpi = 400,
       units = 'in') 

# Reps separated
# Making new column without the us v genewiz info, for the facet
table_input$rep_info <- table_input$sample
table_input$rep_info <- substr(table_input$rep_info, 1, nchar(table_input$rep_info) - 3)

table_input$usvgenewiz <- factor(table_input$usvgenewiz, levels = c("Us", "Genewiz"))

metaplot_separate <- ggplot(table_input, aes(x = bins, y = as.numeric(score), color = sample)) +
  geom_vline(xintercept = c(100, 235), 
             linetype = 'dotted',
             size = 1) +
  geom_line(size = 1) +
  scale_color_manual(values = c(rep(plot_colors, 5))) +
  scale_x_continuous(breaks = table_input$bins,
                     labels = table_input$bin.labels) +
  labs(title = "Library coverage across all transcripts",
       x = "Relative transcript position",
       y = "Normalized coverage") +
  facet_grid(vars(usvgenewiz), vars(rep_info)) +
  # facet_wrap(~sample, nrow = 2, ncol = 5) +
  theme_bw() +
  theme(axis.title = element_text(size = 26, face = 'bold'),
        axis.text = element_text(size = 22, face = 'bold', color = 'black'),
        axis.text.x = element_text(size = 10, face = 'bold', color = 'black'),
        plot.title = element_text(size = 28, face = 'bold', margin = margin(0, 0, 10, 0)),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
        panel.border = element_blank(),
        axis.line = element_line(size = 1, color = 'black'),
        axis.ticks = element_line(size = 1, color = 'black'),
        axis.ticks.length = unit(8, 'pt'),
        axis.ticks.x = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'cm'),
        panel.grid = element_blank(),
        legend.position = 'none',
        panel.spacing = unit(2, "lines"),
        strip.background = element_rect(color = "white", fill = "white"),
        strip.text.x = element_text(size = 12, face = 'bold', color = "black"),
        strip.text.y = element_text(size = 16, face = 'bold', color = "black"))

ggsave(filename = './plots/metaplot_seperate.png',
       plot = metaplot_separate,
       device = 'png',
       width = 16,
       height = 8,
       dpi = 400,
       units = 'in') 





