# Introduction ------------------------------------------------------------
# Attempting to replicate my collaborators MDS plot using edgeR.

library(edgeR)
library(tximport)
library(tximeta)
library(splicejam)

# Importing the data ------------------------------------------------------
salmon_counts <- read.table(file = file.path(getwd(),
                                             "data",
                                             "summarized_salmon_counts.txt"),
                            sep = '\t',
                            header = TRUE)

sample_groups <- rep(c("R1", "R2"), 5)
names(sample_groups) <- sub("_$", "", colnames(salmon_counts)) 

salmon_DGEList = DGEList(counts = salmon_counts,
                         group = sample_groups,
                         remove.zeros = TRUE)

mds <- plotMDS(salmon_DGEList)

mds <- mds@.Data[[3]] %>%
  as.data.frame() %>%
  # colnames(c("Dim1", "Dim2")) %>%
  rownames_to_column("SampleID")

mds$SampleID <- sub("_$", "", colnames(salmon_counts)) 

mds <- mds %>%
  separate(SampleID, c("accession", "tissue", "time", "temp", "prep"), "\\.") %>%
  unite("short_name", c("time", "temp"), sep = '-')

ggplot(mds, aes(V1, V2, color = short_name, shape = prep)) +
  geom_point(size = 6) +
  labs(title = "edgeR MDS",
       x = paste0("Leading logFC dim1"),
       y = paste0("Leading logFC dim2")) +
  # scale_color_manual(values = c("#2F69FF", "#DC267F", "#FFB000")) +
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
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'cm'),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 20, face = 'bold'),
        legend.position = 'right',
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 26, face = 'bold'))

ggsave(filename = './plots/MDS.png',
       device = 'png',
       width = 9,
       height = 6,
       dpi = 400,
       units = 'in')

# Flipping it so it matches the PCA
ggplot(mds, aes(V1, V2, color = short_name, shape = prep)) +
  geom_point(size = 6) +
  labs(title = "edgeR MDS (axes reversed)",
       x = paste0("Leading logFC dim1"),
       y = paste0("Leading logFC dim2")) +
  # scale_color_manual(values = c("#2F69FF", "#DC267F", "#FFB000")) +
  scale_y_reverse() +
  scale_x_reverse() +
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
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'cm'),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 20, face = 'bold'),
        legend.position = 'right',
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 26, face = 'bold'))

ggsave(filename = './plots/reverse_MDS.png',
       device = 'png',
       width = 9,
       height = 6,
       dpi = 400,
       units = 'in')






