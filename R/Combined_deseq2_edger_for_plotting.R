# Introduction ------------------------------------------------------------
# In this script I will do some dimensionality reduction with Salmon output
# and look at differential gene expression using DESeq2.

# library(tidyverse)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tximport)
library(edgeR)
library(tximeta)
library(splicejam)

# Importing the data ------------------------------------------------------
salmon_counts <- read.table(file = file.path(getwd(), 
                                             "data",
                                             "salmon_output.txt"),
                            sep = '\t',
                            header = TRUE)


# Making tx2gene file -----------------------------------------------------
# tximport needs a tx2gene file to summarize transcript abundance to gene
# abundance. I'll make one here using a function from the splicejam package.
library(splicejam)
tx2gene <- makeTx2geneFromGtf(file.path(getwd(), "data", "ITAG4.0_gene_models.gtf"))
tx2gene$transcript_id <- sub("mRNA:", "", tx2gene$transcript_id)
tx2gene$gene_id <- sub("gene:", "", tx2gene$gene_id)


# Importing the data with tximport ----------------------------------------
salmon_directory <- file.path(getwd(), "data", "salmon_output")
salmon_file_list <- list.files(salmon_directory)

samples <- data.frame(row.names = salmon_file_list,
                      library_prep = factor(rep(c("R1", "R2"), 5)),
                      time = c("0hr", "0hr", 
                               rep("3hr", 4),
                               rep("8hr", 4)),
                      temp = c(rep("25C", 4),
                               rep("37C", 2),
                               rep("25C", 2),
                               rep("37C", 2)),
                      quant_file_path = file.path(salmon_directory, salmon_file_list, "quant.sf"))

txi <- tximport(samples$quant_file_path, type = "salmon", tx2gene = tx2gene)

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ library_prep)

dds <- DESeq(ddsTxi)

vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c("library_prep"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData <- pcaData %>%
  separate(name, c("accession", "tissue", "time", "temp", "prep"), "-") %>%
  unite("short_name", c("time", "temp"), sep = '-')

deseq2_plot <- ggplot(pcaData, aes(PC1, PC2, color = short_name, shape = library_prep)) +
  geom_point(size = 6) +
  labs(title = "DESeq2 PCA",
       x = paste0("PC1: ",percentVar[1],"% variance"),
       y = paste0("PC2: ",percentVar[2],"% variance")) +
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




# edgeR -------------------------------------------------------------------
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

edger_plot <- ggplot(mds, aes(V1, V2, color = short_name, shape = prep)) +
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



# Flipping it so it matches the PCA
edger_plot_reverse <- ggplot(mds, aes(V1, V2, color = short_name, shape = prep)) +
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



# Align and save ----------------------------------------------------------
library(cowplot)
plots_aligned <- align_plots(deseq2_plot,
                             edger_plot,
                             edger_plot_reverse,
                             align = "v")

ggsave(filename = './plots/PCA.png',
       plot = plots_aligned[[1]],
       device = 'png',
       width = 9,
       height = 6,
       dpi = 400,
       units = 'in')

ggsave(filename = './plots/MDS.png',
       plot = plots_aligned[[2]],
       device = 'png',
       width = 9,
       height = 6,
       dpi = 400,
       units = 'in')

ggsave(filename = './plots/reverse_MDS.png',
       plot = plots_aligned[[3]],
       device = 'png',
       width = 9,
       height = 6,
       dpi = 400,
       units = 'in')




