# Introduction ------------------------------------------------------------
# The MultiQC plots are ok but not great for presentation. I'll make some 
# nicer plots here.

library(tidyverse)


# Importing and tidying the data ------------------------------------------
# Starting with the general statistics
general_stats <- read_table(file.path(getwd(), "data", "multiqc", "general_stats.tsv"),
                            col_names = TRUE)
general_stats$percent_dups <- as.numeric(str_remove(general_stats$percent_dups, "%"))
general_stats$percent_GC <- as.numeric(str_remove(general_stats$percent_GC, "%"))
general_stats$sample_name <- str_remove(general_stats$sample_name, "_001")

general_stats <- general_stats %>% 
  separate(sample_name, c("sample_name", "read_side"), sep = "_")
general_stats$sample_name <- str_sub(general_stats$sample_name, end=-4)
general_stats$sample_name <- str_remove(general_stats$sample_name, "Tamaulipas-pistils-")

general_stats$library_preparer <- rep(c("Us", "Us", "Genewiz", "Genewiz"), 5)
general_stats$library_preparer <- factor(general_stats$library_preparer, levels = c("Us", "Genewiz"))

# Now the mapping stats
mapping_stats <- read_table(file.path(getwd(), "data", "multiqc", "mapping_stats.tsv"),
                            col_names = TRUE)
mapping_stats$percent_aligned_salmon <- as.numeric(str_remove(mapping_stats$percent_aligned_salmon, "%"))
mapping_stats$percent_aligned_star <- as.numeric(str_remove(mapping_stats$percent_aligned_star, "%"))
mapping_stats$library_preparer <- rep(c("Us", "Genewiz"), 5)
mapping_stats$library_preparer <- factor(mapping_stats$library_preparer, levels = c("Us", "Genewiz"))
mapping_stats$sample_name <- str_sub(mapping_stats$sample_name, end=-4)
mapping_stats$sample_name <- str_remove(mapping_stats$sample_name, "Tamaulipas-pistils-")


# Making some plots -------------------------------------------------------
# First, a plot for the number of reads
# color_vec <- c("#2F69FF", "#DC267F")
color_vec <- c("#DC267F", "#2F69FF")

ggplot(general_stats, aes(x = sample_name, 
                          y = million_sequences,
                          fill = library_preparer)) +
  geom_bar(color = "black",
           position = "dodge",
           stat = "identity",
           size = 0.8) +
  geom_text(aes(label=million_sequences,
                fontface = "bold"), 
            position = position_dodge(width = 0.9),
            vjust = -0.5,
            size = 5) +
  scale_fill_manual(values = color_vec) +
  labs(title = "Number of reads per library",
       x = "Tamaulipas pistils",
       y = "Million reads") +
  scale_y_continuous(limits = c(0, 70),
                     expand = expansion(mult = c(0, .05))) +
  theme_bw() +
  theme(axis.title = element_text(size = 26, face = 'bold'),
        axis.text = element_text(size = 22, face = 'bold', color = 'black'),
        axis.text.x = element_text(size = 20, face = 'bold', color = 'black'),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
        plot.title = element_text(size = 32, face = 'bold', margin = margin(0, 0, 10, 0)),
        # axis.title.x = element_blank(),
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

ggsave(filename = './plots/number_of_reads.png',
       device = 'png',
       width = 9,
       height = 6,
       dpi = 400,
       units = 'in')

# Percent duplicates
percent_dups <- general_stats %>%
  select(sample_name, percent_dups, library_preparer) %>%
  group_by(sample_name, library_preparer) %>%
  summarize(percent_dups_mean = mean(percent_dups))
  
ggplot(percent_dups, aes(x = sample_name, 
                          y = percent_dups_mean,
                          fill = library_preparer)) +
  geom_bar(color = "black",
           position = "dodge",
           stat = "identity",
           size = 0.8) +
  geom_text(aes(label =format(round(percent_dups_mean, 1), nsmall = 1),
                fontface = "bold"), 
            position = position_dodge(width = 0.9),
            vjust = -0.5,
            size = 5) +
  scale_fill_manual(values = color_vec) +
  labs(title = "Percent duplicated sequences",
       x = "Tamaulipas pistils",
       y = "Percent") +
  scale_y_continuous(limits = c(0, 100),
                     expand = expansion(mult = c(0, .05))) +
  theme_bw() +
  theme(axis.title = element_text(size = 26, face = 'bold'),
        axis.text = element_text(size = 22, face = 'bold', color = 'black'),
        axis.text.x = element_text(size = 20, face = 'bold', color = 'black'),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
        plot.title = element_text(size = 32, face = 'bold', margin = margin(0, 0, 10, 0)),
        # axis.title.x = element_blank(),
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

ggsave(filename = './plots/percent_dups.png',
       device = 'png',
       width = 9,
       height = 6,
       dpi = 400,
       units = 'in')

# GC content
gc_content <- general_stats %>%
  select(sample_name, percent_GC, library_preparer) %>%
  group_by(sample_name, library_preparer) %>%
  summarize(gc_content_mean = mean(percent_GC))

ggplot(gc_content, aes(x = sample_name, 
                       y = gc_content_mean,
                       fill = library_preparer)) +
  geom_bar(color = "black",
           position = "dodge",
           stat = "identity",
           size = 0.8) +
  geom_text(aes(label =format(round(gc_content_mean, 1), nsmall = 1),
                fontface = "bold"), 
            position = position_dodge(width = 0.9),
            vjust = -0.5,
            size = 5) +
  scale_fill_manual(values = color_vec) +
  labs(title = "Percent GC content",
       x = "Tamaulipas pistils",
       y = "Percent") +
  scale_y_continuous(limits = c(0, 100),
                     expand = expansion(mult = c(0, .05))) +
  theme_bw() +
  theme(axis.title = element_text(size = 26, face = 'bold'),
        axis.text = element_text(size = 22, face = 'bold', color = 'black'),
        axis.text.x = element_text(size = 20, face = 'bold', color = 'black'),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
        plot.title = element_text(size = 32, face = 'bold', margin = margin(0, 0, 10, 0)),
        # axis.title.x = element_blank(),
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

ggsave(filename = './plots/gc_content.png',
       device = 'png',
       width = 9,
       height = 6,
       dpi = 400,
       units = 'in')

# STAR percent mapped
ggplot(mapping_stats, aes(x = sample_name, 
                       y = percent_aligned_star,
                       fill = library_preparer)) +
  geom_bar(color = "black",
           position = "dodge",
           stat = "identity",
           size = 0.8) +
  geom_text(aes(label =format(round(percent_aligned_star, 1), nsmall = 1),
                fontface = "bold"), 
            position = position_dodge(width = 0.9),
            vjust = -0.5,
            size = 5) +
  scale_fill_manual(values = color_vec) +
  labs(title = "STAR mapped reads",
       x = "Tamaulipas pistils",
       y = "Percent") +
  scale_y_continuous(limits = c(0, 100),
                     expand = expansion(mult = c(0, .05))) +
  theme_bw() +
  theme(axis.title = element_text(size = 26, face = 'bold'),
        axis.text = element_text(size = 22, face = 'bold', color = 'black'),
        axis.text.x = element_text(size = 20, face = 'bold', color = 'black'),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
        plot.title = element_text(size = 32, face = 'bold', margin = margin(0, 0, 10, 0)),
        # axis.title.x = element_blank(),
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

ggsave(filename = './plots/star_percent_mapped.png',
       device = 'png',
       width = 9,
       height = 6,
       dpi = 400,
       units = 'in')

# Salmon percent mapped
ggplot(mapping_stats, aes(x = sample_name, 
                          y = percent_aligned_salmon,
                          fill = library_preparer)) +
  geom_bar(color = "black",
           position = "dodge",
           stat = "identity",
           size = 0.8) +
  geom_text(aes(label =format(round(percent_aligned_salmon, 1), nsmall = 1),
                fontface = "bold"), 
            position = position_dodge(width = 0.9),
            vjust = -0.5,
            size = 5) +
  scale_fill_manual(values = color_vec) +
  labs(title = "Salmon mapped reads",
       x = "Tamaulipas pistils",
       y = "Percent") +
  scale_y_continuous(limits = c(0, 100),
                     expand = expansion(mult = c(0, .05))) +
  theme_bw() +
  theme(axis.title = element_text(size = 26, face = 'bold'),
        axis.text = element_text(size = 22, face = 'bold', color = 'black'),
        axis.text.x = element_text(size = 20, face = 'bold', color = 'black'),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
        plot.title = element_text(size = 32, face = 'bold', margin = margin(0, 0, 10, 0)),
        # axis.title.x = element_blank(),
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

ggsave(filename = './plots/salmon_percent_mapped.png',
       device = 'png',
       width = 9,
       height = 6,
       dpi = 400,
       units = 'in')
