####################
##### Plotting #####
####################

library(ggplot2)
library(dplyr)
library(ggrepel)


##### Volcano Plot #####
results_FTDvControl = read.csv("FTD_vs_Control_DCTP.csv")
rownames(results_FTDvControl) = results_FTDvControl$X
direction_colors <- list(up_sig="#A63B00", up="#EC530080", down_sig="#477397", down="#6497BF80")

# Define significant genes
results_FTDvControl$significant <- (abs(results_FTDvControl$logFC) > 0.01) & 
  (results_FTDvControl$adj.P.Val < 0.05)
results_FTDvControl$gene <- rownames(results_FTDvControl)  # Add gene column


ggplot(results_FTDvControl, aes(x=logFC, y=-log10(adj.P.Val))) +
  # Plot all points in grey first
  geom_point(size=4, color='grey', alpha=0.5) +
  # Overlay significant upregulated genes
  geom_point(data=results_FTDvControl %>% filter(logFC > 0.02),
             size=4,
             color=direction_colors[['up']]) +
  geom_point(data=results_FTDvControl %>% filter(significant & logFC > 0.01),
             size=4,
             color=direction_colors[['up_sig']]) +
  # Overlay significant downregulated genes
  geom_point(data=results_FTDvControl %>% filter(logFC < -0.02),
             size=4,
             color=direction_colors[['down']]) +
  geom_point(data=results_FTDvControl %>% filter(significant & logFC < -0.02),
             size=4,
             color=direction_colors[['down_sig']]) +
  # Add circle around genes with padj < 0.05
  geom_point(data=results_FTDvControl %>% filter(adj.P.Val < 0.05), 
             col="black", pch=1, size=5) +
  # Add threshold lines
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=0.01, linetype="dashed") +
  geom_vline(xintercept=-0.01, linetype="dashed") +
  # Add labels for genes with padj < 0.05
  ggrepel::geom_label_repel(data=results_FTDvControl %>% filter(adj.P.Val < 0.05),
                            aes(label=gene),
                            max.overlaps=Inf,
                            fontface='italic',
                            size=6,
                            segment.color='black',
                            segment.size=0.5,
                            box.padding=0.5,
                            point.padding=0.3) +
  #xlim(-0.15, 0.15) +
  #ylim(-0, 3) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text = element_text(size=20),
        axis.title = element_text(size=20)) +
  labs(x="Log2 Fold Change", 
       y="-Log10(padj)",
       title="FTD vs. Control")


##### Stacked barplot #####
proportions = read.csv("Proportions.csv")
#proportions$study_disease = factor(proportions$study_disease, levels = c("Control-Biogen", "FTD-Biogen", "Control-Marsan", "FTD-Marsan", "Control-Mathys", "AD-Mathys"))
proportions$clusters = factor(proportions$clusters, levels = c("Neuron", "Lymphocyte", "Myeloid", "Astrocyte", "OL", "OPC", "BVC", "CP"))
proportions$clusters

proportions_pct_disease <- proportions %>%
  group_by(condition, clusters) %>%
  summarise(Proportions = sum(Proportions), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(
    percentage = Proportions / sum(Proportions) * 100,
    label = ifelse(percentage > 3, sprintf("%.1f%%", percentage), "")  # Only show if >3%
  )


celltype_color = c(
  "Astrocyte"= "#abd158", 
  "BVC"= "#74858c",  
  "CP"= "#d5c18e",   
  "Neuron"= "#21409a",
  "Lymphocyte"= "#68151f",
  "Myeloid"= "#d24a75",
  "OL"= "#fdb515",
  "OPC"= "#ed7020"
)

ggplot(proportions_pct_disease, aes(x=condition, y=Proportions, fill=clusters)) +
  geom_bar(position="fill", stat="identity") +
  geom_text(aes(label=label), 
            position=position_fill(vjust=0.5), 
            size=3, 
            color="white") +
  theme_classic() +
  labs(y="Proportion") +
  scale_fill_manual(values = celltype_color) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))

## Transformed proportions
proportions_pct_disease <- proportions %>%
  group_by(condition, clusters) %>%
  summarise(TransformedProportions = sum(TransformedProportions), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(
    percentage = TransformedProportions / sum(TransformedProportions) * 100,
    label = ifelse(percentage > 5, sprintf("%.0f%%", percentage), "")  # Only show if >3%
  )

ggplot(proportions_pct_disease, aes(x=condition, y=TransformedProportions, fill=clusters)) +
  geom_bar(position="fill", stat="identity", width = 0.7) +
  geom_text(aes(label=label), 
            position=position_fill(vjust=0.5), 
            size=3, 
            color="white") +
  theme_classic() +
  labs(y="Transformed Proportion") +
  scale_fill_manual(values = celltype_color) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))





















