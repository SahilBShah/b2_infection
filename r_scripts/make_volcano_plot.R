library(ggplot2)
library(cowplot)
library(dplyr)

setwd("/Users/sahilshah/Documents/b2_infection/data/DE_analysis")

subcategory <- read.delim("mm10.aligned.sorted_classification.txt", sep = "\t")
subcategory <- subcategory[c("isoform", "subcategory")]
colnames(subcategory) <- c("Transcript", "subcategory")

fragments <- c("3prime_fragment", "5prime_fragment", "internal_fragment", "mono-exon", "mono-exon_by_intron_retention")

brf1_48_mock <- read.delim("glaunsinger_Brf1_48hpi_v_Brf1_Mock_DE_analysis_mm10.txt", sep = "\t")
brf1_48_mock <- left_join(x=brf1_48_mock, y=subcategory, by="Transcript")
brf1_48_mock <- brf1_48_mock %>% filter(!subcategory %in% fragments)
brf1_ctrl_inf <- read.delim("glaunsinger_Brf1_48hpi_v_Ctrl_48_DE_analysis_mm10.txt", sep = "\t")
brf1_ctrl_inf <- left_join(x=brf1_ctrl_inf, y=subcategory, by="Transcript")
brf1_ctrl_inf <- brf1_ctrl_inf %>% filter(!subcategory %in% fragments)
brf1_ctrl_mock <- read.delim("glaunsinger_Brf1_Mock_v_Ctrl_Mock_DE_analysis_mm10.txt", sep = "\t")
brf1_ctrl_mock <- left_join(x=brf1_ctrl_mock, y=subcategory, by="Transcript")
brf1_ctrl_mock <- brf1_ctrl_mock %>% filter(!subcategory %in% fragments)
ctrl_48_mock <- read.delim("glaunsinger_Ctrl_48_v_Ctrl_Mock_DE_analysis_mm10.txt", sep = "\t")
ctrl_48_mock <- left_join(x=ctrl_48_mock, y=subcategory, by="Transcript")
ctrl_48_mock <- ctrl_48_mock %>% filter(!subcategory %in% fragments)

ggsave <- function(..., bg = 'white') ggplot2::ggsave(..., bg = bg)

#Volcano plot of Brf1 conditions
ggplot(brf1_48_mock, aes(x=logFC, y=-log(adj.P.Val,10))) + # -log10 conversion  
  geom_point(size = 2/5, aes(color=-log(adj.P.Val,10)>=1.3 & (logFC>=1.5 | logFC<=-1.5)), show.legend=FALSE) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  ggtitle("Brf1 48hpi vs. Brf1 Mock") +
  geom_hline(yintercept=1.3, linetype="dashed", color="red") +
  geom_vline(xintercept=1.5, linetype="dashed", color="red") +
  geom_vline(xintercept=-1.5, linetype="dashed", color="red") +
  scale_color_manual(values=c("black", "blue")) +
  theme_cowplot()

ggsave("../../figures/brf1_48_mock_DE.png")

#Volcano plot of Ctrl conditions
ggplot(brf1_ctrl_inf, aes(x=logFC, y=-log(adj.P.Val,10))) + # -log10 conversion  
  geom_point(size = 2/5, aes(color=-log(adj.P.Val,10)>=1.3 & (logFC>=1.5 | logFC<=-1.5)), show.legend=FALSE) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  ggtitle("Brf1 48hpi vs. Ctrl 48hpi") +
  geom_hline(yintercept=1.3, linetype="dashed", color="red") +
  geom_vline(xintercept=1.5, linetype="dashed", color="red") +
  geom_vline(xintercept=-1.5, linetype="dashed", color="red") +
  scale_color_manual(values=c("black", "blue")) +
  theme_cowplot()

ggsave("../../figures/brf1_ctrl_48hpi.png")

ggplot(brf1_ctrl_mock, aes(x=logFC, y=-log(adj.P.Val,10))) + # -log10 conversion  
  geom_point(size = 2/5, aes(color=-log(adj.P.Val,10)>=1.3 & (logFC>=1.5 | logFC<=-1.5)), show.legend=FALSE) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  ggtitle("Brf1 Mock vs. Ctrl Mock") +
  geom_hline(yintercept=1.3, linetype="dashed", color="red") +
  geom_vline(xintercept=1.5, linetype="dashed", color="red") +
  geom_vline(xintercept=-1.5, linetype="dashed", color="red") +
  scale_color_manual(values=c("black", "blue")) +
  theme_cowplot()

ggsave("../../figures/brf1_ctrl_mock_DE.png")

ggplot(ctrl_48_mock, aes(x=logFC, y=-log(adj.P.Val,10))) + # -log10 conversion  
  geom_point(size = 2/5, aes(color=-log(adj.P.Val,10)>=1.3 & (logFC>=1.5 | logFC<=-1.5)), show.legend=FALSE) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  ggtitle("Ctrl 48hpi vs. Ctrl Mock") +
  geom_hline(yintercept=1.3, linetype="dashed", color="red") +
  geom_vline(xintercept=1.5, linetype="dashed", color="red") +
  geom_vline(xintercept=-1.5, linetype="dashed", color="red") +
  scale_color_manual(values=c("black", "blue")) +
  theme_cowplot()

ggsave("../../figures/ctrl_48_mock_DE.png")

#MHV68
mhv68 <- read.delim("glaunsinger_Brf1_48hpi_v_Ctrl_48_DE_analysis_mhv68.txt", sep = "\t")

ggplot(mhv68, aes(x=logFC, y=-log(adj.P.Val,10))) + # -log10 conversion  
  geom_point(size = 2/5, aes(color=-log(adj.P.Val,10)>=1.3 & (logFC>=1.5 | logFC<=-1.5)), show.legend=FALSE) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  ggtitle("Ctrl 48hpi vs. Ctrl Mock") +
  geom_hline(yintercept=1.3, linetype="dashed", color="red") +
  geom_vline(xintercept=1.5, linetype="dashed", color="red") +
  geom_vline(xintercept=-1.5, linetype="dashed", color="red") +
  scale_color_manual(values=c("black", "blue")) +
  theme_cowplot()

ggsave("../../figures/mhv68_DE.png")

#Gene level
brf1_48_mock <- read.delim("glaunsinger_Brf1_48hpi_v_Brf1_Mock_DE_analysis_mm10_genes.txt", sep = "\t")
brf1_48_mock <- left_join(x=brf1_48_mock, y=subcategory, by="Transcript")
brf1_48_mock <- brf1_48_mock %>% filter(!subcategory %in% fragments)
brf1_ctrl_inf <- read.delim("glaunsinger_Brf1_48hpi_v_Ctrl_48_DE_analysis_mm10_genes.txt", sep = "\t")
brf1_ctrl_inf <- left_join(x=brf1_ctrl_inf, y=subcategory, by="Transcript")
brf1_ctrl_inf <- brf1_ctrl_inf %>% filter(!subcategory %in% fragments)
brf1_ctrl_mock <- read.delim("glaunsinger_Brf1_Mock_v_Ctrl_Mock_DE_analysis_mm10_genes.txt", sep = "\t")
brf1_ctrl_mock <- left_join(x=brf1_ctrl_mock, y=subcategory, by="Transcript")
brf1_ctrl_mock <- brf1_ctrl_mock %>% filter(!subcategory %in% fragments)
ctrl_48_mock <- read.delim("glaunsinger_Ctrl_48_v_Ctrl_Mock_DE_analysis_mm10_genes.txt", sep = "\t")
ctrl_48_mock <- left_join(x=ctrl_48_mock, y=subcategory, by="Transcript")
ctrl_48_mock <- ctrl_48_mock %>% filter(!subcategory %in% fragments)

#Volcano plot of Brf1 conditions
ggplot(brf1_48_mock, aes(x=logFC, y=-log(adj.P.Val,10))) + # -log10 conversion  
  geom_point(size = 2/5, aes(color=-log(adj.P.Val,10)>=1.3 & (logFC>=1.5 | logFC<=-1.5)), show.legend=FALSE) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  ggtitle("Brf1 48hpi vs. Brf1 Mock") +
  geom_hline(yintercept=1.3, linetype="dashed", color="red") +
  geom_vline(xintercept=1.5, linetype="dashed", color="red") +
  geom_vline(xintercept=-1.5, linetype="dashed", color="red") +
  scale_color_manual(values=c("black", "blue")) +
  theme_cowplot()

ggsave("../../figures/brf1_48_mock_DE_genes.png")

#Volcano plot of Ctrl conditions
ggplot(brf1_ctrl_inf, aes(x=logFC, y=-log(adj.P.Val,10))) + # -log10 conversion  
  geom_point(size = 2/5, aes(color=-log(adj.P.Val,10)>=1.3 & (logFC>=1.5 | logFC<=-1.5)), show.legend=FALSE) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  ggtitle("Brf1 48hpi vs. Ctrl 48hpi") +
  geom_hline(yintercept=1.3, linetype="dashed", color="red") +
  geom_vline(xintercept=1.5, linetype="dashed", color="red") +
  geom_vline(xintercept=-1.5, linetype="dashed", color="red") +
  scale_color_manual(values=c("black", "blue")) +
  theme_cowplot()

ggsave("../../figures/brf1_ctrl_48hpi_genes.png")

ggplot(brf1_ctrl_mock, aes(x=logFC, y=-log(adj.P.Val,10))) + # -log10 conversion  
  geom_point(size = 2/5, aes(color=-log(adj.P.Val,10)>=1.3 & (logFC>=1.5 | logFC<=-1.5)), show.legend=FALSE) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  ggtitle("Brf1 Mock vs. Ctrl Mock") +
  geom_hline(yintercept=1.3, linetype="dashed", color="red") +
  geom_vline(xintercept=1.5, linetype="dashed", color="red") +
  geom_vline(xintercept=-1.5, linetype="dashed", color="red") +
  scale_color_manual(values=c("black", "blue")) +
  theme_cowplot()

ggsave("../../figures/brf1_ctrl_mock_DE_genes.png")

ggplot(ctrl_48_mock, aes(x=logFC, y=-log(adj.P.Val,10))) + # -log10 conversion  
  geom_point(size = 2/5, aes(color=-log(adj.P.Val,10)>=1.3 & (logFC>=1.5 | logFC<=-1.5)), show.legend=FALSE) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  ggtitle("Ctrl 48hpi vs. Ctrl Mock") +
  geom_hline(yintercept=1.3, linetype="dashed", color="red") +
  geom_vline(xintercept=1.5, linetype="dashed", color="red") +
  geom_vline(xintercept=-1.5, linetype="dashed", color="red") +
  scale_color_manual(values=c("black", "blue")) +
  theme_cowplot()

ggsave("../../figures/ctrl_48_mock_DE_genes.png")
