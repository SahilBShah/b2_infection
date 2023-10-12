library(ggplot2)
library(cowplot)
library(dplyr)
library(stringr)

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


#Volcano plot of Ctrl conditions colored based on B2 location
isoform_b2_loc <- read.csv("intact_b2_mapped_to_isoforms.csv")
colnames(isoform_b2_loc) <- c("Index", "Transcript", "B2_location")
isoform_b2_loc <- subset(isoform_b2_loc, select = -c(Index))
isoform_b2_loc <- isoform_b2_loc %>% distinct()

##Brf1 48hpi vs. Ctrl 48hpi
brf1_ctrl_inf_b2 <- left_join(brf1_ctrl_inf, isoform_b2_loc, by=join_by(Transcript == Transcript))
brf1_ctrl_inf_b2[is.na(brf1_ctrl_inf_b2)] <- "None/Undetermined"
ggplot(brf1_ctrl_inf_b2, aes(x=logFC, y=-log(adj.P.Val,10), label=B2_location)) + # -log10 conversion  
  geom_point(size = 2/5, aes(color=B2_location), show.legend=TRUE) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  ggtitle("Brf1 48hpi vs. Ctrl 48hpi") +
  geom_hline(yintercept=1.3, linetype="dashed", color="red") +
  geom_vline(xintercept=1.5, linetype="dashed", color="red") +
  geom_vline(xintercept=-1.5, linetype="dashed", color="red") +
  scale_color_manual(values=c("black", "blue", "green", "orange")) +
  theme_cowplot()

ggsave("../../figures/brf1_ctrl_inf_b2_mapped_to_DE_isoforms.png")

#Intact B2s mapped to all isoforms
all_iso_introns <- read.delim("../bed_files/All_isoforms_intactB2_mapped_to_introns.bed", sep = "\t", header=FALSE)
all_iso_introns <- subset(all_iso_introns, select = c(V4))
all_iso_introns <- as.data.frame(str_split_fixed(all_iso_introns$V4, '_', 2))
all_iso_introns <- subset(all_iso_introns, select = c(V1))
all_iso_introns <- all_iso_introns %>% distinct()
colnames(all_iso_introns) <- c("Transcript")

all_iso_exons <- read.delim("../bed_files/All_isoforms_intactB2_mapped_to_exons.bed", sep = "\t", header=FALSE)
all_iso_exons <- as.data.frame(all_iso_exons$V4)
colnames(all_iso_exons) <- c("Transcript")

all_iso_5prime <- read.delim("../bed_files/All_isoforms_intactB2_mapped_to_5prime_UTR.bed", sep = "\t", header=FALSE)
all_iso_5prime <- as.data.frame(all_iso_5prime$V4)
colnames(all_iso_5prime) <- c("Transcript")

all_iso_3prime <- read.delim("../bed_files/All_isoforms_intactB2_mapped_to_3prime_UTR.bed", sep = "\t", header=FALSE)
all_iso_3prime <- as.data.frame(all_iso_3prime$V4)
colnames(all_iso_3prime) <- c("Transcript")

all_iso_prom <- read.delim("../bed_files/All_isoforms_intactB2_mapped_to_promoter_up4kb.bed", sep = "\t", header=FALSE)
all_iso_prom <- as.data.frame(all_iso_prom$V4)
colnames(all_iso_prom) <- c("Transcript")

B2_location <- rep(c("Intron", "Exon", "5primeUTR", "3primeUTR", "Promoter"), 
                   times=c(length(all_iso_introns$Transcript), length(all_iso_exons$Transcript),
                           length(all_iso_5prime$Transcript), length(all_iso_3prime$Transcript),
                           length(all_iso_prom$Transcript)))
Transcript <- rbind(all_iso_introns, all_iso_exons, all_iso_5prime, all_iso_3prime, all_iso_prom)

all_iso_b2_loc <- data.frame(Transcript, B2_location)
all_iso_b2_loc <- all_iso_b2_loc %>%  group_by(Transcript) %>% mutate(B2_location = paste0(B2_location, collapse = "_"))
all_iso_b2_loc$B2_location <- ifelse(grepl("_",all_iso_b2_loc$B2_location), "Multiple B2s", all_iso_b2_loc$B2_location)
all_iso_b2_loc <- all_iso_b2_loc %>% distinct()

##Brf1 48hpi vs. Ctrl 48hpi
brf1_ctrl_inf_b2 <- left_join(brf1_ctrl_inf, all_iso_b2_loc, by=join_by(Transcript == Transcript))
brf1_ctrl_inf_b2[is.na(brf1_ctrl_inf_b2)] <- "None"
ggplot(brf1_ctrl_inf_b2, aes(x=logFC, y=-log(adj.P.Val,10), label=B2_location)) + # -log10 conversion  
  geom_point(size = 2/5, aes(color=B2_location), show.legend=TRUE) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  ggtitle("Brf1 48hpi vs. Ctrl 48hpi") +
  geom_hline(yintercept=1.3, linetype="dashed", color="red") +
  geom_vline(xintercept=1.5, linetype="dashed", color="red") +
  geom_vline(xintercept=-1.5, linetype="dashed", color="red") +
  scale_color_manual(values=c("red", "blue", "green", "orange", "purple", "black", "brown")) +
  theme_cowplot()

ggsave("../../figures/brf1_ctrl_inf_b2_mapped_to_ALL_isoforms.png")

# ##Brfl mock vs. Ctrl mock
# brf1_ctrl_mock_b2 <- left_join(brf1_ctrl_mock, all_iso_b2_loc, by=join_by(Transcript == Transcript))
# brf1_ctrl_mock_b2[is.na(brf1_ctrl_mock_b2)] <- "None"
# ggplot(brf1_ctrl_mock_b2, aes(x=logFC, y=-log(adj.P.Val,10), label=B2_location)) + # -log10 conversion  
#   geom_point(size = 2/5, aes(color=B2_location), show.legend=TRUE) +
#   xlab(expression("log"[2]*"FC")) + 
#   ylab(expression("-log"[10]*"FDR")) +
#   ggtitle("Brf1 mock vs. Ctrl mock") +
#   geom_hline(yintercept=1.3, linetype="dashed", color="red") +
#   geom_vline(xintercept=1.5, linetype="dashed", color="red") +
#   geom_vline(xintercept=-1.5, linetype="dashed", color="red") +
#   scale_color_manual(values=c("black", "blue", "green", "orange")) +
#   theme_cowplot()
# 
# ggsave("../../figures/brf1_ctrl_mock_b2_mapped_to_ALL_isoforms.png")

# ##Brf1 48hpi vs. Brf1 mock
# brf1_48_mock_b2 <- left_join(brf1_48_mock, all_iso_b2_loc, by=join_by(Transcript == Transcript))
# brf1_48_mock_b2[is.na(brf1_48_mock_b2)] <- "None"
# ggplot(brf1_48_mock_b2, aes(x=logFC, y=-log(adj.P.Val,10), label=B2_location)) + # -log10 conversion  
#   geom_point(size = 2/5, aes(color=B2_location), show.legend=TRUE) +
#   xlab(expression("log"[2]*"FC")) + 
#   ylab(expression("-log"[10]*"FDR")) +
#   ggtitle("Brf1 48hpi vs. Brf1 mock") +
#   geom_hline(yintercept=1.3, linetype="dashed", color="red") +
#   geom_vline(xintercept=1.5, linetype="dashed", color="red") +
#   geom_vline(xintercept=-1.5, linetype="dashed", color="red") +
#   scale_color_manual(values=c("black", "blue", "green", "orange")) +
#   theme_cowplot()
# 
# ggsave("../../figures/brf1_48_mock_b2_mapped_to_ALL_isoforms.png")
# 
# ##Ctrl 48hpi vs. Ctrl mock
# ctrl_48_mock_b2 <- left_join(ctrl_48_mock, all_iso_b2_loc, by=join_by(Transcript == Transcript))
# ctrl_48_mock_b2[is.na(ctrl_48_mock_b2)] <- "None"
# ggplot(ctrl_48_mock_b2, aes(x=logFC, y=-log(adj.P.Val,10), label=B2_location)) + # -log10 conversion  
#   geom_point(size = 2/5, aes(color=B2_location), show.legend=TRUE) +
#   xlab(expression("log"[2]*"FC")) + 
#   ylab(expression("-log"[10]*"FDR")) +
#   ggtitle("Ctrl 48hpi vs. Ctrl mock") +
#   geom_hline(yintercept=1.3, linetype="dashed", color="red") +
#   geom_vline(xintercept=1.5, linetype="dashed", color="red") +
#   geom_vline(xintercept=-1.5, linetype="dashed", color="red") +
#   scale_color_manual(values=c("black", "blue", "green", "orange")) +
#   theme_cowplot()
# 
# ggsave("../../figures/ctrl_48_mock_b2_mapped_to_ALL_isoforms.png")
