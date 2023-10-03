install.packages("manhattanly")

library(manhattanly)

setwd("/Users/sahil/Documents/b2_sines/archive_2022_07_18__12_12_14/Genes_DE")

#Gene-level
brf1_48_mock <- read.delim("glaunsinger_Brf1_48hpi_v_Brf1_Mock_genes_DE_analysis_be.txt", sep = "\t")
ctrl_48_mock <- read.delim("glaunsinger_Ctrl_48_v_Ctrl_Mock_genes_DE_analysis_be.txt", sep = "\t")
ctrl_brf1_mock <- read.delim("glaunsinger_Brf1_Mock_v_Ctrl_Mock_genes_DE_analysis_be.txt", sep = "\t")
ctrl_brf1_48 <- read.delim("glaunsinger_Brf1_48hpi_v_Ctrl_48_genes_DE_analysis_be.txt", sep = "\t")
#full_df <- rbind(brf1_48_mock, ctrl_48_mock)

#Brf1 conditions
volcanorObj <- volcanor(brf1_48_mock,
                        p = "adj.P.Val",
                        effect_size = "logFC",
                        snp = "Transcript"
)

volcanoly(volcanorObj, gene = "Transcript", genomewideline = -log10(5e-2), effect_size_line = FALSE)

#Ctrl conditions
volcanorObj <- volcanor(ctrl_48_mock,
                        p = "adj.P.Val",
                        effect_size = "logFC",
                        snp = "Transcript"
)

volcanoly(volcanorObj, gene = "Transcript", genomewideline = -log10(5e-2), effect_size_line = FALSE)

#Brf1 Mock vs Ctrl Mock
volcanorObj <- volcanor(ctrl_brf1_mock,
                        p = "adj.P.Val",
                        effect_size = "logFC",
                        snp = "Transcript"
)

volcanoly(volcanorObj, gene = "Transcript", genomewideline = -log10(1e-2), effect_size_line = FALSE)

#Brf1 vs Ctrl 48hpi conditions
volcanorObj <- volcanor(ctrl_brf1_48,
                        p = "adj.P.Val",
                        effect_size = "logFC",
                        snp = "Transcript"
)

volcanoly(volcanorObj, gene = "Transcript", genomewideline = -log10(5e-2), effect_size_line = FALSE)


###############################################MHV68
setwd("/Users/sahil/Documents/b2_sines/archive_2022_07_18__12_12_14/Genes_DE/mhv68")
mhv68_genes <- read.delim("glaunsinger_Brf1_48hpi_v_Ctrl_48_genes_DE_mhv68_analysis_be.txt", sep = "\t")
mhv68_isoforms <- read.delim("glaunsinger_Brf1_48hpi_v_Ctrl_48_DE_mhv68_analysis_be.txt", sep = "\t")

#Brf1 conditions
volcanorObj <- volcanor(mhv68_genes,
                        p = "adj.P.Val",
                        effect_size = "logFC",
                        snp = "Transcript"
)

volcanoly(volcanorObj, gene = "Transcript", genomewideline = -log10(5e-2), effect_size_line = FALSE)

#Ctrl conditions
volcanorObj <- volcanor(mhv68_isoforms,
                        p = "adj.P.Val",
                        effect_size = "logFC",
                        snp = "Transcript"
)

volcanoly(volcanorObj, gene = "Transcript", genomewideline = -log10(5e-2), effect_size_line = FALSE)



####################################################################
#Filtered genes
setwd("/Users/sahil/Documents/b2_sines/archive_2022_07_18__12_12_14/Genes_DE/filteredGenes")
brf1_48_mock_filtered <- read.delim("glaunsinger_Brf1_48hpi_v_Brf1_Mock_DE_filteredGenes_analysis.txt", sep = "\t")
ctrl_48_mock_filtered <- read.delim("glaunsinger_Ctrl_48_v_Ctrl_Mock_DE_filteredGenes_analysis.txt", sep = "\t")

#Brf1 conditions
volcanorObj <- volcanor(brf1_48_mock_filtered,
                        p = "adj.P.Val",
                        effect_size = "logFC",
                        snp = "Gene"
)

volcanoly(volcanorObj, gene = "Gene", genomewideline = -log10(5e-2), effect_size_line = FALSE)

#Ctrl conditions
volcanorObj <- volcanor(ctrl_48_mock_filtered,
                        p = "adj.P.Val",
                        effect_size = "logFC",
                        snp = "Gene"
)

volcanoly(volcanorObj, gene = "Gene", genomewideline = -log10(5e-2), effect_size_line = FALSE)

#######################################################################
#Transcript-level
brf1_48_mock_iso <- read.delim("glaunsinger_Brf1_48hpi_v_Brf1_Mock_DE_analysis_mm10_be.txt", sep = "\t")
ctrl_48_mock_iso <- read.delim("glaunsinger_Ctrl_48_v_Ctrl_Mock_DE_analysis_mm10_be.txt", sep = "\t")
brf1_ctrl_48_iso <- read.delim("glaunsinger_Brf1_48hpi_v_Ctrl_48_DE_analysis_mm10_be.txt", sep = "\t")
brf1_ctrl_mock_iso <- read.delim("glaunsinger_Brf1_Mock_v_Ctrl_Mock_DE_analysis_mm10_be.txt", sep = "\t")
#full_df <- rbind(brf1_48_mock, ctrl_48_mock)

#Brf1 vs Ctrl 48hpi condition
volcanorObj <- volcanor(brf1_ctrl_48_iso,
                        p = "adj.P.Val",
                        effect_size = "logFC",
                        snp = "Transcript"
)

volcanoly(volcanorObj, gene = "Transcript", genomewideline = -log10(5e-2), effect_size_line = FALSE)

#Brf1 Mock vs Ctrl Mock
volcanorObj <- volcanor(brf1_ctrl_mock_iso,
                        p = "adj.P.Val",
                        effect_size = "logFC",
                        snp = "Transcript"
)

volcanoly(volcanorObj, gene = "Transcript", genomewideline = -log10(1e-2), effect_size_line = FALSE)

#Brf1 conditions
volcanorObj <- volcanor(brf1_48_mock,
                        p = "adj.P.Val",
                        effect_size = "logFC",
                        snp = "Transcript"
)

volcanoly(volcanorObj, gene = "Transcript", genomewideline = -log10(5e-2), effect_size_line = FALSE)

#Ctrl conditions
volcanorObj <- volcanor(ctrl_48_mock_iso,
                        p = "adj.P.Val",
                        effect_size = "logFC",
                        snp = "Transcript"
)

volcanoly(volcanorObj, gene = "Transcript", genomewideline = -log10(5e-2), effect_size_line = FALSE)

