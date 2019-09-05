library(tidyverse)
library(ggrepel)

# --- Read in pheno data ------------------------------------------------------
info <- read_csv("metadata/Run_824_mmyers.csv") %>%
  bind_rows(read_csv("metadata/Run_1176_mmyers.csv"))

pData <- select(info, SampleId, Description) %>%
  filter(!duplicated(SampleId)) %>%
  mutate(SampleId = as.character(SampleId)) %>%
  mutate(cells = ifelse(grepl('LepR', Description), 'bead','sup')) %>%
  mutate(cells = factor(cells, levels = c('sup','bead'))) %>%
  mutate(geno = ifelse(grepl('Ob', Description), 'ob',
                       ifelse(grepl('Stat3', Description), 'stat3ko', 'control'))) %>%
  mutate(geno = factor(geno, levels = c('control','stat3ko','ob')))

# remove ob/ob with leptin treatment
pData <- filter(pData, !(grepl('^Ob10hrLep', Description)))

# add a row of sample number to the pData for pairing
pData <- pData %>%
  mutate(pair = rep(1:(nrow(pData)/2), each = 2)) %>%
  mutate(pair = factor(pair))


# - Read in expression data (STAR) --------------------------------------------
files <- list.files(path = "counts",
                    pattern = 'ReadsPerGene.out.tab',
                    recursive = TRUE,
                    full.names = TRUE)
samples <- str_extract(files, "(?<=Sample_)([\\d]+)")

star <- map(files, ~read_tsv(.x, col_names = FALSE, skip = 4)) %>%
  bind_cols(.) %>%
  select(X1, starts_with('X2')) %>%
  set_names("gene", samples)

# - Make gene_ID to gene_name mapping
tx2gene <- read_csv('data/tx2gene.csv') %>%
  filter(!duplicated(gene_name))

star <- star %>%
  left_join(., select(tx2gene, gene_id, gene_name), by = c("gene" = "gene_id"))

star[is.na(star$gene_name), "gene_name"] <- star[is.na(star$gene_name), "gene"]

star <- star %>%
  select("gene_name", pData$SampleId)

# - Calculate pulldown efficiency for each experiment -------------------------
star %>%
  filter(gene_name %in% c("Lepr", "Cre", "Gfp_L10a")) %>%
  gather(-gene_name, key = "sample", value = "counts") %>%
  inner_join(., select(pData, SampleId, cells, geno, pair),
            by = c("sample" = "SampleId")) %>%
  ggplot(., aes(x = cells, y = counts, group = pair, color = geno)) + 
  geom_point() +
  geom_line() +
  facet_wrap(~ gene_name) +
  labs(y = 'Average counts', x = element_blank()) +
  theme_classic() +
  scale_y_continuous(trans = "log2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# - Heatmap -------------------------------------------------------------------
counts <-
  star %>%
  as.data.frame() %>%
  column_to_rownames("gene_name")

cpm <-
  star %>%
  mutate_if(is.numeric, ~.x /sum(.x) * 10^6)

expressed <-
  cpm %>%
  mutate_if(is.numeric, ~ .x > 1)
expressed <- rowSums(expressed[,-1]) > 4

annotation <- data.frame(
  cells = pData$cells,
  geno = pData$geno
)
rownames(annotation) <- pData$SampleId

annotation_colors <- list(
  cells = c("bead" = "darkgreen", "sup" = "gray"),
  geno = c("control" = "black", "stat3ko" = "blue", "ob" = "red")
)

cpm %>%
  filter(expressed) %>%
  inner_join(., select(tx2gene, gene_name, gene_biotype), by = "gene_name") %>%
  filter(gene_biotype == "protein_coding") %>%
  select(-gene_name, -gene_biotype) %>%
  mutate_all(., ~ log2(.x + 1)) %>%
  pheatmap::pheatmap(., show_rownames = FALSE,
           annotation_col = annotation,
           annotation_colors = annotation_colors)


# - PCA -----------------------------------------------------------------------
pca <- 
  cpm %>%
  filter(expressed) %>%
  inner_join(., select(tx2gene, gene_name, gene_biotype), by = "gene_name") %>%
  filter(gene_biotype == "protein_coding") %>%
  select(-gene_name, -gene_biotype) %>%
  mutate_all(., ~ log2(.x + 1)) %>%
  t() %>%
  prcomp(., center = TRUE, scale = TRUE)

pca$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample") %>%
  left_join(., select(pData, SampleId, cells, geno), 
            by = c("sample" = "SampleId")) %>%
  ggplot(., aes(x = PC1, y = PC2, color = cells, shape = geno)) +
  geom_point() +
  scale_color_manual(values = c("gray", "forestgreen")) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black"))


# - DESeq ---------------------------------------------------------------------
suppressPackageStartupMessages(library(DESeq2))
# Find enriched genes (bead vs. sup) across all datasets
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = pData, 
                              design = ~ pair + cells)
dds <- DESeq(dds)

enrich <-
  dds %>%
  results(name = c("cells_bead_vs_sup")) %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  filter(gene %in% filter(tx2gene, gene_biotype == "protein_coding")$gene_name) %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  arrange(desc(padj < 0.05), desc(log2FoldChange))

# Additional enriched genes that are not in broader model
pData <- pData %>%
  mutate(smash = paste(geno, cells, sep = '_')) %>%
  mutate(smash = factor(smash, levels = c('control_sup','control_bead',
                                          'stat3ko_sup','stat3ko_bead',
                                          'ob_sup','ob_bead')))

dds_s <- DESeqDataSetFromMatrix(countData = counts, 
                                colData = pData, 
                                design = ~ smash)
dds_s <- DESeq(dds_s)

enrich_ctrl <- 
  dds_s %>%
  results(contrast = c("smash", "control_bead", "control_sup")) %>%
  as.data.frame() %>% 
  rownames_to_column("gene") %>%
  filter(gene %in% filter(tx2gene, gene_biotype == "protein_coding")$gene_name) %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  arrange(desc(padj < 0.05), desc(log2FoldChange))

enrich_ko <- 
  dds_s %>%
  results(contrast = c('smash','stat3ko_bead','stat3ko_sup')) %>%
  as.data.frame() %>% 
  rownames_to_column("gene") %>%
  filter(gene %in% filter(tx2gene, gene_biotype == "protein_coding")$gene_name) %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  arrange(desc(padj < 0.05), desc(log2FoldChange))

enrich_ob <- 
  dds_s %>%
  results(contrast = c('smash','ob_bead','ob_sup')) %>%
  as.data.frame() %>% 
  rownames_to_column("gene") %>%
  filter(gene %in% filter(tx2gene, gene_biotype == "protein_coding")$gene_name) %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  arrange(desc(padj < 0.05), desc(log2FoldChange))

enriched <- filter(enrich, log2FoldChange >= log2(1.5) & padj < 0.05)$gene %>%
  union(., filter(enrich_ctrl, log2FoldChange >= log2(1.5) & padj < 0.05)$gene) %>%
  union(., filter(enrich_ko, log2FoldChange >= log2(1.5) & padj < 0.05)$gene) %>%
  union(., filter(enrich_ob, log2FoldChange >= log2(1.5) & padj < 0.05)$gene)



# find genes regulated in stat3 vs. control or ob vs. control
reg_ctrl_v_ko <- 
  dds_s %>%
  results(contrast = c('smash', 'stat3ko_bead', 'control_bead')) %>%
  as.data.frame() %>% 
  rownames_to_column('gene') %>%
  filter(gene %in% enriched) %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  arrange(desc(padj < 0.05), desc(log2FoldChange))

reg_ctrl_v_ob <- 
  dds_s %>%
  results(contrast = c('smash', 'ob_bead', 'control_bead')) %>%
  as.data.frame(.) %>% 
  tibble::rownames_to_column('gene') %>%
  filter(gene %in% enriched) %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  arrange(desc(padj < 0.05), desc(log2FoldChange))

reg_ob_v_ko <- 
  dds_s %>%
  results(contrast = c('smash', 'ob_bead', 'stat3ko_bead')) %>%
  as.data.frame(.) %>% 
  tibble::rownames_to_column('gene') %>%
  filter(gene %in% enriched) %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  arrange(desc(padj < 0.05), desc(log2FoldChange))

# gene lists
reg_ctrl_v_ko_gene <- filter(reg_ctrl_v_ko, padj < 0.05)$gene
reg_ctrl_v_ob_gene <- filter(reg_ctrl_v_ob, padj < 0.05)$gene
reg_ob_v_ko_gene <- filter(reg_ob_v_ko, padj < 0.05)$gene

regulated <- union(reg_ctrl_v_ko_gene, reg_ctrl_v_ob_gene) %>%
  union(., reg_ob_v_ko_gene)

fpm <-
  dds %>%
  fpm() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene")

fpm_means <-
  fpm %>%
  gather(-gene, key = "sample", value = "fpm") %>%
  left_join(., select(pData, SampleId, smash), by = c("sample"= "SampleId")) %>%
  group_by(gene, smash) %>%
  summarize(mean = mean(fpm)) %>%
  spread(key = smash, value = mean) %>%
  ungroup()


# load Stat1 target genes from Satoh and Tabunoki, Gene Regul Syst Biol 2013
# Supplmental Table 1
stat1 <- readxl::read_xls("data/stat1_genes.xls", skip = 1) %>%
  select(`Gene Symbol`) %>% 
  dplyr::rename(gene = `Gene Symbol`) %>% 
  mutate(gene = str_to_title(gene))

genes_to_label <- c("Agrp", "Npy", "Pomc", "Cartpt", "Stat1")


# - Plots ---------------------------------------------------------------------
# WT vs. Stat3KO
wt_v_stat3_plot <-
  fpm_means %>% 
  filter(gene %in% enriched) %>%
  ggplot(aes(x = control_bead, y = stat3ko_bead, 
             color = gene %in% reg_ctrl_v_ko_gene)) + 
  geom_point(show.legend = FALSE) + 
  scale_color_manual(values = c("gray", "black")) +
  geom_text_repel(data = filter(fpm_means, gene %in% genes_to_label &
                                  gene %in% reg_ctrl_v_ko_gene),
                  aes(x = control_bead, y = stat3ko_bead, label = gene), 
                  inherit.aes = FALSE, color = 'black') +
  theme_bw() + 
  scale_x_continuous(trans = 'log2', limits = c(1, 2048), expand = c(0,0)) + 
  scale_y_continuous(trans = 'log2', limits = c(1, 2048), expand = c(0,0)) +
  geom_abline(aes(slope = 1, intercept = 0), linetype = 'dashed') +
  labs(x = "Control bead (FPM)", 
       y = expression("STAT3"^"LepRb"*"KO bead (FPM)"),
       title = expression("WT vs. STAT3"^"LepRb"*"KO")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


# ob/ob vs. Stat3KO plot
ob_v_stat3_plot <-
  fpm_means %>% 
  filter(gene %in% enriched) %>%
  ggplot(aes(x = ob_bead, y = stat3ko_bead, 
             color = gene %in% reg_ob_v_ko_gene)) + 
  geom_point(show.legend = FALSE) + 
  scale_color_manual(values = c("gray", "black")) +
  geom_text_repel(data = filter(fpm_means, gene %in% genes_to_label &
                                  gene %in% reg_ob_v_ko_gene),
                  aes(x = ob_bead, y = stat3ko_bead, label = gene), 
                  inherit.aes = FALSE, color = 'black') +
  theme_bw() + 
  scale_x_continuous(trans = 'log2', limits = c(1, 2048), expand = c(0,0)) + 
  scale_y_continuous(trans = 'log2', limits = c(1, 2048), expand = c(0,0)) +
  geom_abline(aes(slope = 1, intercept = 0), linetype = 'dashed') +
  labs(x = expression(italic("ob/ob")*" bead (FPM)"),
       y = expression("STAT3"^"LepRb"*"KO bead (FPM)"), 
       title = expression(italic("ob/ob")*" vs. STAT3"^"LepRb"*"KO")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
  
# dot plot comparing values
ob_v_stat_plot <-
  reg_ctrl_v_ko %>%
  select(gene, log2FoldChange) %>%
  dplyr::rename("stat3ko" = log2FoldChange) %>%
  left_join(., select(reg_ctrl_v_ob, gene, log2FoldChange), by = "gene") %>%
  dplyr::rename("ob" = log2FoldChange) %>%
  mutate(group = case_when(
    gene %in% union(reg_ctrl_v_ko_gene, reg_ob_v_ko_gene) &
      gene %in% stat1$gene ~ "STAT1 target",
    gene %in% union(reg_ctrl_v_ko_gene, reg_ob_v_ko_gene) ~ "Regulated",
    TRUE ~ "not Regulated"
  )) %>%
  filter(gene %in% enriched) %>%
  mutate(reg = (abs(stat3ko) + abs(ob))/2) %>%
  arrange(desc(reg)) %>%
  mutate(gene = factor(gene, levels = gene)) %>%
  ggplot(aes(x = ob, y = stat3ko, color = group)) +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  geom_vline(aes(xintercept = -log2(1.5)), linetype = "dashed") +
  geom_vline(aes(xintercept = log2(1.5)), linetype = "dashed") +
  geom_hline(aes(yintercept = -log2(1.5)), linetype = "dashed") +
  geom_hline(aes(yintercept = log2(1.5)), linetype = "dashed") +
  geom_point(show.legend = FALSE) +
  scale_color_manual(values = c("gray", "black", "red")) +
  geom_text(data = data.frame(
    label = c("I", "I'", "II", "II'", "III", "III'", "IV", "IV'"),
    x = c(8.2, -2.65, 8.2, -2.65, -0.3, -0.3, -2.65, 8.2),
    y = c(7.5, -2.2, 0.3, 0.3, 7.5, -2.2, 7.5, -2.2)
    ),
    aes(x = x, y = y, label = label), inherit.aes = FALSE) +
  labs(y = expression("log"["2"]*" Fold Change in STAT3"^"LepRb"*"KO"), 
       x = expression("log"["2"]*" Fold Change in "*italic("ob/ob"))) +
  scale_x_continuous(expand = c(0, 0), limits = c(-3, 8.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-2.5, 8)) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(vjust = 0),
        legend.position = c(0.8, 0.8),
        legend.title = element_blank(),
        legend.background = element_blank())


# Figure 1 image
cowplot::plot_grid(wt_v_stat3_plot, ob_v_stat3_plot, ob_v_stat_plot,
                   labels = c("B", "C", "D"))

ggsave("Figure1.png", width = 7.5, height = 7.5, units = "in", dpi = 600)


# - Supplemental tables -------------------------------------------------------
# Table 1: Expression values
table1 <-
  fpm_means %>%
  select("gene", contains("bead")) %>%
  filter(gene %in% enriched) %>%
  left_join(., select(enrich_ctrl, gene, log2FoldChange, padj), by = "gene") %>%
  dplyr::rename("Enrichment (WT)" = log2FoldChange,
                "P (WT)" = padj) %>%
  left_join(., select(enrich_ko, gene, log2FoldChange, padj), by = "gene") %>%
  dplyr::rename("Enrichment (STAT3KO)" = log2FoldChange,
                "P (STAT3KO)" = padj) %>%
  left_join(., select(enrich_ob, gene, log2FoldChange, padj), by = "gene") %>%
  dplyr::rename("Gene" = gene,
                "FPM (WT)" = control_bead,
                "FPM (STAT3KO)" = stat3ko_bead,
                "FPM (ob/ob)" = ob_bead,
                "Enrichment (ob/ob)" = log2FoldChange,
                "P (ob/ob)" = padj) %>%
  mutate(enrich = (`Enrichment (WT)` + `Enrichment (STAT3KO)` + 
                     `Enrichment (ob/ob)`)/3) %>%
  arrange(desc(enrich)) %>%
  select(-enrich,
         Gene, starts_with("FPM"), starts_with("Enrichment"), starts_with("P"))
  
write.csv(table1, "Table1.csv", row.names = FALSE, na = "")


# Table 2: Regulation values
table2 <-
  fpm_means %>%
  select("gene", "control_bead") %>%
  left_join(., select(enrich_ctrl, gene, log2FoldChange), by = "gene") %>%
  dplyr::rename("Enrichment (WT)" = log2FoldChange,
                "FPM (WT)" = control_bead) %>%
  select(gene, `Enrichment (WT)`, `FPM (WT)`) %>%
  left_join(., select(reg_ctrl_v_ob, gene, log2FoldChange, padj), by = "gene") %>%
  dplyr::rename("ob/ob vs. Control" = log2FoldChange,
                "P (ob-ctrl)" = padj) %>%
  left_join(., select(reg_ctrl_v_ko, gene, log2FoldChange, padj), by = "gene") %>%
  dplyr::rename("STAT3KO vs. Control" = log2FoldChange,
                "P (ko-ctrl)" = padj) %>%
  left_join(., select(reg_ob_v_ko, gene, log2FoldChange, padj), by = "gene") %>%
  dplyr::rename("STAT3KO vs. ob/ob" = log2FoldChange,
                "P (ko-ob)" = padj) %>%
  filter(gene %in% regulated) %>%
  mutate(group = case_when(
    `ob/ob vs. Control` >= log2(1.5) & 
      `STAT3KO vs. Control` >= log2(1.5) ~ "I",
    `ob/ob vs. Control` <= -log2(1.5) & 
      `STAT3KO vs. Control` <= -log2(1.5) ~ "I'",
    `ob/ob vs. Control` >= log2(1.5) & 
      `STAT3KO vs. Control` < log2(1.5) & 
      `STAT3KO vs. Control` > -log2(1.5) ~ "II",
    `ob/ob vs. Control` <= -log2(1.5) & 
      `STAT3KO vs. Control` < log2(1.5) & 
      `STAT3KO vs. Control` > -log2(1.5) ~ "II'",
    `ob/ob vs. Control` < log2(1.5) & 
      `ob/ob vs. Control` > -log2(1.5) &
      `STAT3KO vs. Control` >= log2(1.5) ~ "III",
    `ob/ob vs. Control` < log2(1.5) & 
      `ob/ob vs. Control` > -log2(1.5) &
      `STAT3KO vs. Control` <= -log2(1.5) ~ "III'",
    `ob/ob vs. Control` <= log2(1.5) & 
      `STAT3KO vs. Control` >= log2(1.5) ~ "IV",
    `ob/ob vs. Control` >= log2(1.5) & 
      `STAT3KO vs. Control` <= log2(1.5) ~ "IV'"
  )) %>%
  mutate(group = factor(group, levels = c("I", "I'", "II", "II'", "III",
                                          "III'", "IV", "IV'"))) %>%
  dplyr::rename("Gene" = gene,
                "Group" = group) %>%
  arrange(Group, Gene) %>%
  filter(!is.na(Group))

write.csv(table2, "Table2.csv", row.names = FALSE, na = "")

# - EnrichR analysis (Table 3) ------------------------------------------------
ob_not_stat3 <-
  reg_ctrl_v_ko %>%
  filter(gene %in% setdiff(reg_ctrl_v_ob_gene, reg_ctrl_v_ko_gene)) %>%
  filter(abs(log2FoldChange) <= log2(1.5)) %>%
  .$gene

stat3_not_ob <-
  reg_ctrl_v_ob %>%
  filter(gene %in% setdiff(reg_ctrl_v_ko_gene, reg_ctrl_v_ob_gene)) %>%
  filter(abs(log2FoldChange) <= log2(1.5)) %>%
  .$gene

reg_ctrl_v_ob %>%
  filter(padj < 0.05) %>% # changed in ob vs. control
  filter(gene %in% reg_ob_v_ko_gene) %>% # different in ob vs. stat3
  filter(!(gene %in% setdiff(reg_ctrl_v_ko_gene, reg_ctrl_v_ob_gene))) %>%
  # not changed in stat3 but not ob
  nrow()
