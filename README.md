# *Stat3*-regulated transcriptional changes in LepRb hypothalamic neurons

Looking for genes specifically regulated by the Leptin-STAT3 signaling pathway in hypothalamic LepRb neurons. Comparing hypothalamic *LepRb<sup>Cre</sup>* TRAP-seq data from *Stat3<sup>LepRb</sup>KO* and *ob/ob* mice compared to controls.

## Files
* `counts`: count data from `STAR` mapping
* `data`
  * `stat1_genes.xls`: *Stat1*-regulated genes
  * `tx2gene.csv`: transcript_id, transcript_name, gene_id, gene_name, and gene_biotype information
* `mapping`: folder containing shell scripts to filter and map original FASTQ data
* `metadata`: folder containing CSV files with sample and run information
* `Stat3KO_v_ob.R`: analysis of count files

## Requirements
* `R`
  * `tidyverse`
  * `DESeq2`
  * `ggrepel`

## Outputs
Data and figures from Pan et al. (see Citation)
* `Figure1.png`: Figure 1 B-D
* `Table1.csv`: Supplemental Table 1
* `Table2.csv`: Supplemental Table 2

## Citation
Pan W et al. Transcriptional and physiological roles for STAT proteins in leptin action. *Mol Metab*. 2019 Apr;**22**:121-131. doi: 10.1016/j.molmet.2019.01.007.
