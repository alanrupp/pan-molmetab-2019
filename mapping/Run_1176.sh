# Analyzing data for Stat3KO mice (LeprCre/Cre;Stat3fl/fl) vs controls
# 50 bp SE reads

# Run FASTQC
mkdir FASTQC
SAMPLES=`ls | grep -E Sample_[0-9]{5}`
for f in $SAMPLES
do
cd ${f}
mkdir ../FASTQC/${f}
echo 'FASTQC on' $f
zcat *.fastq.gz | fastqc -t 8 -o ../FASTQC/${f} stdin
cd ../
done

# Remove low quality reads (Phred 20)
for f in $SAMPLES
do
date
cd ${f}
echo 'Filtering' ${f}
zcat *.fastq.gz | fastq_quality_filter -q 20 -z -o ${f}_filtered.fastq.gz
cd ../
done

# Map filtered reads to STAR
for f in $SAMPLES
do
date
cd ${f}
echo 'Aligning' ${f}
~/STAR/bin/Linux_x86_64/STAR --runThreadN 8 --genomeDir ~/STAR/genome_mm_Cre_GfpL10a --sjdbGTFfile ~/STAR/mm_Cre_GfpL10a.gtf --readFilesIn ${f}_filtered.fastq.gz --readFilesCommand zcat --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate
cd ../
done

## ----- rest of analysis in R ------ 
