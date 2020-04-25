#!/bin/sh

cd /home/shriyas/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/rep2_oct2019/
#mkdir reads
#cd reads

#gunzip -c ~/cifs-lab/sequencing_dataset_archive/Chloroplast_Hi-C_IPOD_MNase-H3_DNAmethyl/ruddle.brcf.med.umich.edu/NovaA-147/wierzbicki/Sample_140767/140767_TGATCCTT-ATTGCATC_S113_R1_001.fastq.gz > 140767_TGATCCTT-ATTGCATC_S113_R1_001.fastq
#gunzip -c ~/cifs-lab/sequencing_dataset_archive/Chloroplast_Hi-C_IPOD_MNase-H3_DNAmethyl/ruddle.brcf.med.umich.edu/NovaA-147/wierzbicki/Sample_140767/140767_TGATCCTT-ATTGCATC_S113_R2_001.fastq.gz > 140767_TGATCCTT-ATTGCATC_S113_R2_001.fastq

#gunzip -c ~/cifs-lab/sequencing_dataset_archive/Chloroplast_Hi-C_IPOD_MNase-H3_DNAmethyl/ruddle.brcf.med.umich.edu/NovaA-147/wierzbicki/Sample_140768/140768_ACTCAGCA-CGCCTAAC_S114_R1_001.fastq.gz > 140768_ACTCAGCA-CGCCTAAC_S114_R1_001.fastq
#gunzip -c ~/cifs-lab/sequencing_dataset_archive/Chloroplast_Hi-C_IPOD_MNase-H3_DNAmethyl/ruddle.brcf.med.umich.edu/NovaA-147/wierzbicki/Sample_140768/140768_ACTCAGCA-CGCCTAAC_S114_R2_001.fastq.gz > 140768_ACTCAGCA-CGCCTAAC_S114_R2_001.fastq

#gunzip -c ~/cifs-lab/sequencing_dataset_archive/Chloroplast_Hi-C_IPOD_MNase-H3_DNAmethyl/ruddle.brcf.med.umich.edu/NovaA-147/wierzbicki/Sample_140769/140769_GTAGTACA-CTCGCTTA_S115_R1_001.fastq.gz > 140769_GTAGTACA-CTCGCTTA_S115_R1_001.fastq
#gunzip -c ~/cifs-lab/sequencing_dataset_archive/Chloroplast_Hi-C_IPOD_MNase-H3_DNAmethyl/ruddle.brcf.med.umich.edu/NovaA-147/wierzbicki/Sample_140769/140769_GTAGTACA-CTCGCTTA_S115_R2_001.fastq.gz > 140769_GTAGTACA-CTCGCTTA_S115_R2_001.fastq

#gunzip -c ~/cifs-lab/sequencing_dataset_archive/Chloroplast_Hi-C_IPOD_MNase-H3_DNAmethyl/ruddle.brcf.med.umich.edu/NovaA-147/wierzbicki/Sample_140770/140770_CGAAGGCG-TTCCGCGA_S116_R1_001.fastq.gz > 140770_CGAAGGCG-TTCCGCGA_S116_R1_001.fastq
#gunzip -c ~/cifs-lab/sequencing_dataset_archive/Chloroplast_Hi-C_IPOD_MNase-H3_DNAmethyl/ruddle.brcf.med.umich.edu/NovaA-147/wierzbicki/Sample_140770/140770_CGAAGGCG-TTCCGCGA_S116_R2_001.fastq.gz > 140770_CGAAGGCG-TTCCGCGA_S116_R2_001.fastq

#gunzip -c ~/cifs-lab/sequencing_dataset_archive/Chloroplast_Hi-C_IPOD_MNase-H3_DNAmethyl/ruddle.brcf.med.umich.edu/NovaA-147/wierzbicki/Sample_140771/140771_AGTACTTA-ACAGGGAA_S117_R1_001.fastq.gz > 140771_AGTACTTA-ACAGGGAA_S117_R1_001.fastq
#gunzip -c ~/cifs-lab/sequencing_dataset_archive/Chloroplast_Hi-C_IPOD_MNase-H3_DNAmethyl/ruddle.brcf.med.umich.edu/NovaA-147/wierzbicki/Sample_140771/140771_AGTACTTA-ACAGGGAA_S117_R2_001.fastq.gz > 140771_AGTACTTA-ACAGGGAA_S117_R2_001.fastq

#gunzip -c ~/cifs-lab/sequencing_dataset_archive/Chloroplast_Hi-C_IPOD_MNase-H3_DNAmethyl/ruddle.brcf.med.umich.edu/NovaA-147/wierzbicki/Sample_140769/140769_GTAGTACA-CTCGCTTA_S115_R2_001.fastq.gz > 140769_GTAGTACA-CTCGCTTA_S115_R2_001.fastq

#fastqc --casava ~/cifs-lab/sequencing_dataset_archive/Chloroplast_Hi-C_IPOD_MNase-H3_DNAmethyl/ruddle.brcf.med.umich.edu/NovaA-147/wierzbicki/Sample_140767/*fastq
#fastqc --casava ~/cifs-lab/sequencing_dataset_archive/Chloroplast_Hi-C_IPOD_MNase-H3_DNAmethyl/ruddle.brcf.med.umich.edu/NovaA-147/wierzbicki/Sample_140768/*fastq
#fastqc --casava ~/cifs-lab/sequencing_dataset_archive/Chloroplast_Hi-C_IPOD_MNase-H3_DNAmethyl/ruddle.brcf.med.umich.edu/NovaA-147/wierzbicki/Sample_140769/*fastq
#fastqc --casava ~/cifs-lab/sequencing_dataset_archive/Chloroplast_Hi-C_IPOD_MNase-H3_DNAmethyl/ruddle.brcf.med.umich.edu/NovaA-147/wierzbicki/Sample_140770/*fastq
#fastqc --casava ~/cifs-lab/sequencing_dataset_archive/Chloroplast_Hi-C_IPOD_MNase-H3_DNAmethyl/ruddle.brcf.med.umich.edu/NovaA-147/wierzbicki/Sample_140771/*fastq

#cd ../
#mkdir trimmed_reads
#cd trimmed_reads

#trim_galore --length 35 --fastqc --paired ~/cifs-lab/sequencing_dataset_archive/Chloroplast_Hi-C_IPOD_MNase-H3_DNAmethyl/ruddle.brcf.med.umich.edu/NovaA-147/wierzbicki/Sample_140767/140767_TGATCCTT-ATTGCATC_S113_R1_001.fastq ~/cifs-lab/sequencing_dataset_archive/Chloroplast_Hi-C_IPOD_MNase-H3_DNAmethyl/ruddle.brcf.med.umich.edu/NovaA-147/wierzbicki/Sample_140767/140767_TGATCCTT-ATTGCATC_S113_R2_001.fastq
#trim_galore --length 35 --fastqc --paired ~/cifs-lab/sequencing_dataset_archive/Chloroplast_Hi-C_IPOD_MNase-H3_DNAmethyl/ruddle.brcf.med.umich.edu/NovaA-147/wierzbicki/Sample_140768/140768_ACTCAGCA-CGCCTAAC_S114_R1_001.fastq ~/cifs-lab/sequencing_dataset_archive/Chloroplast_Hi-C_IPOD_MNase-H3_DNAmethyl/ruddle.brcf.med.umich.edu/NovaA-147/wierzbicki/Sample_140768/140768_ACTCAGCA-CGCCTAAC_S114_R2_001.fastq
#trim_galore --length 35 --fastqc --paired ~/cifs-lab/sequencing_dataset_archive/Chloroplast_Hi-C_IPOD_MNase-H3_DNAmethyl/ruddle.brcf.med.umich.edu/NovaA-147/wierzbicki/Sample_140769/140769_GTAGTACA-CTCGCTTA_S115_R1_001.fastq ~/cifs-lab/sequencing_dataset_archive/Chloroplast_Hi-C_IPOD_MNase-H3_DNAmethyl/ruddle.brcf.med.umich.edu/NovaA-147/wierzbicki/Sample_140769/140769_GTAGTACA-CTCGCTTA_S115_R2_001.fastq

#trim_galore --length 35 --fastqc --paired ~/cifs-lab/sequencing_dataset_archive/Chloroplast_Hi-C_IPOD_MNase-H3_DNAmethyl/ruddle.brcf.med.umich.edu/NovaA-147/wierzbicki/Sample_140770/140770_CGAAGGCG-TTCCGCGA_S116_R1_001.fastq ~/cifs-lab/sequencing_dataset_archive/Chloroplast_Hi-C_IPOD_MNase-H3_DNAmethyl/ruddle.brcf.med.umich.edu/NovaA-147/wierzbicki/Sample_140770/140770_CGAAGGCG-TTCCGCGA_S116_R2_001.fastq
#trim_galore --length 35 --fastqc --paired ~/cifs-lab/sequencing_dataset_archive/Chloroplast_Hi-C_IPOD_MNase-H3_DNAmethyl/ruddle.brcf.med.umich.edu/NovaA-147/wierzbicki/Sample_140771/140771_AGTACTTA-ACAGGGAA_S117_R1_001.fastq ~/cifs-lab/sequencing_dataset_archive/Chloroplast_Hi-C_IPOD_MNase-H3_DNAmethyl/ruddle.brcf.med.umich.edu/NovaA-147/wierzbicki/Sample_140771/140771_AGTACTTA-ACAGGGAA_S117_R2_001.fastq

#cd ../
#mkdir mapped_reads
cd mapped_reads/

#bowtie2 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -p 6 -q --no-mixed --no-discordant --no-overlap --no-contain --fr -1 ../trimmed_reads/140767_TGATCCTT-ATTGCATC_S113_R1_001_val_1.fq -2 ../trimmed_reads/140767_TGATCCTT-ATTGCATC_S113_R2_001_val_2.fq -S col0_rep2_methyl_MNase.sam
#26452847 reads; of these:
#  26452847 (100.00%) were paired; of these:
#    18959470 (71.67%) aligned concordantly 0 times
#    2035210 (7.69%) aligned concordantly exactly 1 time
#    5458167 (20.63%) aligned concordantly >1 times
#28.33% overall alignment rate

#bowtie2 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -p 6 -q --no-mixed --no-discordant --no-overlap --no-contain --fr -1 ../trimmed_reads/140768_ACTCAGCA-CGCCTAAC_S114_R1_001_val_1.fq -2 ../trimmed_reads/140768_ACTCAGCA-CGCCTAAC_S114_R2_001_val_2.fq -S met13_rep2_methyl_MNase.sam
#21530384 reads; of these:
#  21530384 (100.00%) were paired; of these:
#    12501096 (58.06%) aligned concordantly 0 times
#    4244023 (19.71%) aligned concordantly exactly 1 time
#    4785265 (22.23%) aligned concordantly >1 times
#41.94% overall alignment rate

#bowtie2 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -p 6 -q --no-mixed --no-discordant --no-overlap --no-contain --fr -1 ../trimmed_reads/140769_GTAGTACA-CTCGCTTA_S115_R1_001_val_1.fq -2 ../trimmed_reads/140769_GTAGTACA-CTCGCTTA_S115_R2_001_val_2.fq -S cmt2_rep2_methyl_MNase.sam
#22509497 reads; of these:
#  22509497 (100.00%) were paired; of these:
#    14249919 (63.31%) aligned concordantly 0 times
#    2950779 (13.11%) aligned concordantly exactly 1 time
#    5308799 (23.58%) aligned concordantly >1 times
#36.69% overall alignment rate

#bowtie2 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -p 6 -q --no-mixed --no-discordant --no-overlap --no-contain --fr -1 ../trimmed_reads/140770_CGAAGGCG-TTCCGCGA_S116_R1_001_val_1.fq -2 ../trimmed_reads/140770_CGAAGGCG-TTCCGCGA_S116_R2_001_val_2.fq -S cmt3_rep2_methyl_MNase.sam
#24749847 reads; of these:
#  24749847 (100.00%) were paired; of these:
#    15968968 (64.52%) aligned concordantly 0 times
#    3491724 (14.11%) aligned concordantly exactly 1 time
#    5289155 (21.37%) aligned concordantly >1 times
#35.48% overall alignment rate

#bowtie2 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -p 6 -q --no-mixed --no-discordant --no-overlap --no-contain --fr -1 ../trimmed_reads/140771_AGTACTTA-ACAGGGAA_S117_R1_001_val_1.fq -2 ../trimmed_reads/140771_AGTACTTA-ACAGGGAA_S117_R2_001_val_2.fq -S drm2_rep2_methyl_MNase.sam
#29255006 reads; of these:
#  29255006 (100.00%) were paired; of these:
#    22405042 (76.59%) aligned concordantly 0 times
#    2676810 (9.15%) aligned concordantly exactly 1 time
#    4173154 (14.26%) aligned concordantly >1 times
#23.41% overall alignment rate

bowtie2 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -p 6 -q --no-mixed --no-discordant --fr -1 ../trimmed_reads/140767_TGATCCTT-ATTGCATC_S113_R1_001_val_1.fq -2 ../trimmed_reads/140767_TGATCCTT-ATTGCATC_S113_R2_001_val_2.fq -S col0_rep2_methyl_MNase.sam
26452847 reads; of these:
  26452847 (100.00%) were paired; of these:
    1488649 (5.63%) aligned concordantly 0 times
    6338884 (23.96%) aligned concordantly exactly 1 time
    18625314 (70.41%) aligned concordantly >1 times
94.37% overall alignment rate

bowtie2 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -p 6 -q --no-mixed --no-discordant --fr -1 ../trimmed_reads/140768_ACTCAGCA-CGCCTAAC_S114_R1_001_val_1.fq -2 ../trimmed_reads/140768_ACTCAGCA-CGCCTAAC_S114_R2_001_val_2.fq -S met13_rep2_methyl_MNase.sam
21530384 reads; of these:
  21530384 (100.00%) were paired; of these:
    4565559 (21.21%) aligned concordantly 0 times
    7279900 (33.81%) aligned concordantly exactly 1 time
    9684925 (44.98%) aligned concordantly >1 times
78.79% overall alignment rate

bowtie2 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -p 6 -q --no-mixed --no-discordant --fr -1 ../trimmed_reads/140769_GTAGTACA-CTCGCTTA_S115_R1_001_val_1.fq -2 ../trimmed_reads/140769_GTAGTACA-CTCGCTTA_S115_R2_001_val_2.fq -S cmt2_rep2_methyl_MNase.sam
22509497 reads; of these:
  22509497 (100.00%) were paired; of these:
    2307039 (10.25%) aligned concordantly 0 times
    6533649 (29.03%) aligned concordantly exactly 1 time
    13668809 (60.72%) aligned concordantly >1 times
89.75% overall alignment rate

bowtie2 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -p 6 -q --no-mixed --no-discordant --fr -1 ../trimmed_reads/140770_CGAAGGCG-TTCCGCGA_S116_R1_001_val_1.fq -2 ../trimmed_reads/140770_CGAAGGCG-TTCCGCGA_S116_R2_001_val_2.fq -S cmt3_rep2_methyl_MNase.sam
24749847 reads; of these:
  24749847 (100.00%) were paired; of these:
    2051762 (8.29%) aligned concordantly 0 times
    8257996 (33.37%) aligned concordantly exactly 1 time
    14440089 (58.34%) aligned concordantly >1 times
91.71% overall alignment rate

bowtie2 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -p 6 -q --no-mixed --no-discordant --fr -1 ../trimmed_reads/140771_AGTACTTA-ACAGGGAA_S117_R1_001_val_1.fq -2 ../trimmed_reads/140771_AGTACTTA-ACAGGGAA_S117_R2_001_val_2.fq -S drm2_rep2_methyl_MNase.sam
29255006 reads; of these:
  29255006 (100.00%) were paired; of these:
    1885338 (6.44%) aligned concordantly 0 times
    10581311 (36.17%) aligned concordantly exactly 1 time
    16788357 (57.39%) aligned concordantly >1 times
93.56% overall alignment rate

mkdir process_map_reads
cd process_map_reads

java -jar $PICARD SortSam INPUT=../col0_rep2_methyl_MNase.sam OUTPUT=col0_rep2_methyl_MNase_sorted.bam SORT_ORDER=coordinate
java -jar $PICARD MarkDuplicates REMOVE_DUPLICATES=true I=col0_rep2_methyl_MNase_sorted.bam O=col0_rep2_methyl_MNase_sorted_dedup.bam M=col0_h3_lo_dup_metrics.txt
java -jar $PICARD BuildBamIndex INPUT=col0_rep2_methyl_MNase_sorted_dedup.bam
samtools sort -n col0_rep2_methyl_MNase_sorted_dedup.bam col0_rep2_methyl_MNase_sorted_dedup_sn
bamToBed -bedpe -i col0_rep2_methyl_MNase_sorted_dedup_sn.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > col0_rep2_methyl_MNase_sorted_dedup.bed

java -jar $PICARD SortSam INPUT=../met13_rep2_methyl_MNase.sam OUTPUT=met13_rep2_methyl_MNase_sorted.bam SORT_ORDER=coordinate
java -jar $PICARD MarkDuplicates REMOVE_DUPLICATES=true I=met13_rep2_methyl_MNase_sorted.bam O=met13_rep2_methyl_MNase_sorted_dedup.bam M=met13_h3_lo_dup_metrics.txt
java -jar $PICARD BuildBamIndex INPUT=met13_rep2_methyl_MNase_sorted_dedup.bam
samtools sort -n met13_rep2_methyl_MNase_sorted_dedup.bam met13_rep2_methyl_MNase_sorted_dedup_sn
bamToBed -bedpe -i met13_rep2_methyl_MNase_sorted_dedup_sn.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > met13_rep2_methyl_MNase_sorted_dedup.bed

java -jar $PICARD SortSam INPUT=../cmt2_rep2_methyl_MNase.sam OUTPUT=cmt2_rep2_methyl_MNase_sorted.bam SORT_ORDER=coordinate
java -jar $PICARD MarkDuplicates REMOVE_DUPLICATES=true I=cmt2_rep2_methyl_MNase_sorted.bam O=cmt2_rep2_methyl_MNase_sorted_dedup.bam M=cmt2_h3_lo_dup_metrics.txt
java -jar $PICARD BuildBamIndex INPUT=cmt2_rep2_methyl_MNase_sorted_dedup.bam
samtools sort -n cmt2_rep2_methyl_MNase_sorted_dedup.bam cmt2_rep2_methyl_MNase_sorted_dedup_sn
bamToBed -bedpe -i cmt2_rep2_methyl_MNase_sorted_dedup_sn.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > cmt2_rep2_methyl_MNase_sorted_dedup.bed

java -jar $PICARD SortSam INPUT=../cmt3_rep2_methyl_MNase.sam OUTPUT=cmt3_rep2_methyl_MNase_sorted.bam SORT_ORDER=coordinate
java -jar $PICARD MarkDuplicates REMOVE_DUPLICATES=true I=cmt3_rep2_methyl_MNase_sorted.bam O=cmt3_rep2_methyl_MNase_sorted_dedup.bam M=cmt3_h3_lo_dup_metrics.txt
java -jar $PICARD BuildBamIndex INPUT=cmt3_rep2_methyl_MNase_sorted_dedup.bam
samtools sort -n cmt3_rep2_methyl_MNase_sorted_dedup.bam cmt3_rep2_methyl_MNase_sorted_dedup_sn
bamToBed -bedpe -i cmt3_rep2_methyl_MNase_sorted_dedup_sn.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > cmt3_rep2_methyl_MNase_sorted_dedup.bed

java -jar $PICARD SortSam INPUT=../drm2_rep2_methyl_MNase.sam OUTPUT=drm2_rep2_methyl_MNase_sorted.bam SORT_ORDER=coordinate
java -jar $PICARD MarkDuplicates REMOVE_DUPLICATES=true I=drm2_rep2_methyl_MNase_sorted.bam O=drm2_rep2_methyl_MNase_sorted_dedup.bam M=drm2_h3_lo_dup_metrics.txt
java -jar $PICARD BuildBamIndex INPUT=drm2_rep2_methyl_MNase_sorted_dedup.bam
samtools sort -n drm2_rep2_methyl_MNase_sorted_dedup.bam drm2_rep2_methyl_MNase_sorted_dedup_sn
bamToBed -bedpe -i drm2_rep2_methyl_MNase_sorted_dedup_sn.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > drm2_rep2_methyl_MNase_sorted_dedup.bed


samtools view -bq 2 col0_rep2_methyl_MNase_sorted_dedup_sn.bam > col0_rep2_methyl_MNase_sorted_dedup_sn_mapq2.bam
samtools view -bq 2 drm2_rep2_methyl_MNase_sorted_dedup_sn.bam > drm2_rep2_methyl_MNase_sorted_dedup_sn_mapq2.bam


awk '{if($5>=2) print $0}' col0_rep2_methyl_MNase_sorted_dedup.bed > col0_rep2_methyl_MNase_sorted_dedup_mapq.bed ### filter for mapq>2
awk '{if($5>=2) print $0}' met13_rep2_methyl_MNase_sorted_dedup.bed > met13_rep2_methyl_MNase_sorted_dedup_mapq.bed ### filter for mapq>2
awk '{if($5>=2) print $0}' cmt2_rep2_methyl_MNase_sorted_dedup.bed > cmt2_rep2_methyl_MNase_sorted_dedup_mapq.bed ### filter for mapq>2
awk '{if($5>=2) print $0}' cmt3_rep2_methyl_MNase_sorted_dedup.bed > cmt3_rep2_methyl_MNase_sorted_dedup_mapq.bed ### filter for mapq>2
awk '{if($5>=2) print $0}' drm2_rep2_methyl_MNase_sorted_dedup.bed > drm2_rep2_methyl_MNase_sorted_dedup_mapq.bed ### filter for mapq>2





samtools view -bS -o col0_rep2_methyl_MNase.bam col0_rep2_methyl_MNase.sam
samtools view -bS -o met13_rep2_methyl_MNase.bam met13_rep2_methyl_MNase.sam
samtools view -bS -o cmt2_rep2_methyl_MNase.bam cmt2_rep2_methyl_MNase.sam
samtools view -bS -o cmt3_rep2_methyl_MNase.bam cmt3_rep2_methyl_MNase.sam
samtools view -bS -o drm2_rep2_methyl_MNase.bam drm2_rep2_methyl_MNase.sam

bamToBed -bedpe -i col0_rep2_methyl_MNase.bam | awk '{if($2>0) print $0}' | sort -k1 -k2 -n | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > col0_rep2_methyl_MNase.bed
bamToBed -bedpe -i met13_rep2_methyl_MNase.bam | awk '{if($2>0) print $0}' | sort -k1 -k2 -n | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > met13_rep2_methyl_MNase.bed
bamToBed -bedpe -i cmt2_rep2_methyl_MNase.bam | awk '{if($2>0) print $0}' | sort -k1 -k2 -n | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > cmt2_rep2_methyl_MNase.bed
bamToBed -bedpe -i cmt3_rep2_methyl_MNase.bam | awk '{if($2>0) print $0}' | sort -k1 -k2 -n | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > cmt3_rep2_methyl_MNase.bed
bamToBed -bedpe -i drm2_rep2_methyl_MNase.bam | awk '{if($2>0) print $0}' | sort -k1 -k2 -n | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > drm2_rep2_methyl_MNase.bed



samtools view -h col0_rep2_methyl_MNase_sorted_dedup_sn.bam | awk -f ~/cifs-lab/Hafiz/scripts/hafiz_scripts/bam_tlen_filter2.awk | samtools view -Sb - > col0_rep2_methyl_MNase_sorted_dedup_sn_f.bam
bamToBed -bedpe -i col0_rep2_methyl_MNase_sorted_dedup_sn_f.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > col0_rep2_methyl_MNase_sorted_dedup_sn_f.bed

samtools view -h cmt2_rep2_methyl_MNase_sorted_dedup_sn.bam | awk -f ~/cifs-lab/Hafiz/scripts/hafiz_scripts/bam_tlen_filter2.awk | samtools view -Sb - > cmt2_rep2_methyl_MNase_sorted_dedup_sn_f.bam
bamToBed -bedpe -i cmt2_rep2_methyl_MNase_sorted_dedup_sn_f.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > cmt2_rep2_methyl_MNase_sorted_dedup_sn_f.bed

samtools view -h cmt3_rep2_methyl_MNase_sorted_dedup_sn.bam | awk -f ~/cifs-lab/Hafiz/scripts/hafiz_scripts/bam_tlen_filter2.awk | samtools view -Sb - > cmt3_rep2_methyl_MNase_sorted_dedup_sn_f.bam
bamToBed -bedpe -i cmt3_rep2_methyl_MNase_sorted_dedup_sn_f.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > cmt3_rep2_methyl_MNase_sorted_dedup_sn_f.bed

samtools view -h drm2_rep2_methyl_MNase_sorted_dedup_sn.bam | awk -f ~/cifs-lab/Hafiz/scripts/hafiz_scripts/bam_tlen_filter2.awk | samtools view -Sb - > drm2_rep2_methyl_MNase_sorted_dedup_sn_f.bam
bamToBed -bedpe -i drm2_rep2_methyl_MNase_sorted_dedup_sn_f.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > drm2_rep2_methyl_MNase_sorted_dedup_sn_f.bed

samtools view -h met13_rep2_methyl_MNase_sorted_dedup_sn.bam | awk -f ~/cifs-lab/Hafiz/scripts/hafiz_scripts/bam_tlen_filter2.awk | samtools view -Sb - > met13_rep2_methyl_MNase_sorted_dedup_sn_f.bam
bamToBed -bedpe -i met13_rep2_methyl_MNase_sorted_dedup_sn_f.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > met13_rep2_methyl_MNase_sorted_dedup_sn_f.bed


awk '{if($5>=2) print $0}' col0_rep2_methyl_MNase_sorted_dedup_sn_f.bed > col0_rep2_methyl_MNase_sorted_dedup_sn_f_mapq.bed
awk '{if($5>=2) print $0}' cmt2_rep2_methyl_MNase_sorted_dedup_sn_f.bed > cmt2_rep2_methyl_MNase_sorted_dedup_sn_f_mapq.bed
awk '{if($5>=2) print $0}' cmt3_rep2_methyl_MNase_sorted_dedup_sn_f.bed > cmt3_rep2_methyl_MNase_sorted_dedup_sn_f_mapq.bed
awk '{if($5>=2) print $0}' drm2_rep2_methyl_MNase_sorted_dedup_sn_f.bed > drm2_rep2_methyl_MNase_sorted_dedup_sn_f_mapq.bed
awk '{if($5>=2) print $0}' met13_rep2_methyl_MNase_sorted_dedup_sn_f.bed > met13_rep2_methyl_MNase_sorted_dedup_sn_f_mapq.bed


