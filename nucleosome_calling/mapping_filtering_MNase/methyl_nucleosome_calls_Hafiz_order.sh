#!/bin/sh

cd ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/meth_MNase_Hafiz_style_read_processing/


#bowtie2 -x ~/cifs-lab/Hafiz/tair10_genome_26/at26 -p 10 --very-sensitive --no-mixed --no-discordant --no-overlap --no-contain --fr -1 ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/trimmed_reads/112997_CGATGT_S1_L004_R1_001_val_1.fq -2 ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/trimmed_reads/112997_CGATGT_S1_L004_R2_001_val_2.fq -S col0_lo_h3_hafizMap.sam

bowtie2 -x ~/cifs-lab/Hafiz/tair10_genome_26/at26 -p 6 --very-sensitive --no-mixed --no-discordant --no-overlap --no-contain --fr -1 ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/trimmed_reads/112998_TTAGGC_S2_L004_R1_001_val_1.fq -2 ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/trimmed_reads/112998_TTAGGC_S2_L004_R2_001_val_2.fq -S met13_lo_h3_hafizMap.sam

bowtie2 -x ~/cifs-lab/Hafiz/tair10_genome_26/at26 -p 6 --very-sensitive --no-mixed --no-discordant --no-overlap --no-contain --fr -1 ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/trimmed_reads/112999_GCCAAT_S3_L004_R1_001_val_1.fq -2 ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/trimmed_reads/112999_GCCAAT_S3_L004_R2_001_val_2.fq -S cmt2_lo_h3_hafizMap.sam

bowtie2 -x ~/cifs-lab/Hafiz/tair10_genome_26/at26 -p 6 --very-sensitive --no-mixed --no-discordant --no-overlap --no-contain --fr -1 ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/trimmed_reads/113000_GATCAG_S4_L004_R1_001_val_1.fq -2 ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/trimmed_reads/113000_GATCAG_S4_L004_R2_001_val_2.fq -S cmt3_lo_h3_hafizMap.sam

bowtie2 -x ~/cifs-lab/Hafiz/tair10_genome_26/at26 -p 6 --very-sensitive --no-mixed --no-discordant --no-overlap --no-contain --fr -1 ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/trimmed_reads/113001_ATGTCA_S5_L004_R1_001_val_1.fq -2 ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/trimmed_reads/113001_ATGTCA_S5_L004_R2_001_val_2.fq -S drm2_lo_h3_hafizMap.sam


java -jar $PICARD SortSam INPUT=col0_lo_h3_hafizMap.sam OUTPUT=col0_lo_h3_hafizMap.bam SORT_ORDER=coordinate ### sam to bam
java -jar $PICARD SortSam INPUT=met13_lo_h3_hafizMap.sam OUTPUT=met13_lo_h3_hafizMap.bam SORT_ORDER=coordinate ### sam to bam
java -jar $PICARD SortSam INPUT=cmt2_lo_h3_hafizMap.sam OUTPUT=cmt2_lo_h3_hafizMap.bam SORT_ORDER=coordinate ### sam to bam
java -jar $PICARD SortSam INPUT=cmt3_lo_h3_hafizMap.sam OUTPUT=cmt3_lo_h3_hafizMap.bam SORT_ORDER=coordinate ### sam to bam
java -jar $PICARD SortSam INPUT=drm2_lo_h3_hafizMap.sam OUTPUT=drm2_lo_h3_hafizMap.bam SORT_ORDER=coordinate ### sam to bam

java -jar $PICARD MarkDuplicates REMOVE_DUPLICATES=true I=col0_lo_h3_hafizMap.bam O=col0_lo_h3_hafizMap_dedup.bam M=col0_h3_lo_dup_metrics.txt ### remove duplicates
java -jar $PICARD BuildBamIndex INPUT=col0_lo_h3_hafizMap_dedup.bam ### Index dedup bam file
java -jar $PICARD MarkDuplicates REMOVE_DUPLICATES=true I=met13_lo_h3_hafizMap.bam O=met13_lo_h3_hafizMap_dedup.bam M=met13_h3_lo_dup_metrics.txt ### remove duplicates
java -jar $PICARD BuildBamIndex INPUT=met13_lo_h3_hafizMap_dedup.bam ### Index dedup bam file
java -jar $PICARD MarkDuplicates REMOVE_DUPLICATES=true I=cmt2_lo_h3_hafizMap.bam O=cmt2_lo_h3_hafizMap_dedup.bam M=cmt2_h3_lo_dup_metrics.txt ### remove duplicates
java -jar $PICARD BuildBamIndex INPUT=cmt2_lo_h3_hafizMap_dedup.bam ### Index dedup bam file
java -jar $PICARD MarkDuplicates REMOVE_DUPLICATES=true I=cmt3_lo_h3_hafizMap.bam O=cmt3_lo_h3_hafizMap_dedup.bam M=cmt3_h3_lo_dup_metrics.txt ### remove duplicates
java -jar $PICARD BuildBamIndex INPUT=cmt3_lo_h3_hafizMap_dedup.bam ### Index dedup bam file
java -jar $PICARD MarkDuplicates REMOVE_DUPLICATES=true I=drm2_lo_h3_hafizMap.bam O=drm2_lo_h3_hafizMap_dedup.bam M=drm2_h3_lo_dup_metrics.txt ### remove duplicates
java -jar $PICARD BuildBamIndex INPUT=drm2_lo_h3_hafizMap_dedup.bam ### Index dedup bam file

samtools sort -n col0_lo_h3_hafizMap_dedup.bam -o col0_lo_h3_hafizMap_dedup_sn.bam ### sort dedup file
samtools sort -n met13_lo_h3_hafizMap_dedup.bam -o met13_lo_h3_hafizMap_dedup_sn.bam ### sort dedup file
samtools sort -n cmt2_lo_h3_hafizMap_dedup.bam -o cmt2_lo_h3_hafizMap_dedup_sn.bam ### sort dedup file
samtools sort -n cmt3_lo_h3_hafizMap_dedup.bam -o cmt3_lo_h3_hafizMap_dedup_sn.bam ### sort dedup file
samtools sort -n drm2_lo_h3_hafizMap_dedup.bam -o drm2_lo_h3_hafizMap_dedup_sn.bam ### sort dedup file

samtools view -h col0_lo_h3_hafizMap_dedup_sn.bam | awk -f ~/cifs-lab/Hafiz/scripts/hafiz_scripts/bam_tlen_filter2.awk | samtools view -Sb - > col0_lo_h3_hafizMap_f.bam ### filter fragment size
bamToBed -bedpe -i col0_lo_h3_hafizMap_f.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > col0_lo_h3_hafizMap_f.bed ### bam to bed : mapped-dedup-sorted-fragSelect file
samtools view -h met13_lo_h3_hafizMap_dedup_sn.bam | awk -f ~/cifs-lab/Hafiz/scripts/hafiz_scripts/bam_tlen_filter2.awk | samtools view -Sb - > met13_lo_h3_hafizMap_f.bam ### filter fragment size
bamToBed -bedpe -i met13_lo_h3_hafizMap_f.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > met13_lo_h3_hafizMap_f.bed ### bam to bed : mapped-dedup-sorted-fragSelect file
samtools view -h cmt2_lo_h3_hafizMap_dedup_sn.bam | awk -f ~/cifs-lab/Hafiz/scripts/hafiz_scripts/bam_tlen_filter2.awk | samtools view -Sb - > cmt2_lo_h3_hafizMap_f.bam ### filter fragment size
bamToBed -bedpe -i cmt2_lo_h3_hafizMap_f.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > cmt2_lo_h3_hafizMap_f.bed ### bam to bed : mapped-dedup-sorted-fragSelect file
samtools view -h cmt3_lo_h3_hafizMap_dedup_sn.bam | awk -f ~/cifs-lab/Hafiz/scripts/hafiz_scripts/bam_tlen_filter2.awk | samtools view -Sb - > cmt3_lo_h3_hafizMap_f.bam ### filter fragment size
bamToBed -bedpe -i cmt3_lo_h3_hafizMap_f.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > cmt3_lo_h3_hafizMap_f.bed ### bam to bed : mapped-dedup-sorted-fragSelect file
samtools view -h drm2_lo_h3_hafizMap_dedup_sn.bam | awk -f ~/cifs-lab/Hafiz/scripts/hafiz_scripts/bam_tlen_filter2.awk | samtools view -Sb - > drm2_lo_h3_hafizMap_f.bam ### filter fragment size
bamToBed -bedpe -i drm2_lo_h3_hafizMap_f.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > drm2_lo_h3_hafizMap_f.bed ### bam to bed : mapped-dedup-sorted-fragSelect file

#bamToBed -i col0_lo_h3_hafizMap_f.bam> col0_lo_h3_hafizMap_danpos_i.bed ### non PE bam to bed

bamToBed -bedpe -i col0_lo_h3_hafizMap_dedup_sn.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > col0_lo_h3_hafizMap.bed
bamToBed -bedpe -i met13_lo_h3_hafizMap_dedup_sn.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > met13_lo_h3_hafizMap.bed
bamToBed -bedpe -i cmt2_lo_h3_hafizMap_dedup_sn.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > cmt2_lo_h3_hafizMap.bed
bamToBed -bedpe -i cmt3_lo_h3_hafizMap_dedup_sn.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > cmt3_lo_h3_hafizMap.bed
bamToBed -bedpe -i drm2_lo_h3_hafizMap_dedup_sn.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > drm2_lo_h3_hafizMap.bed



awk '{if($5>=2) print $0}' col0_lo_h3_hafizMap.bed > col0_lo_h3_hafizMap_mapq.bed ### filter for mapq>2
awk '{if($5>=2) print $0}' met13_lo_h3_hafizMap.bed > met13_lo_h3_hafizMap_mapq.bed ### filter for mapq>2
awk '{if($5>=2) print $0}' cmt2_lo_h3_hafizMap.bed > cmt2_lo_h3_hafizMap_mapq.bed ### filter for mapq>2
awk '{if($5>=2) print $0}' cmt3_lo_h3_hafizMap.bed > cmt3_lo_h3_hafizMap_mapq.bed ### filter for mapq>2
awk '{if($5>=2) print $0}' drm2_lo_h3_hafizMap.bed > drm2_lo_h3_hafizMap_mapq.bed ### filter for mapq>2


