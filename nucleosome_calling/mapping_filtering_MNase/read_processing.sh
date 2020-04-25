cd /home/shriyas/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/reads

gunzip -c ~/cifs-lab/sequencing_dataset_archive/meth_mutants_MNase/ruddle.brcf.med.umich.edu/Run_2399/wierzbicki/Sample_112997/112997_CGATGT_S1_L004_R1_001.fastq.gz > ./112997_CGATGT_S1_L004_R1_001.fastq
gunzip -c ~/cifs-lab/sequencing_dataset_archive/meth_mutants_MNase/ruddle.brcf.med.umich.edu/Run_2399/wierzbicki/Sample_112997/112997_CGATGT_S1_L004_R2_001.fastq.gz > ./112997_CGATGT_S1_L004_R2_001.fastq

gunzip -c ~/cifs-lab/sequencing_dataset_archive/meth_mutants_MNase/ruddle.brcf.med.umich.edu/Run_2399/wierzbicki/Sample_112998/112998_TTAGGC_S2_L004_R1_001.fastq.gz > 112998_TTAGGC_S2_L004_R1_001.fastq
gunzip -c ~/cifs-lab/sequencing_dataset_archive/meth_mutants_MNase/ruddle.brcf.med.umich.edu/Run_2399/wierzbicki/Sample_112998/112998_TTAGGC_S2_L004_R2_001.fastq.gz > 112998_TTAGGC_S2_L004_R2_001.fastq

gunzip -c ~/cifs-lab/sequencing_dataset_archive/meth_mutants_MNase/ruddle.brcf.med.umich.edu/Run_2399/wierzbicki/Sample_112999/112999_GCCAAT_S3_L004_R1_001.fastq.gz > 112999_GCCAAT_S3_L004_R1_001.fastq
gunzip -c ~/cifs-lab/sequencing_dataset_archive/meth_mutants_MNase/ruddle.brcf.med.umich.edu/Run_2399/wierzbicki/Sample_112999/112999_GCCAAT_S3_L004_R2_001.fastq.gz > 112999_GCCAAT_S3_L004_R2_001.fastq

gunzip -c ~/cifs-lab/sequencing_dataset_archive/meth_mutants_MNase/ruddle.brcf.med.umich.edu/Run_2399/wierzbicki/Sample_113000/113000_GATCAG_S4_L004_R1_001.fastq.gz > 113000_GATCAG_S4_L004_R1_001.fastq
gunzip -c ~/cifs-lab/sequencing_dataset_archive/meth_mutants_MNase/ruddle.brcf.med.umich.edu/Run_2399/wierzbicki/Sample_113000/113000_GATCAG_S4_L004_R2_001.fastq.gz > 113000_GATCAG_S4_L004_R2_001.fastq

gunzip -c ~/cifs-lab/sequencing_dataset_archive/meth_mutants_MNase/ruddle.brcf.med.umich.edu/Run_2399/wierzbicki/Sample_113001/113001_ATGTCA_S5_L004_R1_001.fastq.gz > 113001_ATGTCA_S5_L004_R1_001.fastq
gunzip -c ~/cifs-lab/sequencing_dataset_archive/meth_mutants_MNase/ruddle.brcf.med.umich.edu/Run_2399/wierzbicki/Sample_113001/113001_ATGTCA_S5_L004_R2_001.fastq.gz > 113001_ATGTCA_S5_L004_R2_001.fastq

fastqc --casava *fastq

cd ../
mkdir trimmed_reads
cd trimmed_reads

trim_galore --length 10 --fastqc --paired ../reads/112997_CGATGT_S1_L004_R1_001.fastq ../reads/112997_CGATGT_S1_L004_R2_001.fastq
trim_galore --length 35 --fastqc --paired ../reads/112998_TTAGGC_S2_L004_R1_001.fastq ../reads/112998_TTAGGC_S2_L004_R2_001.fastq
trim_galore --length 35 --fastqc --paired ../reads/112999_GCCAAT_S3_L004_R1_001.fastq ../reads/112999_GCCAAT_S3_L004_R2_001.fastq

trim_galore --length 35 --fastqc --paired ../reads/113000_GATCAG_S4_L004_R1_001.fastq ../reads/113000_GATCAG_S4_L004_R2_001.fastq
trim_galore --length 35 --fastqc --paired ../reads/113001_ATGTCA_S5_L004_R1_001.fastq ../reads/113001_ATGTCA_S5_L004_R2_001.fastq

cd ../
mkdir mapped_reads
cd mapped_reads/

bowtie2 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -p 6 -q --no-mixed --no-discordant --no-overlap --no-contain --fr -1 ../trimmed_reads/112997_CGATGT_S1_L004_R1_001_val_1.fq -2 ../trimmed_reads/112997_CGATGT_S1_L004_R2_001_val_2.fq -S col0_methyl_MNase.sam

75655748 reads; of these:
  75655748 (100.00%) were paired; of these:
    11484489 (15.18%) aligned concordantly 0 times
    36799480 (48.64%) aligned concordantly exactly 1 time
    27371779 (36.18%) aligned concordantly >1 times
84.82% overall alignment rate

bowtie2 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -p 6 -q --no-mixed --no-discordant --no-overlap --no-contain --fr -1 ../trimmed_reads/112998_TTAGGC_S2_L004_R1_001_val_1.fq -2 ../trimmed_reads/112998_TTAGGC_S2_L004_R2_001_val_2.fq -S met13_methyl_MNase.sam

52972840 reads; of these:
  52972840 (100.00%) were paired; of these:
    11010437 (20.79%) aligned concordantly 0 times
    24803503 (46.82%) aligned concordantly exactly 1 time
    17158900 (32.39%) aligned concordantly >1 times
79.21% overall alignment rate


bowtie2 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -p 6 -q --no-mixed --no-discordant --no-overlap --no-contain --fr -1 ../trimmed_reads/112999_GCCAAT_S3_L004_R1_001_val_1.fq -2 ../trimmed_reads/112999_GCCAAT_S3_L004_R2_001_val_2.fq -S cmt2_methyl_MNase.sam

61023254 reads; of these:
  61023254 (100.00%) were paired; of these:
    10421213 (17.08%) aligned concordantly 0 times
    29699501 (48.67%) aligned concordantly exactly 1 time
    20902540 (34.25%) aligned concordantly >1 times
82.92% overall alignment rate


bowtie2 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -p 6 -q --no-mixed --no-discordant --no-overlap --no-contain --fr -1 ../trimmed_reads/113000_GATCAG_S4_L004_R1_001_val_1.fq -2 ../trimmed_reads/113000_GATCAG_S4_L004_R2_001_val_2.fq -S cmt3_methyl_MNase.sam

64337144 reads; of these:
  64337144 (100.00%) were paired; of these:
    11186142 (17.39%) aligned concordantly 0 times
    30925112 (48.07%) aligned concordantly exactly 1 time
    22225890 (34.55%) aligned concordantly >1 times
82.61% overall alignment rate


bowtie2 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -p 6 -q --no-mixed --no-discordant --no-overlap --no-contain --fr -1 ../trimmed_reads/113001_ATGTCA_S5_L004_R1_001_val_1.fq -2 ../trimmed_reads/113001_ATGTCA_S5_L004_R2_001_val_2.fq -S drm2_methyl_MNase.sam

51359788 reads; of these:
  51359788 (100.00%) were paired; of these:
    14022396 (27.30%) aligned concordantly 0 times
    10467217 (20.38%) aligned concordantly exactly 1 time
    26870175 (52.32%) aligned concordantly >1 times
72.70% overall alignment rate


java -jar $PICARD SortSam INPUT=../col0_methyl_MNase.sam OUTPUT=col0_methyl_MNase_sorted.bam SORT_ORDER=coordinate
java -jar $PICARD MarkDuplicates REMOVE_DUPLICATES=true I=col0_methyl_MNase_sorted.bam O=col0_methyl_MNase_sorted_dedup.bam M=col0_h3_lo_dup_metrics.txt
java -jar $PICARD BuildBamIndex INPUT=col0_methyl_MNase_sorted_dedup.bam
samtools sort -n col0_methyl_MNase_sorted_dedup.bam col0_methyl_MNase_sorted_dedup_sn
bamToBed -bedpe -i col0_methyl_MNase_sorted_dedup_sn.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > col0_methyl_MNase_sorted_dedup.bed

java -jar $PICARD SortSam INPUT=../met13_methyl_MNase.sam OUTPUT=met13_methyl_MNase_sorted.bam SORT_ORDER=coordinate
java -jar $PICARD MarkDuplicates REMOVE_DUPLICATES=true I=met13_methyl_MNase_sorted.bam O=met13_methyl_MNase_sorted_dedup.bam M=met13_h3_lo_dup_metrics.txt
java -jar $PICARD BuildBamIndex INPUT=met13_methyl_MNase_sorted_dedup.bam
samtools sort -n met13_methyl_MNase_sorted_dedup.bam met13_methyl_MNase_sorted_dedup_sn
bamToBed -bedpe -i met13_methyl_MNase_sorted_dedup_sn.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > met13_methyl_MNase_sorted_dedup.bed

java -jar $PICARD SortSam INPUT=../cmt2_methyl_MNase.sam OUTPUT=cmt2_methyl_MNase_sorted.bam SORT_ORDER=coordinate
java -jar $PICARD MarkDuplicates REMOVE_DUPLICATES=true I=cmt2_methyl_MNase_sorted.bam O=cmt2_methyl_MNase_sorted_dedup.bam M=cmt2_h3_lo_dup_metrics.txt
java -jar $PICARD BuildBamIndex INPUT=cmt2_methyl_MNase_sorted_dedup.bam
samtools sort -n cmt2_methyl_MNase_sorted_dedup.bam cmt2_methyl_MNase_sorted_dedup_sn
bamToBed -bedpe -i cmt2_methyl_MNase_sorted_dedup_sn.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > cmt2_methyl_MNase_sorted_dedup.bed

java -jar $PICARD SortSam INPUT=../cmt3_methyl_MNase.sam OUTPUT=cmt3_methyl_MNase_sorted.bam SORT_ORDER=coordinate
java -jar $PICARD MarkDuplicates REMOVE_DUPLICATES=true I=cmt3_methyl_MNase_sorted.bam O=cmt3_methyl_MNase_sorted_dedup.bam M=cmt3_h3_lo_dup_metrics.txt
java -jar $PICARD BuildBamIndex INPUT=cmt3_methyl_MNase_sorted_dedup.bam
samtools sort -n cmt3_methyl_MNase_sorted_dedup.bam cmt3_methyl_MNase_sorted_dedup_sn
bamToBed -bedpe -i cmt3_methyl_MNase_sorted_dedup_sn.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > cmt3_methyl_MNase_sorted_dedup.bed

java -jar $PICARD SortSam INPUT=../drm2_methyl_MNase.sam OUTPUT=drm2_methyl_MNase_sorted.bam SORT_ORDER=coordinate
java -jar $PICARD MarkDuplicates REMOVE_DUPLICATES=true I=drm2_methyl_MNase_sorted.bam O=drm2_methyl_MNase_sorted_dedup.bam M=drm2_h3_lo_dup_metrics.txt
java -jar $PICARD BuildBamIndex INPUT=drm2_methyl_MNase_sorted_dedup.bam
samtools sort -n drm2_methyl_MNase_sorted_dedup.bam drm2_methyl_MNase_sorted_dedup_sn
bamToBed -bedpe -i drm2_methyl_MNase_sorted_dedup_sn.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > drm2_methyl_MNase_sorted_dedup.bed


samtools view -bq 2 col0_methyl_MNase_sorted_dedup_sn.bam > col0_methyl_MNase_sorted_dedup_sn_mapq2.bam
samtools view -bq 2 drm2_methyl_MNase_sorted_dedup_sn.bam > drm2_methyl_MNase_sorted_dedup_sn_mapq2.bam


awk '{if($5>=2) print $0}' col0_methyl_MNase_sorted_dedup.bed > col0_methyl_MNase_sorted_dedup_mapq.bed ### filter for mapq>2
awk '{if($5>=2) print $0}' met13_methyl_MNase_sorted_dedup.bed > met13_methyl_MNase_sorted_dedup_mapq.bed ### filter for mapq>2
awk '{if($5>=2) print $0}' cmt2_methyl_MNase_sorted_dedup.bed > cmt2_methyl_MNase_sorted_dedup_mapq.bed ### filter for mapq>2
awk '{if($5>=2) print $0}' cmt3_methyl_MNase_sorted_dedup.bed > cmt3_methyl_MNase_sorted_dedup_mapq.bed ### filter for mapq>2
awk '{if($5>=2) print $0}' drm2_methyl_MNase_sorted_dedup.bed > drm2_methyl_MNase_sorted_dedup_mapq.bed ### filter for mapq>2





samtools view -bS -o col0_methyl_MNase.bam col0_methyl_MNase.sam
samtools view -bS -o met13_methyl_MNase.bam met13_methyl_MNase.sam
samtools view -bS -o cmt2_methyl_MNase.bam cmt2_methyl_MNase.sam
samtools view -bS -o cmt3_methyl_MNase.bam cmt3_methyl_MNase.sam
samtools view -bS -o drm2_methyl_MNase.bam drm2_methyl_MNase.sam

bamToBed -bedpe -i col0_methyl_MNase.bam | awk '{if($2>0) print $0}' | sort -k1 -k2 -n | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > col0_methyl_MNase.bed
bamToBed -bedpe -i met13_methyl_MNase.bam | awk '{if($2>0) print $0}' | sort -k1 -k2 -n | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > met13_methyl_MNase.bed
bamToBed -bedpe -i cmt2_methyl_MNase.bam | awk '{if($2>0) print $0}' | sort -k1 -k2 -n | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > cmt2_methyl_MNase.bed
bamToBed -bedpe -i cmt3_methyl_MNase.bam | awk '{if($2>0) print $0}' | sort -k1 -k2 -n | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > cmt3_methyl_MNase.bed
bamToBed -bedpe -i drm2_methyl_MNase.bam | awk '{if($2>0) print $0}' | sort -k1 -k2 -n | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > drm2_methyl_MNase.bed



samtools view -h col0_methyl_MNase_sorted_dedup_sn.bam | awk -f ~/cifs-lab/Hafiz/scripts/hafiz_scripts/bam_tlen_filter2.awk | samtools view -Sb - > col0_methyl_MNase_sorted_dedup_sn_f.bam
bamToBed -bedpe -i col0_methyl_MNase_sorted_dedup_sn_f.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > col0_methyl_MNase_sorted_dedup_sn_f.bed

samtools view -h cmt2_methyl_MNase_sorted_dedup_sn.bam | awk -f ~/cifs-lab/Hafiz/scripts/hafiz_scripts/bam_tlen_filter2.awk | samtools view -Sb - > cmt2_methyl_MNase_sorted_dedup_sn_f.bam
bamToBed -bedpe -i cmt2_methyl_MNase_sorted_dedup_sn_f.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > cmt2_methyl_MNase_sorted_dedup_sn_f.bed

samtools view -h cmt3_methyl_MNase_sorted_dedup_sn.bam | awk -f ~/cifs-lab/Hafiz/scripts/hafiz_scripts/bam_tlen_filter2.awk | samtools view -Sb - > cmt3_methyl_MNase_sorted_dedup_sn_f.bam
bamToBed -bedpe -i cmt3_methyl_MNase_sorted_dedup_sn_f.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > cmt3_methyl_MNase_sorted_dedup_sn_f.bed

samtools view -h drm2_methyl_MNase_sorted_dedup_sn.bam | awk -f ~/cifs-lab/Hafiz/scripts/hafiz_scripts/bam_tlen_filter2.awk | samtools view -Sb - > drm2_methyl_MNase_sorted_dedup_sn_f.bam
bamToBed -bedpe -i drm2_methyl_MNase_sorted_dedup_sn_f.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > drm2_methyl_MNase_sorted_dedup_sn_f.bed

samtools view -h met13_methyl_MNase_sorted_dedup_sn.bam | awk -f ~/cifs-lab/Hafiz/scripts/hafiz_scripts/bam_tlen_filter2.awk | samtools view -Sb - > met13_methyl_MNase_sorted_dedup_sn_f.bam
bamToBed -bedpe -i met13_methyl_MNase_sorted_dedup_sn_f.bam | awk '{if($2>0) print $0}' | sortBed -i stdin  | awk '{print $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }' > met13_methyl_MNase_sorted_dedup_sn_f.bed


awk '{if($5>=2) print $0}' col0_methyl_MNase_sorted_dedup_sn_f.bed > col0_methyl_MNase_sorted_dedup_sn_f_mapq.bed
awk '{if($5>=2) print $0}' cmt2_methyl_MNase_sorted_dedup_sn_f.bed > cmt2_methyl_MNase_sorted_dedup_sn_f_mapq.bed
awk '{if($5>=2) print $0}' cmt3_methyl_MNase_sorted_dedup_sn_f.bed > cmt3_methyl_MNase_sorted_dedup_sn_f_mapq.bed
awk '{if($5>=2) print $0}' drm2_methyl_MNase_sorted_dedup_sn_f.bed > drm2_methyl_MNase_sorted_dedup_sn_f_mapq.bed
awk '{if($5>=2) print $0}' met13_methyl_MNase_sorted_dedup_sn_f.bed > met13_methyl_MNase_sorted_dedup_sn_f_mapq.bed

