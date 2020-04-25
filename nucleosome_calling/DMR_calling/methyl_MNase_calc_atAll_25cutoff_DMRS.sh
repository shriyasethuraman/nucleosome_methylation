#!/bin/sh

############### DRM2 ###############

#cd ~/cifs-lab/Shriya/Samples/manuscript_prep/nucleosome_methylation/working/methyl_mutants_DMRs/All_final_DMRs/drm2_DMRs/
#sed 's/"//g' ../mydiff_reg_drm2_col0_d25_q01.txt | sed 's/chr//g' | awk '$1==1 || $1==2 || $1==3 || $1==4 || $1==5' > mydiff_reg_drm2_col0_d25_q01.bed
#cd ~/cifs-lab/Shriya/Samples/manuscript_prep/nucleosome_methylation/working/methyl_mutants_DMRs/All_final_DMRs/cmt3_DMRs
#sed 's/"//g' ../mydiff_CHG_reg_cmt3_col0_d25_q01.txt | sed 's/chr//g' | awk '$1==1 || $1==2 || $1==3 || $1==4 || $1==5' > mydiff_CHG_reg_cmt3_col0_d25_q01.bed
#cd ~/cifs-lab/Shriya/Samples/manuscript_prep/nucleosome_methylation/working/methyl_mutants_DMRs/All_final_DMRs/cmt2_DMRs
#sed 's/"//g' ../mydiff_reg_cmt2_col0_d25_q01.txt | sed 's/chr//g' | awk '$1==1 || $1==2 || $1==3 || $1==4 || $1==5' > mydiff_reg_cmt2_col0_d25_q01.bed
#cd ~/cifs-lab/Shriya/Samples/manuscript_prep/nucleosome_methylation/working/methyl_mutants_DMRs/All_final_DMRs/met1_DMRs
#sed 's/"//g' ../mydiff_CG_reg_met1_col0_d25_q01.txt | sed 's/chr//g' | awk '$1==1 || $1==2 || $1==3 || $1==4 || $1==5' > mydiff_CG_reg_met1_col0_d25_q01.bed


### CHH-DMRs #######
cd ~/cifs-lab/Shriya/Samples/manuscript_prep/nucleosome_methylation/working/methyl_mutants_DMRs/All_final_DMRs/drm2_DMRs/

#awk '{OFS="\t"} $2>500 {print $1,$2-500,$3+500}' mydiff_reg_drm2_col0_d25_q01.bed | coverageBed -d -a stdin -b ~/prog_run/mock_mapped.bed | awk '{OFS="\t"}{print $1,$2+$4-1,$2+$4, $2,$3,$4}' > nucleotideWise_drm2_CHH_DMR_d25_500bp.bed


#intersectBed -c -a nucleotideWise_drm2_CHH_DMR_d25_500bp.bed -b ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/col0_methyl_MNase_sorted_dedup_sn_f_mapq.bed > col0_MNase_drm2_CHH_DMR_counts.bed
#intersectBed -c -a nucleotideWise_drm2_CHH_DMR_d25_500bp.bed -b ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/drm2_methyl_MNase_sorted_dedup_sn_f_mapq.bed > drm2_MNase_drm2_CHH_DMR_counts.bed

#awk '{OFS="\t"}{if($6==1) s+=1; print $0,s}' col0_MNase_drm2_CHH_DMR_counts.bed > col0_MNase_drm2_CHH_DMR_counts_2.bed
#awk '{OFS="\t"}{if($6==1) s+=1; print $0,s}' drm2_MNase_drm2_CHH_DMR_counts.bed > drm2_MNase_drm2_CHH_DMR_counts_2.bed

#perl ~/cifs-lab/Shriya/Samples/manuscript_prep/nucleosome_methylation/working/methyl_mutants_DMRs/convert_2d_2_MNase_2.pl drm2_MNase_drm2_CHH_DMR_counts_2.bed 6 0 > drm2_MNase_drm2_CHH_DMR_counts_2.csv
#perl ~/cifs-lab/Shriya/Samples/manuscript_prep/nucleosome_methylation/working/methyl_mutants_DMRs/convert_2d_2_MNase_2.pl col0_MNase_drm2_CHH_DMR_counts_2.bed 6 0 > col0_MNase_drm2_CHH_DMR_counts_2.csv

#wc -l ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/col0_methyl_MNase_sorted_dedup_sn_f_mapq.bed 
#987044 /home/shriyas/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/col0_methyl_MNase_sorted_dedup_sn_f_mapq.bed

#wc -l ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/drm2_methyl_MNase_sorted_dedup_sn_f_mapq.bed 
#245458 /home/shriyas/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/drm2_methyl_MNase_sorted_dedup_sn_f_mapq.bed

intersectBed -c -a nucleotideWise_drm2_CHH_DMR_d25_500bp.bed -b ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/rep2_oct2019/mapped_reads/process_map_reads/col0_rep2_methyl_MNase_sorted_dedup_sn_f_mapq.bed > col0_rep2_MNase_drm2_d25_CHH_DMR_counts.bed
intersectBed -c -a nucleotideWise_drm2_CHH_DMR_d25_500bp.bed -b ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/rep2_oct2019/mapped_reads/process_map_reads/drm2_rep2_methyl_MNase_sorted_dedup_sn_f_mapq.bed > drm2_rep2_MNase_drm2_d25_CHH_DMR_counts.bed

awk '{OFS="\t"}{if($6==1) s+=1; print $0,s}' col0_rep2_MNase_drm2_d25_CHH_DMR_counts.bed > col0_rep2_MNase_drm2_d25_CHH_DMR_counts_2.bed
awk '{OFS="\t"}{if($6==1) s+=1; print $0,s}' drm2_rep2_MNase_drm2_d25_CHH_DMR_counts.bed > drm2_rep2_MNase_drm2_d25_CHH_DMR_counts_2.bed

perl ../../convert_2d_2_MNase_2.pl drm2_rep2_MNase_drm2_d25_CHH_DMR_counts_2.bed 6 0 > drm2_rep2_MNase_drm2_d25_CHH_DMR_counts_2.csv
perl ../../convert_2d_2_MNase_2.pl col0_rep2_MNase_drm2_d25_CHH_DMR_counts_2.bed 6 0 > col0_rep2_MNase_drm2_d25_CHH_DMR_counts_2.csv


#wc -l ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/rep2_oct2019/mapped_reads/process_map_reads/drm2_rep2_methyl_MNase_sorted_dedup_sn_f_mapq.bed
#3115678 /home/shriyas/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/rep2_oct2019/mapped_reads/process_map_reads/drm2_rep2_methyl_MNase_sorted_dedup_sn_f_mapq.bed


#intersectBed -wao -a nucleotideWise_drm2_CHH_DMR_d25_500bp.bed -b ~/cifs-lab/Shriya/Samples/Jacobsen/Col0-merged/CHH_col0_merged_12.bed > drm2_CHH_DMR_col0_meth.bed
awk '{OFS="\t"}{if($9=="-1") print $1,$2,$3,$4,$5,$6,"NA"; else {print $1,$2,$3,$4,$5,$6,$10/$11}}' drm2_CHH_DMR_col0_meth.bed > drm2_CHH_DMR_col0_meth_temp2.bed
awk '{OFS="\t"}{if($6==1 && $1==1 || $1==2 || $1==3 || $1==4 || $1==5) s+=1; print $0, s}' drm2_CHH_DMR_col0_meth_temp2.bed > WT_CHH_dataset.bed


intersectBed -wao -a nucleotideWise_drm2_CHH_DMR_d25_500bp.bed -b ~/cifs-lab/Shriya/Samples/Jacobsen/Drm2/meth_extract/drm2_rep12_CHHmeth.bed > drm2_CHH_DMR_drm2_meth.bed
awk '{OFS="\t"}{if($9=="-1") print $1,$2,$3,$4,$5,$6,"NA"; else {print $1,$2,$3,$4,$5,$6,$10/$11}}' drm2_CHH_DMR_drm2_meth.bed > drm2_CHH_DMR_drm2_meth_temp2.bed
awk '{OFS="\t"}{if($6==1 && $1==1 || $1==2 || $1==3 || $1==4 || $1==5) s+=1; print $0, s}' drm2_CHH_DMR_drm2_meth_temp2.bed > drm2_CHH_dataset.bed

perl ~/cifs-lab/Shriya/Samples/manuscript_prep/nucleosome_methylation/working/methyl_mutants_DMRs/convert_2d_2_MNase_2.pl WT_CHH_dataset.bed 6 0 > WT_CHH_dataset.csv
perl ~/cifs-lab/Shriya/Samples/manuscript_prep/nucleosome_methylation/working/methyl_mutants_DMRs/convert_2d_2_MNase_2.pl drm2_CHH_dataset.bed 6 0 > drm2_CHH_dataset.csv



#intersectBed -c -a nucleotideWise_drm2_CG_DMR_d25_500bp.bed -b ~/cifs-lab/Shriya/Samples/RIP_dms3_drm2/bowtie2_mapping/old_RIP_test/AllCol.bed > col0_old_RIP_RCs_drm2_CG_DMR.bed
#intersectBed -c -a nucleotideWise_drm2_CHH_DMR_d25_500bp.bed -b ~/cifs-lab/Shriya/Samples/RIP_dms3_drm2/bowtie2_mapping/old_RIP_test/AllCol.bed > col0_old_RIP_RCs_drm2_CHH_DMR.bed


paste col0_MNase_drm2_CHH_DMR_counts_2.bed col0_rep2_MNase_drm2_d25_CHH_DMR_counts_2.bed | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7+$15,$8}' > col0_rep1+rep2_MNase_drm2_d25_CHH_DMR_counts.bed
paste drm2_MNase_drm2_CHH_DMR_counts_2.bed drm2_rep2_MNase_drm2_d25_CHH_DMR_counts_2.bed | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7+$15,$8}' > drm2_rep1+rep2_MNase_drm2_d25_CHH_DMR_counts.bed

cd supp_DMRs
paste col0_MNase_drm2_CG_DMR_counts_2.bed col0_rep2_MNase_drm2_d25_CG_DMR_counts_2.bed | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7+$15,$8}' > col0_rep1+rep2_MNase_drm2_d25_CG_DMR_counts.bed
paste drm2_MNase_drm2_CG_DMR_counts_2.bed drm2_rep2_MNase_drm2_d25_CG_DMR_counts_2.bed | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7+$15,$8}' > drm2_rep1+rep2_MNase_drm2_d25_CG_DMR_counts.bed
paste col0_MNase_drm2_CHG_DMR_counts_2.bed col0_rep2_MNase_drm2_d25_CHG_DMR_counts_2.bed | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7+$15,$8}' > col0_rep1+rep2_MNase_drm2_d25_CHG_DMR_counts.bed
paste drm2_MNase_drm2_CHG_DMR_counts_2.bed drm2_rep2_MNase_drm2_d25_CHG_DMR_counts_2.bed | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7+$15,$8}' > drm2_rep1+rep2_MNase_drm2_d25_CHG_DMR_counts.bed

############### met13 ###############

### CG-DMRs #######

cd ~/cifs-lab/Shriya/Samples/manuscript_prep/nucleosome_methylation/working/methyl_mutants_DMRs/All_final_DMRs/met1_DMRs

awk '{OFS="\t"} $2>500 {print $1,$2-500,$3+500}' mydiff_CG_reg_met1_col0_d25_q01.bed | coverageBed -d -a stdin -b ~/prog_run/mock_mapped.bed | awk '{OFS="\t"}{print $1,$2+$4-1,$2+$4, $2,$3,$4}' > nucleotideWise_met13_CG_DMR_d25_500bp.bed


intersectBed -c -a nucleotideWise_met13_CG_DMR_d25_500bp.bed -b ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/col0_methyl_MNase_sorted_dedup_sn_f_mapq.bed > col0_MNase_met13_CG_DMR_counts.bed
intersectBed -c -a nucleotideWise_met13_CG_DMR_d25_500bp.bed -b ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/met13_methyl_MNase_sorted_dedup_sn_f_mapq.bed > met13_MNase_met13_CG_DMR_counts.bed

awk '{OFS="\t"}{if($6==1) s+=1; print $0,s}' col0_MNase_met13_CG_DMR_counts.bed > col0_MNase_met13_CG_DMR_counts_2.bed
awk '{OFS="\t"}{if($6==1) s+=1; print $0,s}' met13_MNase_met13_CG_DMR_counts.bed > met13_MNase_met13_CG_DMR_counts_2.bed

perl ../../convert_2d_2_MNase_2.pl met13_MNase_met13_CG_DMR_counts_2.bed 6 0 > met13_MNase_met13_CG_DMR_counts_2.csv
perl ../../convert_2d_2_MNase_2.pl col0_MNase_met13_CG_DMR_counts_2.bed 6 0 > col0_MNase_met13_CG_DMR_counts_2.csv

#wc -l ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/col0_methyl_MNase_sorted_dedup_sn_f_mapq.bed 
#987044 /home/shriyas/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/col0_methyl_MNase_sorted_dedup_sn_f_mapq.bed

#wc -l ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/met13_methyl_MNase_sorted_dedup_sn_f_mapq.bed 
#675495 /home/shriyas/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/met13_methyl_MNase_sorted_dedup_sn_f_mapq.bed


intersectBed -wao -a nucleotideWise_met13_CG_DMR_d25_500bp.bed -b /home/shriyas/cifs-lab/Shriya/Samples/Jacobsen/met1/meth_extract/col0_met1_CGmeth.bed > met13_CG_DMR_col0_meth.bed
awk '{OFS="\t"}{if($9=="-1") print $1,$2,$3,$4,$5,$6,"NA"; else {print $1,$2,$3,$4,$5,$6,$10/$11}}' met13_CG_DMR_col0_meth.bed > met13_CG_DMR_col0_meth_temp2.bed
awk '{OFS="\t"}{if($6==1 && $1==1 || $1==2 || $1==3 || $1==4 || $1==5) s+=1; print $0, s}' met13_CG_DMR_col0_meth_temp2.bed > WT_CG_dataset.bed


intersectBed -wao -a nucleotideWise_met13_CG_DMR_d25_500bp.bed -b /home/shriyas/cifs-lab/Shriya/Samples/Jacobsen/met1/meth_extract/met1_CGmeth.bed > met13_CG_DMR_met13_meth.bed
awk '{OFS="\t"}{if($9=="-1") print $1,$2,$3,$4,$5,$6,"NA"; else {print $1,$2,$3,$4,$5,$6,$10/$11}}' met13_CG_DMR_met13_meth.bed > met13_CG_DMR_met13_meth_temp2.bed
awk '{OFS="\t"}{if($6==1 && $1==1 || $1==2 || $1==3 || $1==4 || $1==5) s+=1; print $0, s}' met13_CG_DMR_met13_meth_temp2.bed > met13_CG_dataset.bed

perl ~/cifs-lab/Shriya/Samples/manuscript_prep/nucleosome_methylation/working/methyl_mutants_DMRs/convert_2d_2_MNase_2.pl WT_CG_dataset.bed 6 0 > WT_CG_dataset.csv
perl ~/cifs-lab/Shriya/Samples/manuscript_prep/nucleosome_methylation/working/methyl_mutants_DMRs/convert_2d_2_MNase_2.pl met13_CG_dataset.bed 6 0 > met13_CG_dataset.csv


intersectBed -c -a nucleotideWise_met13_CG_DMR_d25_500bp.bed -b ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/rep2_oct2019/mapped_reads/process_map_reads/col0_rep2_methyl_MNase_sorted_dedup_sn_f_mapq.bed > col0_rep2_MNase_met13_d25_CG_DMR_counts.bed
intersectBed -c -a nucleotideWise_met13_CG_DMR_d25_500bp.bed -b ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/rep2_oct2019/mapped_reads/process_map_reads/met13_rep2_methyl_MNase_sorted_dedup_sn_f_mapq.bed > met13_rep2_MNase_met13_d25_CG_DMR_counts.bed

awk '{OFS="\t"}{if($6==1) s+=1; print $0,s}' col0_rep2_MNase_met13_d25_CG_DMR_counts.bed > col0_rep2_MNase_met13_d25_CG_DMR_counts_2.bed
awk '{OFS="\t"}{if($6==1) s+=1; print $0,s}' met13_rep2_MNase_met13_d25_CG_DMR_counts.bed > met13_rep2_MNase_met13_d25_CG_DMR_counts_2.bed

perl ../../convert_2d_2_MNase_2.pl met13_rep2_MNase_met13_d25_CG_DMR_counts_2.bed 6 0 > met13_rep2_MNase_met13_d25_CG_DMR_counts_2.csv
perl ../../convert_2d_2_MNase_2.pl col0_rep2_MNase_met13_d25_CG_DMR_counts_2.bed 6 0 > col0_rep2_MNase_met13_d25_CG_DMR_counts_2.csv

#wc -l ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/rep2_oct2019/mapped_reads/process_map_reads/col0_rep2_methyl_MNase_sorted_dedup_sn_f_mapq.bed
#1415365 /home/shriyas/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/rep2_oct2019/mapped_reads/process_map_reads/col0_rep2_methyl_MNase_sorted_dedup_sn_f_mapq.bed

#wc -l ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/rep2_oct2019/mapped_reads/process_map_reads/met13_rep2_methyl_MNase_sorted_dedup_sn_f_mapq.bed 
#568632 /home/shriyas/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/rep2_oct2019/mapped_reads/process_map_reads/met13_rep2_methyl_MNase_sorted_dedup_sn_f_mapq.bed


paste col0_MNase_met13_CG_DMR_counts_2.bed col0_rep2_MNase_met13_d25_CG_DMR_counts_2.bed | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7+$15,$8}' > col0_rep1+rep2_MNase_met13_d25_CG_DMR_counts.bed
paste met13_MNase_met13_CG_DMR_counts_2.bed met13_rep2_MNase_met13_d25_CG_DMR_counts_2.bed | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7+$15,$8}' > met13_rep1+rep2_MNase_met13_d25_CG_DMR_counts.bed

cd supp_DMRs
paste col0_MNase_met1_CHH_DMR_counts_2.bed col0_rep2_MNase_met1_d25_CHH_DMR_counts_2.bed | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7+$15,$8}' > col0_rep1+rep2_MNase_met13_d25_CHH_DMR_counts.bed
paste met1_MNase_met1_CHH_DMR_counts_2.bed met1_rep2_MNase_met1_d25_CHH_DMR_counts_2.bed | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7+$15,$8}' > met13_rep1+rep2_MNase_met13_d25_CHH_DMR_counts.bed
paste col0_MNase_met1_CHG_DMR_counts_2.bed col0_rep2_MNase_met1_d25_CHG_DMR_counts_2.bed | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7+$15,$8}' > col0_rep1+rep2_MNase_met13_d25_CHG_DMR_counts.bed
paste met1_MNase_met1_CHG_DMR_counts_2.bed met1_rep2_MNase_met1_d25_CHG_DMR_counts_2.bed | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7+$15,$8}' > met13_rep1+rep2_MNase_met13_d25_CHG_DMR_counts.bed



############### CMT2 ###############

### CHH-DMRs #######
cd ~/cifs-lab/Shriya/Samples/manuscript_prep/nucleosome_methylation/working/methyl_mutants_DMRs/All_final_DMRs/cmt2_DMRs


awk '{OFS="\t"} $2>500 {print $1,$2-500,$3+500}' mydiff_reg_cmt2_col0_d25_q01.bed | coverageBed -d -a stdin -b ~/prog_run/mock_mapped.bed | awk '{OFS="\t"}{print $1,$2+$4-1,$2+$4, $2,$3,$4}' > nucleotideWise_cmt2_CHH_DMR_d25_500bp.bed


intersectBed -c -a nucleotideWise_cmt2_CHH_DMR_d25_500bp.bed -b ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/col0_methyl_MNase_sorted_dedup_sn_f_mapq.bed > col0_MNase_cmt2_CHH_DMR_counts.bed
intersectBed -c -a nucleotideWise_cmt2_CHH_DMR_d25_500bp.bed -b ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/cmt2_methyl_MNase_sorted_dedup_sn_f_mapq.bed > cmt2_MNase_cmt2_CHH_DMR_counts.bed

awk '{OFS="\t"}{if($6==1) s+=1; print $0,s}' col0_MNase_cmt2_CHH_DMR_counts.bed > col0_MNase_cmt2_CHH_DMR_counts_2.bed
awk '{OFS="\t"}{if($6==1) s+=1; print $0,s}' cmt2_MNase_cmt2_CHH_DMR_counts.bed > cmt2_MNase_cmt2_CHH_DMR_counts_2.bed

perl ../../convert_2d_2_MNase_2.pl cmt2_MNase_cmt2_CHH_DMR_counts_2.bed 6 0 > cmt2_MNase_cmt2_CHH_DMR_counts_2.csv
perl ../../convert_2d_2_MNase_2.pl col0_MNase_cmt2_CHH_DMR_counts_2.bed 6 0 > col0_MNase_cmt2_CHH_DMR_counts_2.csv

#wc -l ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/col0_methyl_MNase_sorted_dedup_sn_f_mapq.bed 
#987044 /home/shriyas/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/col0_methyl_MNase_sorted_dedup_sn_f_mapq.bed

#wc -l ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/cmt2_methyl_MNase_sorted_dedup_sn_f_mapq.bed 
#749728 /home/shriyas/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/cmt2_methyl_MNase_sorted_dedup_sn_f_mapq.bed


intersectBed -wao -a nucleotideWise_cmt2_CHH_DMR_d25_500bp.bed -b ~/cifs-lab/Shriya/Samples/Jacobsen/Col0-merged/CHH_col0_merged_12.bed > cmt2_CHH_DMR_col0_meth.bed
awk '{OFS="\t"}{if($9=="-1") print $1,$2,$3,$4,$5,$6,"NA"; else {print $1,$2,$3,$4,$5,$6,$10/$11}}' cmt2_CHH_DMR_col0_meth.bed > cmt2_CHH_DMR_col0_meth_temp2.bed
awk '{OFS="\t"}{if($6==1 && $1==1 || $1==2 || $1==3 || $1==4 || $1==5) s+=1; print $0, s}' cmt2_CHH_DMR_col0_meth_temp2.bed > WT_CHH_dataset.bed


intersectBed -wao -a nucleotideWise_cmt2_CHH_DMR_d25_500bp.bed -b ~/cifs-lab/Shriya/Samples/Jacobsen/cmt2/meth_extract/cmt2_CHHmeth.bed > cmt2_CHH_DMR_cmt2_meth.bed
awk '{OFS="\t"}{if($9=="-1") print $1,$2,$3,$4,$5,$6,"NA"; else {print $1,$2,$3,$4,$5,$6,$10/$11}}' cmt2_CHH_DMR_cmt2_meth.bed > cmt2_CHH_DMR_cmt2_meth_temp2.bed
awk '{OFS="\t"}{if($6==1 && $1==1 || $1==2 || $1==3 || $1==4 || $1==5) s+=1; print $0, s}' cmt2_CHH_DMR_cmt2_meth_temp2.bed > cmt2_CHH_dataset.bed

perl ../../convert_2d_2_MNase_2.pl WT_CHH_dataset.bed 6 0 > WT_CHH_dataset.csv
perl ../../convert_2d_2_MNase_2.pl cmt2_CHH_dataset.bed 6 0 > cmt2_CHH_dataset.csv



intersectBed -c -a nucleotideWise_cmt2_CHH_DMR_d25_500bp.bed -b ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/rep2_oct2019/mapped_reads/process_map_reads/col0_rep2_methyl_MNase_sorted_dedup_sn_f_mapq.bed > col0_rep2_MNase_cmt2_d25_CHH_DMR_counts.bed
intersectBed -c -a nucleotideWise_cmt2_CHH_DMR_d25_500bp.bed -b ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/rep2_oct2019/mapped_reads/process_map_reads/cmt2_rep2_methyl_MNase_sorted_dedup_sn_f_mapq.bed > cmt2_rep2_MNase_cmt2_d25_CHH_DMR_counts.bed

awk '{OFS="\t"}{if($6==1) s+=1; print $0,s}' col0_rep2_MNase_cmt2_d25_CHH_DMR_counts.bed > col0_rep2_MNase_cmt2_d25_CHH_DMR_counts_2.bed
awk '{OFS="\t"}{if($6==1) s+=1; print $0,s}' cmt2_rep2_MNase_cmt2_d25_CHH_DMR_counts.bed > cmt2_rep2_MNase_cmt2_d25_CHH_DMR_counts_2.bed

perl ../../convert_2d_2_MNase_2.pl cmt2_rep2_MNase_cmt2_d25_CHH_DMR_counts_2.bed 6 0 > cmt2_rep2_MNase_cmt2_d25_CHH_DMR_counts_2.csv
perl ../../convert_2d_2_MNase_2.pl col0_rep2_MNase_cmt2_d25_CHH_DMR_counts_2.bed 6 0 > col0_rep2_MNase_cmt2_d25_CHH_DMR_counts_2.csv

#wc -l ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/rep2_oct2019/mapped_reads/process_map_reads/col0_rep2_methyl_MNase_sorted_dedup_sn_f_mapq.bed 
#1415365 /home/shriyas/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/rep2_oct2019/mapped_reads/process_map_reads/col0_rep2_methyl_MNase_sorted_dedup_sn_f_mapq.bed

#wc -l ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/rep2_oct2019/mapped_reads/process_map_reads/cmt2_rep2_methyl_MNase_sorted_dedup_sn_f_mapq.bed 
#1067444 /home/shriyas/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/rep2_oct2019/mapped_reads/process_map_reads/cmt2_rep2_methyl_MNase_sorted_dedup_sn_f_mapq.bed




############### CMT3 ###############

### CHG-DMRs #######
cd ~/cifs-lab/Shriya/Samples/manuscript_prep/nucleosome_methylation/working/methyl_mutants_DMRs/All_final_DMRs/cmt3_DMRs

awk '{OFS="\t"} $2>500 {print $1,$2-500,$3+500}' mydiff_CHG_reg_cmt3_col0_d25_q01.bed | coverageBed -d -a stdin -b ~/prog_run/mock_mapped.bed | awk '{OFS="\t"}{print $1,$2+$4-1,$2+$4, $2,$3,$4}' > nucleotideWise_cmt3_CHG_DMR_d25_500bp.bed


intersectBed -c -a nucleotideWise_cmt3_CHG_DMR_d25_500bp.bed -b ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/col0_methyl_MNase_sorted_dedup_sn_f_mapq.bed > col0_MNase_cmt3_CHG_DMR_counts.bed
intersectBed -c -a nucleotideWise_cmt3_CHG_DMR_d25_500bp.bed -b ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/cmt3_methyl_MNase_sorted_dedup_sn_f_mapq.bed > cmt3_MNase_cmt3_CHG_DMR_counts.bed

awk '{OFS="\t"}{if($6==1) s+=1; print $0,s}' col0_MNase_cmt3_CHG_DMR_counts.bed > col0_MNase_cmt3_CHG_DMR_counts_2.bed
awk '{OFS="\t"}{if($6==1) s+=1; print $0,s}' cmt3_MNase_cmt3_CHG_DMR_counts.bed > cmt3_MNase_cmt3_CHG_DMR_counts_2.bed

perl ../../convert_2d_2_MNase_2.pl cmt3_MNase_cmt3_CHG_DMR_counts_2.bed 6 0 > cmt3_MNase_cmt3_CHG_DMR_counts_2.csv
perl ../../convert_2d_2_MNase_2.pl col0_MNase_cmt3_CHG_DMR_counts_2.bed 6 0 > col0_MNase_cmt3_CHG_DMR_counts_2.csv

#wc -l ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/col0_methyl_MNase_sorted_dedup_sn_f_mapq.bed 
#987044 /home/shriyas/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/col0_methyl_MNase_sorted_dedup_sn_f_mapq.bed

#wc -l ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/cmt3_methyl_MNase_sorted_dedup_sn_f_mapq.bed 
#639260 /home/shriyas/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/cmt3_methyl_MNase_sorted_dedup_sn_f_mapq.bed


intersectBed -wao -a nucleotideWise_cmt3_CHG_DMR_d25_500bp.bed -b ~/cifs-lab/Shriya/Samples/Jacobsen/Col0-merged/CHG_col0_merged_12.bed > cmt3_CHG_DMR_col0_meth.bed
awk '{OFS="\t"}{if($9=="-1") print $1,$2,$3,$4,$5,$6,"NA"; else {print $1,$2,$3,$4,$5,$6,$10/$11}}' cmt3_CHG_DMR_col0_meth.bed > cmt3_CHG_DMR_col0_meth_temp2.bed
awk '{OFS="\t"}{if($6==1 && $1==1 || $1==2 || $1==3 || $1==4 || $1==5) s+=1; print $0, s}' cmt3_CHG_DMR_col0_meth_temp2.bed > WT_CHG_dataset.bed


intersectBed -wao -a nucleotideWise_cmt3_CHG_DMR_d25_500bp.bed -b ~/cifs-lab/Shriya/Samples/Jacobsen/cmt3/meth_extract/cmt3_CHGmeth.bed > cmt3_CHG_DMR_cmt3_meth.bed
awk '{OFS="\t"}{if($9=="-1") print $1,$2,$3,$4,$5,$6,"NA"; else {print $1,$2,$3,$4,$5,$6,$10/$11}}' cmt3_CHG_DMR_cmt3_meth.bed > cmt3_CHG_DMR_cmt3_meth_temp2.bed
awk '{OFS="\t"}{if($6==1 && $1==1 || $1==2 || $1==3 || $1==4 || $1==5) s+=1; print $0, s}' cmt3_CHG_DMR_cmt3_meth_temp2.bed > cmt3_CHG_dataset.bed

perl ../../convert_2d_2_MNase_2.pl WT_CHG_dataset.bed 6 0 > WT_CHG_dataset.csv
perl ../../convert_2d_2_MNase_2.pl cmt3_CHG_dataset.bed 6 0 > cmt3_CHG_dataset.csv


#wc -l ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/col0_methyl_MNase_sorted_dedup_sn_f_mapq.bed 
#987044 /home/shriyas/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/col0_methyl_MNase_sorted_dedup_sn_f_mapq.bed

#wc -l ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/cmt3_methyl_MNase_sorted_dedup_sn_f_mapq.bed 
#639260 /home/shriyas/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/mapped_reads/read_processing/cmt3_methyl_MNase_sorted_dedup_sn_f_mapq.bed


intersectBed -c -a nucleotideWise_cmt3_CHG_DMR_d25_500bp.bed -b ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/rep2_oct2019/mapped_reads/process_map_reads/col0_rep2_methyl_MNase_sorted_dedup_sn_f_mapq.bed > col0_rep2_MNase_cmt3_d25_CHG_DMR_counts.bed
intersectBed -c -a nucleotideWise_cmt3_CHG_DMR_d25_500bp.bed -b ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/rep2_oct2019/mapped_reads/process_map_reads/cmt3_rep2_methyl_MNase_sorted_dedup_sn_f_mapq.bed > cmt3_rep2_MNase_cmt3_d25_CHG_DMR_counts.bed

awk '{OFS="\t"}{if($6==1) s+=1; print $0,s}' col0_rep2_MNase_cmt3_d25_CHG_DMR_counts.bed > col0_rep2_MNase_cmt3_d25_CHG_DMR_counts_2.bed
awk '{OFS="\t"}{if($6==1) s+=1; print $0,s}' cmt3_rep2_MNase_cmt3_d25_CHG_DMR_counts.bed > cmt3_rep2_MNase_cmt3_d25_CHG_DMR_counts_2.bed

perl ../../convert_2d_2_MNase_2.pl cmt3_rep2_MNase_cmt3_d25_CHG_DMR_counts_2.bed 6 0 > cmt3_rep2_MNase_cmt3_d25_CHG_DMR_counts_2.csv
perl ../../convert_2d_2_MNase_2.pl col0_rep2_MNase_cmt3_d25_CHG_DMR_counts_2.bed 6 0 > col0_rep2_MNase_cmt3_d25_CHG_DMR_counts_2.csv

#wc -l ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/rep2_oct2019/mapped_reads/process_map_reads/col0_rep2_methyl_MNase_sorted_dedup_sn_f_mapq.bed 
#1415365 /home/shriyas/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/rep2_oct2019/mapped_reads/process_map_reads/col0_rep2_methyl_MNase_sorted_dedup_sn_f_mapq.bed

#wc -l ~/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/rep2_oct2019/mapped_reads/process_map_reads/cmt3_rep2_methyl_MNase_sorted_dedup_sn_f_mapq.bed 
#1295598 /home/shriyas/cifs-lab/Shriya/Samples/methylation_mutants_nucleosomes/rep2_oct2019/mapped_reads/process_map_reads/cmt3_rep2_methyl_MNase_sorted_dedup_sn_f_mapq.bed


paste col0_MNase_cmt3_CHG_DMR_counts_2.bed col0_rep2_MNase_cmt3_d25_CHG_DMR_counts_2.bed | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7+$15,$8}' > col0_rep1+rep2_MNase_cmt3_d25_CHG_DMR_counts.bed
paste cmt3_MNase_cmt3_CHG_DMR_counts_2.bed cmt3_rep2_MNase_cmt3_d25_CHG_DMR_counts_2.bed | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7+$15,$8}' > cmt3_rep1+rep2_MNase_cmt3_d25_CHG_DMR_counts.bed

cd supp_DMRs
paste col0_MNase_cmt3_CG_DMR_counts_2.bed col0_rep2_MNase_cmt3_d25_CG_DMR_counts_2.bed | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7+$15,$8}' > col0_rep1+rep2_MNase_cmt3_d25_CG_DMR_counts.bed
paste cmt3_MNase_cmt3_CG_DMR_counts_2.bed cmt3_rep2_MNase_cmt3_d25_CG_DMR_counts_2.bed | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7+$15,$8}' > cmt3_rep1+rep2_MNase_cmt3_d25_CG_DMR_counts.bed
paste col0_MNase_cmt3_CHH_DMR_counts_2.bed col0_rep2_MNase_cmt3_d25_CHH_DMR_counts_2.bed | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7+$15,$8}' > col0_rep1+rep2_MNase_cmt3_d25_CHH_DMR_counts.bed
paste cmt3_MNase_cmt3_CHH_DMR_counts_2.bed cmt3_rep2_MNase_cmt3_d25_CHH_DMR_counts_2.bed | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7+$15,$8}' > cmt3_rep1+rep2_MNase_cmt3_d25_CHH_DMR_counts.bed
