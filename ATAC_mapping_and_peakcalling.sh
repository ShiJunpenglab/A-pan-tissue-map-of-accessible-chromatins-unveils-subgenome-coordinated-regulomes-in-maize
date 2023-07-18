#!/bin/bash

# Mapping & QC
python atac-pipe.py --MappingQC -t 25 -i ./samples_direc -o ./mapping_output -r B73v4

# Merge replicate
python atac-pipeB73.py --Merge --bam ./mapping_output -o ./merge_replicate_direc -r B73v4 --group samples_group.txt

# peak calling
python atac-pipeB73.py --PeakCalling --bed ./merge_replicate_direc -o ./peak_direc -r B73v4 --group samples_group.txt

# Peak filter (using Tn5 integration site densities)
cd ./peak_direc
ls *pe.q10.sort.rmdup.shift.bed |while read bed
do
sample_name=$(echo $bed |cut -d "." -f 1)
input_narrowpeak=${sample_name}_peaks.narrowPeak
genome_fai=$1
fdr=$2

# extract Tn5 site of read 
awk 'BEGIN {OFS = "\t"} ; {print $1,($2+$3)/2,($2+$3)/2+1}' $bed > $bed.Tn5

# calculate tn5 site density of ACR（tn5 count/ ocr length)
bedtools coverage -a $input_narrowpeak -b $bed.Tn5 | awk '{print $1,$2,$3,$4,$11}' > input_${sample_name}_peaks_and_tn5count.bed

##bedtools shuffle specifically excluding ACRs from the randomized selection space
bedtools shuffle -i $input_narrowpeak -g $genome_fai |awk '{print $1,$2,$3}' | tr ' ' '\t' > ${sample_name}.shuffle_randomized.bed

# calculate tn5 density of randomized region
bedtools coverage -a ${sample_name}.shuffle_randomized.bed -b $bed.Tn5 | awk '{print $1,$2,$3,$4}' > input_${sample_name}_randomized_and_tn5count.bed #four column chr start end tn5count
Rscript ./filter_eFDR.R input_${sample_name}_peaks_and_tn5count.bed input_${sample_name}_randomized_and_tn5count.bed $fdr ${sample_name}_peaks_filted_by_tn5density.fdr$fdr.bed
##usage example: bash filter_ocr_by_tn5_density.sh B73v4.fa.fai 0.01

done

###add more information (like pvalue)
ls *_filted_by_tn5density.fdr0.01.bed |while read id
do
name=$(echo $id |cut -d "_" -f 1)
awk '{print $4}' $id >$id.fdr0.01.idex
grep -w -f  $id.fdr0.01.idex ${name}_peaks.narrowPeak >filter_OCR/${name}_peaks.narrowPeak
done
