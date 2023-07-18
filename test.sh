source activate /data21/liuyt/miniconda/envs/featurecounts
# ###get SAF file###
# awk -vOFS='\t' '{if($6~/-/)print $4,$1,$3,$3+2000,$6;else print $4,$1,$2-2000,$2,$6}' Zea_mays.chr.B73_RefGen_v4.49_rmnoncoding_gene.closest.gtf >B73V4.gene_up2kb.saf

##get read count matrix
featureCounts -T 30 -F SAF -a B73V4.gene_up2kb.saf -o gene_up2kb.read.count.PE -p -B -C  *sort.rmdup.bam

# featureCounts -T 30 -F SAF -a B73V4.gene_up2kb.saf -o gene_up2kb.read.count.PE -p -B -C -M -O *sort.rmdup.bam
featureCounts -T 30 -F SAF -a B73V4.gene_up2kb.saf -o B73leafmesophyll.gene_up2kb.read.count  B73leafmesophyll.pe.q10.sort.rmdup.bam ##手动删除第一行 
cut -f 7 B73leafmesophyll.gene_up2kb.read.count >B73leafmesophyll.gene_up2kb.read.count.temp
paste gene_up2kb.read.count.PE B73leafmesophyll.gene_up2kb.read.count.temp >gene_up2kb.read.count #删除单端数据 添加新的

Rscript FPKM.R ####可自行修改表头 去除NA 得到FPKM TPM CPM file
grep -v B7312DAPP coutsgene_fpkm.csv |sort -k 1,1 |awk -vOFS="\t" 'BEGIN{print "gene","B7312DAPP","B7312DAPPR","B7312DAPZP","B73ear1cm","B73ear4mm","B73leaf7days","B73root6days","B73SAM","B73seedingFresh","B73seedingFrozen","B73seedingleaf2weeks","B73silk","B73stem","B73tassel4mm","B73leafmesophyll"}{print $0}' |sed 's/\"//g' >coutsgene_fpkm_sort.csv


##get RNA FPKM matrix ### 
##RNA matrix

ls *.rmnoncoding.fpkm |while read id;
do
# ACRname=$(echo $id |cut -d ' ' -f 2 )
RNAname=$(echo $id |cut -d '.' -f 1 )
cut -f 1 B73V4.TF.list.txt | \
grep -w -f - $RNAname.rmnoncoding.fpkm | \
sort -k 6,6 | \
awk -vOFS='\t' -v RNAname=$RNAname 'BEGIN{print "gene",RNAname}{print $NF,$4}' >$RNAname.RNA.Rinput
cut -f 2 $RNAname.RNA.Rinput >$RNAname.RNA_fpkm.Rinput.temp
done
# paste coutsgene_fpkm_sort.csv *.RNA_fpkm.Rinput.temp >ACR_RNA_couts_fpkm_sort.csv
paste geneid.txt *.RNA_fpkm.Rinput.temp >RNA_couts_fpkm_sort.csv
rm *temp
#####3对TF例子的热图

head -1 RNA_couts_fpkm_sort.csv >TFexample.heatmap.Rinput
cat TFexample.list.txt |while read id ;do grep -w $id RNA_couts_fpkm_sort.csv;done >>TFexample.heatmap.Rinput
# grep -w -f TFexample.list.txt RNA_couts_fpkm_sort.csv >>TFexample.heatmap.Rinput
####GLK基因的表达
head -1 RNA_couts_fpkm_sort.csv >RNA_couts_fpkm_sort.header
cut -f 1 G2-like.list.txt |grep -w -f - M1.subgenome.geneid |grep -w -f - RNA_couts_fpkm_sort.csv >GLK_M1.fpkm.temp
cat RNA_couts_fpkm_sort.header GLK_M1.fpkm.temp >GLK_M1.fpkm
cut -f 1 G2-like.list.txt |grep -w -f - M2.subgenome.geneid |grep -w -f - RNA_couts_fpkm_sort.csv >GLK_M2.fpkm.temp
cat RNA_couts_fpkm_sort.header GLK_M2.fpkm.temp >GLK_M2.fpkm
rm *temp


#####全部基因的表达矩阵和计算香农熵
cut -f 6 seedling_14days.rmnoncoding.fpkm |sort -k 6,6  >totalgeneid.txt
sed -i '1i gene' totalgeneid.txt
ls *.rmnoncoding.fpkm |while read id ;
do
RNAname=$(echo $id |cut -d '.' -f 1 )
sort -k 6,6 $RNAname.rmnoncoding.fpkm | \
awk -vOFS='\t' -v RNAname=$RNAname 'BEGIN{print "gene",RNAname}{print $NF,$4}' >$RNAname.totalgene_RNA.Rinput
cut -f 2 $RNAname.totalgene_RNA.Rinput >$RNAname.totalgene_RNA_fpkm.Rinput.temp
done
paste totalgeneid.txt *totalgene_RNA_fpkm.Rinput.temp >totalgene_RNA_couts_fpkm_sort.csv
rm *temp
Rscript Entropy.R #计算每个基因的香农熵
python TF_gene.anno.py -i CS_Entropy.txt -o control_TF.Entropy.txt -id B73V4.TF.list.txt
sed -i '1i gene\tentropy\ttype' control_TF.Entropy.txt


#ACR matrix
head -1 coutsgene_fpkm_sort.csv >ACR_couts_fpkm_sort.csv
cut -f 1 B73V4.TF.list.txt | \
grep -w -f - coutsgene_fpkm_sort.csv >>ACR_couts_fpkm_sort.csv

###compare expression of TF gene and others ###
ls *rmnoncoding.fpkm |while read id;
do
RNAname=$(echo $id |cut -d '.' -f 1 )
awk -vOFS='\t' -v RNAname=$RNAname '{print $NF,$4,RNAname}' $id >$id.temp
done
ls *rmnoncoding.fpkm.temp |while read id ;do
RNAname=$(echo $id |cut -d '.' -f 1 )
python TF_gene.anno.py -i $id -id B73V4.TF.list.txt -o $RNAname.anno.RNA_fpkm
done
cat *.anno.RNA_fpkm >tissue.anno.RNA_fpkm.Rinput
ls *.anno.RNA_fpkm |while read id ;do name=$(echo $id |cut -d "." -f 1) ;grep -w 'Other' $id |cut -f 1 |grep -w -f - Zea_mays.chr.B73_RefGen_v4.49_rmnoncoding_gene.closest.gtf >$name.other_id.bed;done
rm *temp
sed -i '1i gene\tFPKM\ttissue\ttype' tissue.anno.RNA_fpkm.Rinput

###compare accessibility of TF gene and others ###
# 提取每个列到一个新文件中，以列名命名文件
awk -vOFS="\t" 'NR==1 {for(i=3;i<=NF;i++) {name[i]=$i}} NR>1{for(i=3;i<=NF;i++) {print $1,$i,name[i] > name[i]".ACR_fpkm.temp"}}' coutsgene_fpkm_sort.csv
ls *.ACR_fpkm.temp |while read id ;do
ACRname=$(echo $id |cut -d '.' -f 1 )
python TF_gene.anno.py -i $id -id B73V4.TF.list.txt -o $ACRname.anno.ACR_fpkm
done
cat *.anno.ACR_fpkm >tissue.anno.ACR_fpkm.Rinput
rm *temp
sed -i '1i gene\tFPKM\ttissue\ttype' tissue.anno.ACR_fpkm.Rinput

# ###ACR gene body enrichment 
source activate lyt_atac
cat ACR_RNA.name.txt|while read id ;
do
ACRname=$(echo $id |cut -d ' ' -f 2 )
RNAname=$(echo $id |cut -d ' ' -f 1 )
computeMatrix scale-regions -p 40 -b 5000 -a 5000 --regionBodyLength 5000 -R $RNAname.other_id.bed B73V4.TF_gene.bed -S $ACRname.pe.q10.sort.rmdup.shift.norm.bw --skipZeros -o $ACRname.randomgene_TF.body.gz
plotProfile -m  $ACRname.randomgene_TF.body.gz -out $ACRname.randomgene_TF.body_Profile.png
done
####################test
computeMatrix scale-regions -p 40 -b 5000 -a 5000 --regionBodyLength 5000 -R B73M1.homo_gene.bed B73M2.homo_gene.bed B73M1.TF_gene_homo.bed B73M2.TF_gene_homo.bed B73M1.gene.bed B73M2.gene.bed -S B73ear4mm.pe.q10.sort.rmdup.shift.norm.bw --skipZeros -o B73ear4mm.all.body.gz
plotProfile -m  B73ear4mm.all.body.gz -out B73ear4mm.all.body_Profile.png


###subgenome ACR enrichment
bedtools intersect -a Zea_mays.chr.B73_RefGen_v4.49_rmnoncoding_gene.closest.gtf -b B73_subgenome1.sorted.merged_5Mb.bed -f 0.5 >B73M1.gene.bed
bedtools intersect -a Zea_mays.chr.B73_RefGen_v4.49_rmnoncoding_gene.closest.gtf -b B73_subgenome2.sorted.merged_5Mb.bed -f 0.5 >B73M2.gene.bed
grep -w -f M1.updata.geneid Zea_mays.chr.B73_RefGen_v4.49_rmnoncoding_gene.closest.gtf >B73M1.homo_gene.bed
grep -w -f M2.updata.geneid Zea_mays.chr.B73_RefGen_v4.49_rmnoncoding_gene.closest.gtf >B73M2.homo_gene.bed

cat ACR_RNA.name.txt|while read id ;
do
ACRname=$(echo $id |cut -d ' ' -f 2 )
computeMatrix scale-regions -p 40 -b 5000 -a 5000 --regionBodyLength 5000 -R B73M1.gene.bed B73M2.gene.bed -S $ACRname.pe.q10.sort.rmdup.shift.norm.bw --skipZeros -o $ACRname.subgenome.body.gz
plotProfile -m  $ACRname.subgenome.body.gz -out $ACRname.subgenome.body_Profile.png
###
computeMatrix scale-regions -p 40 -b 5000 -a 5000 --regionBodyLength 5000 -R B73M1.homo_gene.bed B73M2.homo_gene.bed -S $ACRname.pe.q10.sort.rmdup.shift.norm.bw --skipZeros -o $ACRname.conserved.body.gz
plotProfile -m  $ACRname.conserved.body.gz -out $ACRname.conserved.body_Profile.png
done

######subgenome TF
python extract_homoTF.py M1_M2_sorghum.geneid B73V4.TF.list.txt >homoTF.geneid
cut -f 1 B73V4.TF.list.txt |grep -w -f - B73M1.gene.bed >B73M1.TF_gene.bed
cut -f 1 B73V4.TF.list.txt |grep -w -f - B73M2.gene.bed >B73M2.TF_gene.bed

computeMatrix scale-regions -p 40 -b 5000 -a 5000 --regionBodyLength 5000 -R B73M1.TF_gene.bed B73M2.TF_gene.bed -S *pe.q10.sort.rmdup.shift.norm.bw --skipZeros -o subgenome_TF.body.gz
plotProfile -m  subgenome_TF.body.gz -out subgenome_TF.body_Profile.png

cut -f 1 homoTF.geneid |grep -w -f - Zea_mays.chr.B73_RefGen_v4.49_rmnoncoding_gene.closest.gtf >B73M1.TF_gene_homo.bed
cut -f 2 homoTF.geneid |grep -w -f - Zea_mays.chr.B73_RefGen_v4.49_rmnoncoding_gene.closest.gtf >B73M2.TF_gene_homo.bed

ls *pe.q10.sort.rmdup.shift.norm.bw |while read id;
do
name=$(echo $id |cut -d '.' -f 1 )
computeMatrix scale-regions -p 40 -b 5000 -a 5000 --regionBodyLength 5000 -R B73M1.TF_gene_homo.bed B73M2.TF_gene_homo.bed -S $id --skipZeros -o $name.conserved_TF.body.gz
plotProfile -m  $name.conserved_TF.body.gz -out $name.conserved_TF.body_Profile.png
computeMatrix scale-regions -p 40 -b 5000 -a 5000 --regionBodyLength 5000 -R B73M1.TF_gene.bed B73M2.TF_gene.bed -S $id --skipZeros -o $name.subgenome_TF.body.gz
plotProfile -m  $name.subgenome_TF.body.gz -out $name.subgenome_TF.body_Profile.png
done



########################################################
###bam文件计数 reads
bedtools makewindows -g B73V4.chromosome.sizes -w 100000 | \
awk -vFS="\t" -vOFS="\t" '{print $1,$2,$3}' | \
bedtools sort -i - > B73V4.chromosome.sizes.100kb
samtools index B73ear4mm.pe.q10.sort.rmdup.bam
bedtools multicov -bams B73ear4mm.pe.q10.sort.rmdup.bam -bed B73V4.chromosome.sizes.100kb  > circos.test.txt
awk -F '\t' -vOFS='\t' '{print $1,$2,$3,$4/100000}' circos.test.txt >circos.test.density.txt
