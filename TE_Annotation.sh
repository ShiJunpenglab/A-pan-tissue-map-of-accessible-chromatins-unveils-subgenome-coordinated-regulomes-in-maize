#!/bin/bash

# EDTA
perl /EDTA/EDTA.pl --genome B73v4.chr.fa --species Maize --anno 1 -t 20

# Further classifies unknown transposons
#step 1
grep "LTR/unknown" B73v4.chr.fa.mod.EDTA.TElib.fa | sed 's/>//' | seqtk subseq B73v4.chr.fa.mod.EDTA.TElib.fa - > B73v4.LTR.unknown.fa

#step2
/public1/home/sc80041/miniconda3/envs/py36-gan/bin/python3.6 /public1/home/sc80041/ganyuanxian/bio/DeepTE-master/DeepTE.py -d ./deepte -i B73v4.LTR.unknown.fa -o ./deepteout -sp P -m_dir /public1/home/sc80041/ganyuanxian/bio/DeepTE-master/Model_dir_self/Plants_model -fam LTR

#step3 extract kown TEs
grep -v "LTR/unknown" B73v4.chr.fa.mod.EDTA.TElib.fa | sed 's/>//' | seqtk subseq B73v4.chr.fa.mod.EDTA.TElib.fa - > B73v4.LTR.known.fa
 
#step4
sed 's/LTR\/unknown__ClassI_LTR_Copia/LTR\/Copia/' ./deepteout/opt_DeepTE.fasta | sed 's/LTR\/unknown__ClassI_LTR_Gypsy/LTR\/Gypsy/'|sed 's/LTR\/unknown__ClassI_LTR/LTR\/unknown/' > ./deepteout/B73v4.opt_DeepTE.LTR.unknown.fa
sed 's/LTR\/unknown__ClassI_LTR_Copia/LTR\/Copia/' ./deepteout/opt_DeepTE.fasta | sed 's/LTR\/unknown__ClassI_LTR_Gypsy/LTR\/Gypsy/' | sed 's/LTR\/unknown__ClassI_LTR/LTR\/unknown/' > ./B73v4.opt_DeepTE.unkown.rename.fa &&
cat ./B73v4.LTR.known.fa ./B73v4.opt_DeepTE.unkown.rename.fa > ./B73v4.fa.mod.EDTA.TElib.new.fa
