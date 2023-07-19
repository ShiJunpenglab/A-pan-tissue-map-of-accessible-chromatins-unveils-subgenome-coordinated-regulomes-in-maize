# De novo motif identification
function getmotif(){
    i=$1
    jaspar=./JASPAR2022_CORE_plants_non-redundant_pfms_meme.txt
    meme-chip -meme-minw 5 -meme-maxw 12 -meme-nmotifs 3 -meme-p 20 -dreme-m 0 -spamo-skip -fimo-skip -norc -oc meme_chip/${file/.fa/} ${file} -db ${jaspar}
    genome=./B73v4.fa
}

# motif enrichment
cat filelist.txt |while read name ;do 
	mkdir $name/
	bedtools getfasta -fi B73v4.fa -fo $name".fa" -s -name -bed $id
	ame --verbose 1 --o $name/ --scoring avg --method fisher --hit-lo-fraction 0.25 --evalue-report-threshold 10.0 --control --shuffle-- --kmer 2 $name.fa ./jaspar_plantfdb_DAP_uniq.meme
done
