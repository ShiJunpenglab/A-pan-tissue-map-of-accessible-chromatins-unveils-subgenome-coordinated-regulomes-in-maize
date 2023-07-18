# Mapping & QC
python atac-pipe.py --MappingQC -t 25 -i ./samples_direc -o ./mapping_output -r B73v4

# merge replicate
python atac-pipeB73_mt.py --Merge --bam ./mapping_output -o ./merge_replicate_direc -r B73v4

# peak calling
