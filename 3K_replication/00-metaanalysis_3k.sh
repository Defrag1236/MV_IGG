lines=$(cat 20200304_3K_IGP1_23_gwas_list_VIKA_wo_header.csv)

for line in $lines
do
echo $line
igp=$([[ $line  =~ ^IGP[0-9]+ ]] && echo ${BASH_REMATCH[0]})
echo $igp
gwas_ids=$(echo $line | sed "s/IGP[0-9]\{1,2\},//")
echo $gwas_ids

python3 ~/code_folder/gwas_v2/gwas_v2/analysis/meta_analysis/run_meta_analysis.py \
--gwas-ids $gwas_ids \
--output-folder /home/ubuntu/polyomica/projects/MV_igg_replication/3k_meta_23traits/meta_results/ \
--output-filename 3K_meta_$igp.csv \
--maf-filter ON \
--maf-threshold 0.01 \
--genomic-control ON
done
