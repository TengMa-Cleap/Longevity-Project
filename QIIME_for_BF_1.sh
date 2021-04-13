#!/bin/bash
# last update 2015.3.23 
# Hou qiangchuan
# unzip core_set_aligned.fasta.zip
# rm core_set_aligned.fasta.zip
python .150923Calculate_the_sequence_numbers_of_each_sample.py -i final.fasta -o reads_numbers_of_each_samples.txt
rm -rf .150923Calculate_the_sequence_numbers_of_each_sample.py
echo "reads of final.fasta: "
grep ">" final.fasta | wc -l
# cp core_set_aligned.fasta.imputed /userdata/tmp/
# cp core_set_aligned.fasta.imputed /tmp/
# parallel_align_seqs_pynast.py -i final.fasta -o pynast_align_1/ -t core_set_aligned.fasta.imputed -d /tmp/core_set_aligned.fasta.imputed -a uclust -e 100 -O 40
# cd pynast_align_1/
# echo "reads of final_aligned.fasta: "
# grep ">" final_aligned.fasta | wc -l
# echo "reads of final_failures.fasta: "
# grep ">" final_failures.fasta | wc -l
# sed 's/-//g' final_aligned.fasta > Seqs_neimeng_all_final_aligned_no_gap_final.fasta
# echo "reads of Seqs_neimeng_all_final_aligned_no_gap_final.fasta: "
# grep ">" Seqs_neimeng_all_final_aligned_no_gap_final.fasta | wc -l
# cp Seqs_neimeng_all_final_aligned_no_gap_final.fasta ../
# cd ..
pick_otus.py -i final.fasta -m uclust -s 1.00 --enable_rev_strand_match -o 0.3uclust_otu_100/
cd 0.3uclust_otu_100
mv final_clusters.uc uclust_otu_100_cluster.uc
mv final_otus.txt uclust_otu_100.txt
mv final_otus.log uclust_otu_100_log.txt
pick_rep_set.py -i uclust_otu_100.txt -f ../final.fasta -o uclust_otu_100_rep_set_first.fna
wc -l uclust_otu_100.txt
wc -l uclust_otu_100_rep_set_first.fna
echo "reads of uclust_otu_100_rep_set_first.fna: "
grep ">" uclust_otu_100_rep_set_first.fna |wc -l
cd ..
pick_otus.py -i 0.3uclust_otu_100/uclust_otu_100_rep_set_first.fna --enable_rev_strand_match -s 1.00 -m uclust -o 0.4uclust_otu_97/  #0.99 for pacbio seqence data
merge_otu_maps.py -i 0.3uclust_otu_100/uclust_otu_100.txt,0.4uclust_otu_97/uclust_otu_100_rep_set_first_otus.txt -o 0.4uclust_otu_97/uclust_otus_100_97_merged.txt
pick_rep_set.py -i 0.4uclust_otu_97/uclust_otus_100_97_merged.txt -f final.fasta -o uclust_otus_100_97_merged_rep_set.fasta
mkdir 0.5rep_set
mv uclust_otus_100_97_merged_rep_set.fasta 0.5rep_set/
cp core_set.fasta.imputed 0.5rep_set/
cd 0.5rep_set/
#identify_chimeric_seqs.py -i uclust_otus_100_97_merged_rep_set_aligned.fasta -a core_set_aligned.fasta.imputed -o chimeric_seq.txt
identify_chimeric_seqs.py -m usearch61 -i uclust_otus_100_97_merged_rep_set.fasta -o usearch61_chimera_checking --suppress_usearch61_ref
#touch chimeric_seq.txt
cp usearch61_chimera_checking/chimeras.txt ./chimeric_seq.txt
cd ..
filter_fasta.py -f 0.5rep_set/uclust_otus_100_97_merged_rep_set.fasta -o uclust_otus_100_97_merged_rep_set_aligned_no_chimera.fasta -s 0.5rep_set/chimeric_seq.txt -n
mv uclust_otus_100_97_merged_rep_set_aligned_no_chimera.fasta 0.5rep_set/
sed 's/-//g' 0.5rep_set/uclust_otus_100_97_merged_rep_set_aligned_no_chimera.fasta > 0.5rep_set/uclust_otus_100_97_merged_rep_set_aligned_no_chimera_no_gap.fasta
wc -l 0.5rep_set/chimeric_seq.txt
echo "reads of uclust_otus_100_97_merged_rep_set_aligned_no_chimera.fasta: "
grep ">" 0.5rep_set/uclust_otus_100_97_merged_rep_set_aligned_no_chimera.fasta |wc -l
echo "reads of uclust_otus_100_97_merged_rep_set_aligned_no_chimera_no_gap.fasta: "
grep ">" 0.5rep_set/uclust_otus_100_97_merged_rep_set_aligned_no_chimera_no_gap.fasta |wc -l
assign_taxonomy.py -i 0.5rep_set/uclust_otus_100_97_merged_rep_set_aligned_no_chimera_no_gap.fasta -r Bifidobacterium_RPSK.fasta -t Bifidobacterium_taxonomy.txt -e 0.00000000001 -m blast -o blast_assigned_taxonomy
mkdir 0.6blast_assigned_taxonomy
awk '{$NF="";print $0}' blast_assigned_taxonomy/uclust_otus_100_97_merged_rep_set_aligned_no_chimera_no_gap_tax_assignments.txt | awk '{$NF="\t1";print $0}' > 0.6blast_assigned_taxonomy/uclust_otus_tax_assignments.txt
sed  's/ /\t/' 0.6blast_assigned_taxonomy/uclust_otus_tax_assignments.txt | sed 's/No blast hit/unclassified/' > 0.6blast_assigned_taxonomy/uclust_otus_tax_assignments_with_unclassified.txt
make_otu_table.py -i 0.4uclust_otu_97/uclust_otus_100_97_merged.txt -o 0.6blast_assigned_taxonomy/otu_table.txt -e 0.5rep_set/chimeric_seq.txt -t 0.6blast_assigned_taxonomy/uclust_otus_tax_assignments_with_unclassified.txt
# cd 0.6rdp_assigned_taxonomy_greengene/
cd 0.6blast_assigned_taxonomy
convert_biom.py -i otu_table.txt -o otu_table1.txt -b --header_key taxonomy
sed 's/[kpcofgs]__//g' otu_table1.txt > otu_table_greengene1.txt
max_reads=`awk -F "\t" 'NR >=3 {for(i=2;i<NF;i++)a[i]+=$i}END{for(i=2;i<NF;i++)print a[i]}' otu_table1.txt | sort -n| tail -1`
min_reads=`awk -F "\t" 'NR >=3 {for(i=2;i<NF;i++)a[i]+=$i}END{for(i=2;i<NF;i++)print a[i]}' otu_table1.txt | sort -n| head -1`
awk -F '\t' '{print $NF}' otu_table1.txt > taxonomy.txt
mv ../.151130Calculate_abundance_of_each_taxonomy_from_otutable.py ./
python .151130Calculate_abundance_of_each_taxonomy_from_otutable.py -i otu_table1.txt -o taxonomy_abundance
rm -rf .151130Calculate_abundance_of_each_taxonomy_from_otutable.py
cd ..
mkdir 0.7fasttree_phylogeny
mafft --auto 0.5rep_set/uclust_otus_100_97_merged_rep_set_aligned_no_chimera.fasta > 0.7fasttree_phylogeny/uclust_otus_100_97_merged_rep_set_aligned_no_chimera_aligned.fasta
filter_alignment.py -i 0.7fasttree_phylogeny/uclust_otus_100_97_merged_rep_set_aligned_no_chimera_aligned.fasta -o 0.7fasttree_phylogeny/ -g 0.9999999 --suppress_lane_mask_filter
mv  0.7fasttree_phylogeny/uclust_otus_100_97_merged_rep_set_aligned_no_chimera_aligned_pfiltered.fasta 0.7fasttree_phylogeny/uclust_otus_100_97_merged_rep_set_aligned_no_chimera_pfiltered.fasta 
#wc -l 0.7fasttree_phylogeny/uclust_otus_100_97_merged_rep_set_aligned_no_chimera_pfiltered.fasta 
echo "reads of uclust_otus_100_97_merged_rep_set_aligned_no_chimera_pfiltered.fasta: "
grep ">" 0.7fasttree_phylogeny/uclust_otus_100_97_merged_rep_set_aligned_no_chimera_pfiltered.fasta -c

make_phylogeny.py -i 0.7fasttree_phylogeny/uclust_otus_100_97_merged_rep_set_aligned_no_chimera_pfiltered.fasta -o 0.7fasttree_phylogeny/fasttree.tre -l 0.7fasttree_phylogeny/fasttree.log -t fasttree

python .161019generate_n_rarefaction_process.py -i 0.6blast_assigned_taxonomy/otu_table.txt -o 0.8multiple_rarefactions -d $min_reads -t 20 
rm -rf .161019generate_n_rarefaction_process.py
python .161020biom2txt.py -i 0.8multiple_rarefactions
rm -rf biom_to_txt.R .161020biom2txt.py
python .161025merge_all_txt.py -i 0.8multiple_rarefactions
rm .161025merge_all_txt.py
cd ./0.8multiple_rarefactions
head -n 2 ./0.8multiple_rarefactions_0/rarefaction_*_0.biom.txt > head.txt
cat head.txt merge.txt > head_merge.txt
paste -d "\t" head_merge.txt ../0.6blast_assigned_taxonomy/taxonomy.txt > otu_table.taxonomy.txt
convert_biom.py -i otu_table.taxonomy.txt -o otu_table.from_txt.biom --biom_table_type="otu table" --process_obs_metadata taxonomy
mv otu_table.from_txt.biom ../
cd ..

beta_diversity_through_plots.py -i otu_table.from_txt.biom -m beta_diversity_mapping_file_20111111.txt -o beta_diversity/ -p beta_diversity_qiime_parameters.txt -t 0.7fasttree_phylogeny/fasttree.tre -a -O 10 -f

#nohup beta_diversity_through_plots.py -i otu_table.from_txt.biom -m mapping_file.txt -o beta_diversity/ -p beta_diversity_qiime_parameters.txt -t 0.7fasttree_phylogeny/fasttree.tre -a -f &
cd beta_diversity/
mkdir pcoa/
mv *pc.txt pcoa/
mkdir emperor_3D_plot
mv *emperor_pcoa_plot ./emperor_3D_plot
mkdir distance_matrix
mv *_dm.txt distance_matrix/
upgma_cluster.py -i distance_matrix/ -o upgma_cluster/
cd ..
parallel_multiple_rarefactions.py -i 0.6blast_assigned_taxonomy/otu_table.txt -o 01.multiple_rare/ -m 10 -x $max_reads -s 100 -n 50 -O 60
parallel_alpha_diversity.py -i 01.multiple_rare/ -o 02.alpha_div/ -m chao1,chao1_confidence,observed_species,shannon,simpson,PD_whole_tree -O 60 -t 0.7fasttree_phylogeny/fasttree.tre
collate_alpha.py -i 02.alpha_div/ -o results/ -e 02.alpha_div/alpha_rarefaction_10_0.txt
cd results/
mv ../.150906Calculate_α_diversity_average.py ./
python .150906Calculate_α_diversity_average.py -i shannon.txt -o shannon_average.txt
python .150906Calculate_α_diversity_average.py -i simpson.txt -o simpson_average.txt
python .150906Calculate_α_diversity_average.py -i observed_species.txt -o observed_species_average.txt
python .150906Calculate_α_diversity_average.py -i chao1.txt -o chao1_average.txt
rm -rf .150906Calculate_α_diversity_average.py