
#!/bin/bash

job_number=$1
Ne=$2
p=$3

# Set the input tree and table.
tree="Tanganyika_b1.tre"

run_number_start=$(( $job_number * 10 ))
run_number_end=$(( $run_number_start + 10 ))

# do ten runs:
for (( run_number=$run_number_start; run_number<$run_number_end; run_number++ )); do

# Perform simulations with different migration matrix subsamples
for mm in {1..20}; do
	migration_matrix=Fmatrix_full_sample_${mm}.txt
	# Set the vcf output file.
        vcf=simulate.${run_number}_${mm}_twoSampleWithRecomb_Ne_${Ne}_p_${p}.vcf

	python3 simulate.py --num_indiv 2 --Ne ${Ne} --chr_length 10000 --gene_flow_period ${p} ${tree} migration_matrices/${migration_matrix} ${vcf}
done

done
