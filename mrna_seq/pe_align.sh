#!/usr/bin/env bash

use_hisat=true
hiseq_dist_heuristic=100
# novaseq_dist_heuristic is a config option, may be used later
# shellcheck disable=SC2034
novaseq_dist_heuristic=2000
dist_heuristic=$hiseq_dist_heuristic
index_dir=$1
prefix=$2
annotation=$3
output=$4
samples=("${@:5}")
threads=$(($(nproc) / 2))

aligned_files=()
used_samples=()

if [[ "${#samples[@]}" == 1 && -d "${samples[0]}" ]]; then
	echo "Got a folder for samples, checking contents..."
	# word splitting is desired here
	# shellcheck disable=SC2207
	samples=($(ls "${samples[0]}" | grep "_R[12].fq.gz" | sed s/"_R[12].fq.gz"// | sort -u));
	echo "Found samples: ${samples[*]}"
fi;

for sample in "${samples[@]}"; do
	aligned_file="${prefix}${sample}.bam"
	if [[ -f "${prefix}${sample}_R1.fq.gz" && -f "${prefix}${sample}_R2.fq.gz" ]]; then
		echo "Raw FASTQ.GZ files exist, adapter trimming, base correcting and aligning with HISAT2";
		fastp --in1 "${prefix}${sample}_R1.fq.gz" --in2 "${prefix}${sample}_R2.fq.gz" \
			--out1 "${prefix}${sample}_R1.trim.fq.gz" --out2 "${prefix}${sample}_R2.trim.fq.gz" \
			-q 20 -u 10 --detect_adapter_for_pe --dont_eval_duplication --correction --thread $threads;

		aligner='bowtie2'
		if $use_hisat; then
			aligner='hisat2'
		fi
		
		$aligner -x "${index_dir}" -1 "${prefix}${sample}_R1.trim.fq.gz" -2 "${prefix}${sample}_R2.trim.fq.gz" -p $threads |
		samtools collate -f -@ $threads -O - |
		samtools fixmate -rcm -@ $threads - - |
		samtools sort -l 9 -@ $threads |
		samtools markdup -l 150 -d "${dist_heuristic}" -rS -@ $threads - - > "$aligned_file";
		samtools index -@ $threads "$aligned_file"
		used_samples+=("$sample")
		aligned_files+=("$aligned_file")
	else
		echo "Unsupported file found, skipping..."
	fi
done;

if [[ "${#used_samples[@]}" > 1 ]]; then
	echo "Running featureCounts...";

	# needs featureCounts from Subread version 2.0.2+ for --countReadPairs
	featureCounts -T $threads -a "$annotation" -p --countReadPairs -C -M -O --fracOverlap 0.01 --fraction --ignoreDup --largestOverlap -o "$output" "${aligned_files[@]}";

	echo "Done, featureCounts quantification saved to $output";
fi