#!/bin/bash
#SBATCH --job-name="intron_containing_reads_filtering"
#SBATCH --partition=PARTITION_NAME       
#SBATCH --nodes=1
#SBATCH --cpus-per-task=CPUS             
#SBATCH --output=/path/to/output/%x_%j.out
#SBATCH --error=/path/to/output/%x_%j.err
##SBATCH --requeue

echo "SLURM_JOBID=$SLURM_JOBID"

# Load conda (if needed)
# source /path/to/miniforge3/etc/profile.d/conda.sh

# Load modules
module load samtools

# --- USER-DEFINED INPUTS ---
WORK_DIR="/path/to/working_directory"
INPUT_BAM="your_input.bam"
OUT_DIR="output_intron_splicing_analysis"

# List of tRNA reference IDs and their corresponding insertion sequences
declare -A tRNA_insertions

tRNA_insertions["Homo_sapiens_tRNA-Arg-TCT-1-1-WI"]="GAGGCTGAAGGC"
tRNA_insertions["Homo_sapiens_tRNA-Arg-TCT-2-1-WI"]="GTGACGAATAGAGCA"
tRNA_insertions["Homo_sapiens_tRNA-Arg-TCT-3-1-WI"]="GCTGAGCCTAGTGTGGTC"
tRNA_insertions["Homo_sapiens_tRNA-Arg-TCT-5-1-WI"]="GCCTAAATCAAGAG"
tRNA_insertions["Homo_sapiens_tRNA-Ile-TAT-1-1-WI"]="TGACAGTGCGAGCGGAGCA"
tRNA_insertions["Homo_sapiens_tRNA-Ile-TAT-2-1-WI"]="CAGCAGTACATGCAGAGCA"
tRNA_insertions["Homo_sapiens_tRNA-Ile-TAT-3-1-WI"]="GACAGTGCACCTGTGAGCAA"
tRNA_insertions["Homo_sapiens_tRNA-Leu-CAA-1-1-WI"]="CTAAGCTTCCTCCGCGGTGGGGA"
tRNA_insertions["Homo_sapiens_tRNA-Leu-CAA-2-1-WI"]="CTTACTGCTTCCTGTGTTCGGGTC"
tRNA_insertions["Homo_sapiens_tRNA-Leu-CAA-3-1-WI"]="GCTACTTCCCAGGTTTGGGGCTT"
tRNA_insertions["Homo_sapiens_tRNA-Leu-CAA-4-1-WI"]="GTAAGCACCTTGCCTGCGGGCT"
tRNA_insertions["Homo_sapiens_tRNA-Tyr-GTA-1-1-WI"]="TTGGCTGTGTCCTTAGAC"
tRNA_insertions["Homo_sapiens_tRNA-Tyr-GTA-2-1-WI"]="TGGATAGGGCGTGGCA"
tRNA_insertions["Homo_sapiens_tRNA-Tyr-GTA-3-1-WI"]="GCTCATTAAGCAAGGT"
tRNA_insertions["Homo_sapiens_tRNA-Tyr-GTA-4-1-WI"]="TGTATAGACATTTGCGGACAT"
tRNA_insertions["Homo_sapiens_tRNA-Tyr-GTA-5-1-WI"]="CTACTTCCTCAGCAGGAGAC"
tRNA_insertions["Homo_sapiens_tRNA-Tyr-GTA-6-1-WI"]="GGGTTTGAATGTGGTC"
tRNA_insertions["Homo_sapiens_tRNA-Tyr-GTA-7-1-WI"]="CTGCGGAAACGTTTGTGGACA"
tRNA_insertions["Homo_sapiens_tRNA-Tyr-GTA-8-1-WI"]="GTTCATTAAACTAAGGC"
tRNA_insertions["Homo_sapiens_tRNA-Tyr-GTA-9-1-WI"]="GTGCACGCCCGTGGCC"

# Create output directory
mkdir -p "$WORK_DIR/$OUT_DIR"

cd "$WORK_DIR"

# Search for the tRNA introns
for tRNA in "${!tRNA_insertions[@]}"; do
    INSERT_SEQ="${tRNA_insertions[$tRNA]}"
    
    echo "Processing $tRNA with insertion $INSERT_SEQ..."

    # Extract reads mapped to the specific tRNA
    samtools view -b -F 4 "$INPUT_BAM" "$tRNA" > "$OUT_DIR/${tRNA}_reads.bam"
	wc -l "$OUT_DIR/${tRNA}_reads.bam"
    # Convert BAM to SAM (keeping the header)
    samtools view -h "$OUT_DIR/${tRNA}_reads.bam" > "$OUT_DIR/${tRNA}_reads.sam"

    # Filter reads containing an insertion ('I' in CIGAR) and the specific insertion sequence
    if [ "$tRNA" = "Homo_sapiens_tRNA-Leu-CAA-1-1-WI" ]; then
    # Add second sequence for this specific tRNA
		awk -v seq1="$INSERT_SEQ" -v seq2="CTTGGCTTCCTCGTGTTGAGGA" \
			'$0 ~ /^@/ || ($6 ~ /[0-9]+I/ && ($0 ~ seq1 || $0 ~ seq2))' \
			"$OUT_DIR/${tRNA}_reads.sam" > "$OUT_DIR/${tRNA}_filtered.sam"
	else
    # Default case
		awk -v seq="$INSERT_SEQ" \
			'$0 ~ /^@/ || ($6 ~ /[0-9]+I/ && $0 ~ seq)' \
			"$OUT_DIR/${tRNA}_reads.sam" > "$OUT_DIR/${tRNA}_filtered.sam"
	fi


	#awk -v seq="$INSERT_SEQ" '$0 ~ /^@/ || ($6 ~ /[0-9]+I/ && $0 ~ seq)' "$OUT_DIR/${tRNA}_reads.sam" > "$OUT_DIR/${tRNA}_filtered.sam"
	
	wc -l "$OUT_DIR/${tRNA}_filtered.sam"
	
    # Convert filtered SAM to BAM
    samtools view -bS "$OUT_DIR/${tRNA}_filtered.sam" > "$OUT_DIR/${tRNA}_filtered.bam"

    # Sort and index BAM
    samtools sort -o "$OUT_DIR/${tRNA}_filtered.sorted.bam" "$OUT_DIR/${tRNA}_filtered.bam"
    samtools index "$OUT_DIR/${tRNA}_filtered.sorted.bam"
	wc -l "$OUT_DIR/${tRNA}_filtered.sorted.bam"
    echo "Done processing $tRNA!"
done

# Merge all filtered BAMs
samtools merge -f -o "$OUT_DIR/all_filtered_tRNAs.bam" "$OUT_DIR"/*_filtered.sorted.bam
samtools index "$OUT_DIR/all_filtered_tRNAs.bam"

echo "All tRNA processing completed. Results in: $OUT_DIR"
