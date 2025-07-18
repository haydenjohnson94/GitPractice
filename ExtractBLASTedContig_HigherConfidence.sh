#!/bin/bash
#SBATCH --job-name=extract_BLAST
#SBATCH --output=extract_BLAST_%j.out
#SBATCH --error=extract_BLAST_%j.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16

# Define input and output paths
BLAST_RESULTS="/nfs4/BPP/Uehling_Lab/johnsh/MortierellaceaeMating_Kyle/BLAST_Files/blast_results.csv"
GENOME_DIR="/nfs4/BPP/Uehling_Lab/johnsh/MortierellaceaeMating_Kyle/Genomes"
OUTPUT_FILE="target_contigs.fna"

# Clear output file
> "$OUTPUT_FILE"

# Process BLAST results: filter for pident > 90, aln_length > 500
awk -F, 'NR>1 && $3 > 90 && $4 > 500' "$BLAST_RESULTS" | while IFS=, read -r qseqid sseqid pident aln_length mismatch gapopen qstart qend sstart send evalue bitscore stitle
do
    echo "Processing: $sseqid"
    
    # Parse accession and contig from sseqid
    accession="${sseqid:0:14}"
    contig="${sseqid:14}"
    contig="${contig#_}"
    
    # Clean accession for file matching
    accession_clean=$(echo "$accession" | awk -F'_' '{print $1 "_" $2}')
    
    echo "Cleaned Accession: $accession_clean"
    echo "Contig: $contig"
    
    # Find genome file matching accession
    genome_file=$(find "$GENOME_DIR" \( -type f -o -type l \) \( -iname "${accession_clean}*" -a \( -iname "*.fna" -o -iname "*.fasta" \) \) | head -1)
    
    if [[ -z "$genome_file" ]]; then
        echo "ERROR: Genome file not found for $accession_clean" >&2
        echo "DEBUG: Tried to find files matching: ${accession_clean}*" >&2
        continue
    fi
    
    echo "Found genome: $genome_file"
    
    # Write FASTA header for full contig
    echo ">${sseqid}|full_contig" >> "$OUTPUT_FILE"
    
    # Extract entire contig sequence
    awk -v contig="$contig" '
        BEGIN {
            RS=">";
            found=0;
            gsub(/[_.]/, "", contig);
        }
        NR > 1 {
            split($0, lines, "\n");
            header = lines[1];
            seq = "";
            for (i = 2; i <= length(lines); i++) {
                seq = seq lines[i];
            }
            contig_header = header;
            gsub(/[ _.:]/, "", contig_header);
            if (index(contig_header, contig)) {
                print seq;
                found = 1;
                exit;
            }
        }
        END {
            if (!found) exit 1;
        }
    ' "$genome_file" >> "$OUTPUT_FILE" 2>/dev/null || echo "ERROR: Failed to extract $contig" >&2
done

# Final processing: index output with samtools faidx
if [[ -s "$OUTPUT_FILE" ]]; then
    samtools faidx "$OUTPUT_FILE"
    echo "Done. Output in $OUTPUT_FILE"
else
    echo "ERROR: No sequences extracted" >&2
    exit 1
fi