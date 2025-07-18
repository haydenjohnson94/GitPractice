#!/bin/bash
#SBATCH --job-name=illumina_assembly_Polish
#SBATCH --output=illumina_assembly_Polish_%j.out
#SBATCH --error=illumina_assembly_Polish_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --partition=uehling_lab
#SBATCH --mem=200G

THREADS=32


# Define directories
READS_DIR="/nfs4/BPP/Uehling_Lab/NewGenomes/CombineLanes_AllReads"
ASSEMBLY_DIR="/nfs4/BPP/Uehling_Lab/NewGenomes/CombineLanes_AllReads/Assembly_Polish"
SAM_BAM_DIR="/nfs4/BPP/Uehling_Lab/NewGenomes/CombineLanes_AllReads/SAMandBAM_Polish"
METABAT_DIR="/nfs4/BPP/Uehling_Lab/NewGenomes/CombineLanes_AllReads/Metabat2Files_Polish"
QUAST_DIR="/nfs4/BPP/Uehling_Lab/NewGenomes/CombineLanes_AllReads/QUAST_Results_Polish"

mkdir -p "$ASSEMBLY_DIR" "$SAM_BAM_DIR" "$METABAT_DIR" "$QUAST_DIR"

# Initialize summary CSV
SUMMARY_CSV="$ASSEMBLY_DIR/assembly_summary.csv"
echo "Sample,Reads_Raw_Total,Reads_Trimmed_Total,PreTrim_Q20_Percent,PreTrim_Q30_Percent,PostTrim_Q20_Percent,PostTrim_Q30_Percent,Reads_Mapped_Percent,PrePolish_Contigs,PrePolish_Largest,PrePolish_Length,PrePolish_GC,PrePolish_N50,PrePolish_L50,PostPolish_Contigs,PostPolish_Largest,PostPolish_Length,PostPolish_GC,PostPolish_N50,PostPolish_L50,Bins_Created,Bin_Stats" > "$SUMMARY_CSV"

# Function to extract key QUAST metrics from report.txt
get_quast_data() {
    local report_file="$1/report.txt"
    if [[ ! -f "$report_file" ]]; then
        echo "0,0,0,0,0,0"
        return 1
    fi
    
    # Initialize default values
    local contigs=0 largest=0 length=0 gc=0 n50=0 l50=0
    
    # Read the report file and extract metrics
    while IFS= read -r line; do
        [[ -z "$line" ]] && continue
        [[ "$line" == *"statistics are based on"* ]] && continue
        
        # Extract metric and value using regex
        if [[ "$line" =~ ^\#\ contigs[[:space:]]+([0-9]+) ]]; then
            contigs="${BASH_REMATCH[1]}"
        elif [[ "$line" =~ ^Largest\ contig[[:space:]]+([0-9]+) ]]; then
            largest="${BASH_REMATCH[1]}"
        elif [[ "$line" =~ ^Total\ length[[:space:]]+([0-9]+) ]]; then
            length="${BASH_REMATCH[1]}"
        elif [[ "$line" =~ ^GC\ \(%\)[[:space:]]+([0-9]+\.[0-9]+) ]]; then
            gc="${BASH_REMATCH[1]}"
        elif [[ "$line" =~ ^N50[[:space:]]+([0-9]+) ]]; then
            n50="${BASH_REMATCH[1]}"
        elif [[ "$line" =~ ^L50[[:space:]]+([0-9]+) ]]; then
            l50="${BASH_REMATCH[1]}"
        fi
    done < "$report_file"
    
    # Return the metrics in the correct order
    echo "$contigs,$largest,$length,$gc,$n50,$l50"
}

# Process bins with QUAST metrics
process_bins() {
    local bin_dir=$1
    local sample_name=$2
    local sample_csv=$3
    
    echo "DEBUG: Checking bin directory: $bin_dir" >&2
    if [[ ! -d "$bin_dir" ]]; then
        echo "DEBUG: Bin directory $bin_dir does not exist" >&2
        echo "0|"
        return
    fi
    
    local bin_count=0
    local bin_stats=""
    local bin_summary_csv="$QUAST_DIR/${sample_name}_bin_summary.csv"
    echo "Bin,Contigs,Largest_Contig,Total_Length,GC_Percent,N50,L50" > "$bin_summary_csv"
    
    echo "DEBUG: Looking for bin files in $bin_dir/bin.*.fa" >&2
    for bin_file in "$bin_dir"/bin.*.fa; do
        if [[ ! -e "$bin_file" ]]; then
            echo "DEBUG: No bin files found in $bin_dir" >&2
            continue
        fi
        ((bin_count++))
        
        bin_num=$(basename "$bin_file" | sed 's/bin\.\([0-9]\+\).fa/\1/')
        bin_quast_dir="${QUAST_DIR}/${sample_name}_bin_${bin_num}_quast"
        
        echo "DEBUG: Processing bin $bin_file (Bin_$bin_num)" >&2
        # Run QUAST on the bin, redirecting all output to /dev/null
        quast.py -t $THREADS -o "$bin_quast_dir" "$bin_file" --silent > /dev/null 2>&1
        if [[ $? -ne 0 ]]; then
            echo "Warning: QUAST failed for bin $bin_num of sample $sample_name" >&2
            continue
        fi
        
        # Extract QUAST metrics
        IFS=',' read -r contigs largest length gc n50 l50 <<< "$(get_quast_data "$bin_quast_dir")"
        
        echo "DEBUG: QUAST metrics for Bin_$bin_num: contigs=$contigs, largest=$largest, length=$length, gc=$gc, n50=$n50, l50=$l50" >&2
        
        # Write to bin summary CSV
        echo "Bin_${bin_num},$contigs,$largest,$length,$gc,$n50,$l50" >> "$bin_summary_csv"
        
        # Write to sample CSV
        echo "$sample_name,Bin_${bin_num}_Contigs,$contigs" >> "$sample_csv"
        echo "$sample_name,Bin_${bin_num}_Largest_Contig,$largest" >> "$sample_csv"
        echo "$sample_name,Bin_${bin_num}_Total_Length,$length" >> "$sample_csv"
        echo "$sample_name,Bin_${bin_num}_GC_Percent,$gc" >> "$sample_csv"
        echo "$sample_name,Bin_${bin_num}_N50,$n50" >> "$sample_csv"
        echo "$sample_name,Bin_${bin_num}_L50,$l50" >> "$sample_csv"
        
        # Collect stats for overall summary
        bin_stats+="Bin_${bin_num}:Contigs=$contigs,Length=$length,N50=$n50;"
    done
    
    if [[ $bin_count -eq 0 ]]; then
        echo "DEBUG: No bins processed in $bin_dir" >&2
    else
        echo "DEBUG: Processed $bin_count bins" >&2
    fi
    
    bin_stats=${bin_stats%;}
    echo "${bin_count}|${bin_stats}"
}

# Main processing loop
for FWD_READS in $READS_DIR/*{_1.fastq.gz,_1P.fastq.gz,_R1.fastq.gz,_R1_001.fastq.gz}; do
    [ -e "$FWD_READS" ] || continue
    
    BASE_NAME=$(basename "$FWD_READS" | sed -E 's/_(1|1P|R1|R1_001)\.fastq\.gz$//')
    REV_READS=""
    for pattern in "_2.fastq.gz" "_2P.fastq.gz" "_R2.fastq.gz" "_R2_001.fastq.gz"; do
        [ -f "${READS_DIR}/${BASE_NAME}${pattern}" ] && REV_READS="${READS_DIR}/${BASE_NAME}${pattern}" && break
    done

    OUTPUT_DIR="${ASSEMBLY_DIR}/${BASE_NAME}_assembly"
    mkdir -p "$OUTPUT_DIR"
    SAMPLE_CSV="$OUTPUT_DIR/${BASE_NAME}_assembly_report.csv"

    echo "Sample,Metric,Value" > "$SAMPLE_CSV"

    # Trimming
    fastp -i "$FWD_READS" -I "$REV_READS" \
          -o "$OUTPUT_DIR/${BASE_NAME}_R1_trimmed.fastq.gz" \
          -O "$OUTPUT_DIR/${BASE_NAME}_R2_trimmed.fastq.gz" \
          -w $THREADS --detect_adapter_for_pe --qualified_quality_phred 20 \
          --trim_poly_g --correction \
          --length_required 50 --html "$OUTPUT_DIR/fastp_report.html" \
          --json "$OUTPUT_DIR/fastp_report.json"


    
    BEFORE_TRIMMING=$(jq -r '.summary.before_filtering.total_reads' "$OUTPUT_DIR/fastp_report.json")
    AFTER_TRIMMING=$(jq -r '.summary.after_filtering.total_reads' "$OUTPUT_DIR/fastp_report.json")
    PRETRIM_Q20_PERCENT=$(jq -r '.summary.before_filtering.q20_rate * 100' "$OUTPUT_DIR/fastp_report.json")
    PRETRIM_Q30_PERCENT=$(jq -r '.summary.before_filtering.q30_rate * 100' "$OUTPUT_DIR/fastp_report.json")
    POSTTRIM_Q20_PERCENT=$(jq -r '.summary.after_filtering.q20_rate * 100' "$OUTPUT_DIR/fastp_report.json")
    POSTTRIM_Q30_PERCENT=$(jq -r '.summary.after_filtering.q30_rate * 100' "$OUTPUT_DIR/fastp_report.json")
    echo "$BASE_NAME,Reads_Before_Trimming,$BEFORE_TRIMMING" >> "$SAMPLE_CSV"
    echo "$BASE_NAME,Reads_After_Trimming,$AFTER_TRIMMING" >> "$SAMPLE_CSV"
    echo "$BASE_NAME,PreTrim_Q20_Percent,$PRETRIM_Q20_PERCENT" >> "$SAMPLE_CSV"
    echo "$BASE_NAME,PreTrim_Q30_Percent,$PRETRIM_Q30_PERCENT" >> "$SAMPLE_CSV"
    echo "$BASE_NAME,PostTrim_Q20_Percent,$POSTTRIM_Q20_PERCENT" >> "$SAMPLE_CSV"
    echo "$BASE_NAME,PostTrim_Q30_Percent,$POSTTRIM_Q30_PERCENT" >> "$SAMPLE_CSV"

    # Assembly
    spades.py -t $THREADS -m 200 \
        -1 "$OUTPUT_DIR/${BASE_NAME}_R1_trimmed.fastq.gz" \
        -2 "$OUTPUT_DIR/${BASE_NAME}_R2_trimmed.fastq.gz" \
        -o "$OUTPUT_DIR/spades_output"
    
    [ -f "$OUTPUT_DIR/spades_output/contigs.fasta" ] && \
        mv "$OUTPUT_DIR/spades_output/contigs.fasta" "$OUTPUT_DIR/${BASE_NAME}_spades.fasta"

    # Pre-polish QUAST
    pre_polish_quast_dir="${QUAST_DIR}/${BASE_NAME}_pre_polish"
    quast.py -t $THREADS -o "$pre_polish_quast_dir" "$OUTPUT_DIR/${BASE_NAME}_spades.fasta" --silent > /dev/null 2>&1
    if [[ $? -ne 0 ]]; then
        echo "Warning: QUAST failed for pre-polish assembly of $BASE_NAME" >&2
    fi
    IFS=',' read -r pre_contigs pre_largest pre_length pre_gc pre_n50 pre_l50 <<< "$(get_quast_data "$pre_polish_quast_dir")"

    echo "$BASE_NAME,PrePolish_Contigs,$pre_contigs" >> "$SAMPLE_CSV"
    echo "$BASE_NAME,PrePolish_Largest_Contig,$pre_largest" >> "$SAMPLE_CSV"
    echo "$BASE_NAME,PrePolish_Total_Length,$pre_length" >> "$SAMPLE_CSV"
    echo "$BASE_NAME,PrePolish_GC_Percent,$pre_gc" >> "$SAMPLE_CSV"
    echo "$BASE_NAME,PrePolish_N50,$pre_n50" >> "$SAMPLE_CSV"
    echo "$BASE_NAME,PrePolish_L50,$pre_l50" >> "$SAMPLE_CSV"

    # Mapping
    bwa index "$OUTPUT_DIR/${BASE_NAME}_spades.fasta"
    bwa mem -t $THREADS "$OUTPUT_DIR/${BASE_NAME}_spades.fasta" \
        "$FWD_READS" "$REV_READS" | \
        samtools sort -@ $THREADS -o "$OUTPUT_DIR/${BASE_NAME}_aligned.bam"
    samtools index "$OUTPUT_DIR/${BASE_NAME}_aligned.bam"
    
    # Commenting out mapping percentage due to undefined get_mapping_stats
    # MAPPING_STATS=($(get_mapping_stats "$OUTPUT_DIR/${BASE_NAME}_aligned.bam" | tr ',' ' '))
    # MAPPED_PCT=$(echo "scale=2; ${MAPPING_STATS[1]}*100/${MAPPING_STATS[0]}" | bc)
    # echo "$BASE_NAME,Mapping_Percentage,$MAPPED_PCT" >> "$SAMPLE_CSV"
    echo "$BASE_NAME,Mapping_Percentage,NA" >> "$SAMPLE_CSV"

    # Polishing
    java -Xmx240g -jar /fs1/local/cqls/software/x86_64/pilon-1.24/envs/pilon/share/pilon-1.24-0/pilon.jar \
         --genome "$OUTPUT_DIR/${BASE_NAME}_spades.fasta" \
         --frags "$OUTPUT_DIR/${BASE_NAME}_aligned.bam" \
         --output "$OUTPUT_DIR/${BASE_NAME}_pilon" \
         --changes

    # Post-polish QUAST
    post_polish_quast_dir="${QUAST_DIR}/${BASE_NAME}_post_polish"
    quast.py -t $THREADS -o "$post_polish_quast_dir" "$OUTPUT_DIR/${BASE_NAME}_pilon.fasta" --silent > /dev/null 2>&1
    if [[ $? -ne 0 ]]; then
        echo "Warning: QUAST failed for post-polish assembly of $BASE_NAME" >&2
    fi
    IFS=',' read -r post_contigs post_largest post_length post_gc post_n50 post_l50 <<< "$(get_quast_data "$post_polish_quast_dir")"

    echo "$BASE_NAME,PostPolish_Contigs,$post_contigs" >> "$SAMPLE_CSV"
    echo "$BASE_NAME,PostPolish_Largest_Contig,$post_largest" >> "$SAMPLE_CSV"
    echo "$BASE_NAME,PostPolish_Total_Length,$post_length" >> "$SAMPLE_CSV"
    echo "$BASE_NAME,PostPolish_GC_Percent,$post_gc" >> "$SAMPLE_CSV"
    echo "$BASE_NAME,PostPolish_N50,$post_n50" >> "$SAMPLE_CSV"
    echo "$BASE_NAME,PostPolish_L50,$post_l50" >> "$SAMPLE_CSV"

    # Binning
    cd "$SAM_BAM_DIR"
    mkdir -p "$BASE_NAME"
    cd "$BASE_NAME"
    minimap2 -ax sr -t $THREADS "$OUTPUT_DIR/${BASE_NAME}_pilon.fasta" \
        "$FWD_READS" "$REV_READS" -o "${BASE_NAME}_SAM.sam"
    samtools view -@ $THREADS -b "${BASE_NAME}_SAM.sam" -o "${BASE_NAME}_BAM.bam"
    samtools sort -@ $THREADS -o "${BASE_NAME}_BAM.sorted.bam" "${BASE_NAME}_BAM.bam"
    samtools index -@ $THREADS "${BASE_NAME}_BAM.sorted.bam"
    
    cd "$METABAT_DIR"
    runMetaBat.sh -t $THREADS "$OUTPUT_DIR/${BASE_NAME}_pilon.fasta" "$SAM_BAM_DIR/$BASE_NAME/${BASE_NAME}_BAM.sorted.bam"
    
    # Find the actual MetaBAT output directory (handles timestamped directories)
    BIN_DIR=$(find "$METABAT_DIR" -maxdepth 1 -type d -name "${BASE_NAME}_pilon.fasta.metabat-bins*" | head -n 1)
    echo "DEBUG: MetaBAT bin directory set to: $BIN_DIR" >&2
    
    # Process bins
    if [[ -z "$BIN_DIR" ]]; then
        echo "DEBUG: No MetaBAT bin directory found for $BASE_NAME" >&2
        BIN_COUNT=0
        BIN_STATS=""
    else
        BIN_INFO=$(process_bins "$BIN_DIR" "$BASE_NAME" "$SAMPLE_CSV")
        BIN_COUNT=$(echo "$BIN_INFO" | cut -d'|' -f1)
        BIN_STATS=$(echo "$BIN_INFO" | cut -d'|' -f2)
    fi

    echo "$BASE_NAME,Bins_Created,$BIN_COUNT" >> "$SAMPLE_CSV"

    # Summary CSV
    echo "$BASE_NAME,$BEFORE_TRIMMING,$AFTER_TRIMMING,$PRETRIM_Q20_PERCENT,$PRETRIM_Q30_PERCENT,$POSTTRIM_Q20_PERCENT,$POSTTRIM_Q30_PERCENT,NA,$pre_contigs,$pre_largest,$pre_length,$pre_gc,$pre_n50,$pre_l50,$post_contigs,$post_largest,$post_length,$post_gc,$post_n50,$post_l50,$BIN_COUNT,\"$BIN_STATS\"" >> "$SUMMARY_CSV"
    #echo "$BASE_NAME,$BEFORE_TRIMMING,$AFTER_TRIMMING,$PRETRIM_Q20_PERCENT,$PRETRIM_Q30_PERCENT,$POSTTRIM_Q20_PERCENT,$POSTTRIM_Q30_PERCENT,NA,$pre_contigs,$pre_largest,$pre_length,$pre_gc,$pre_n50,$pre_l50" >> "$SUMMARY_CSV"

    echo "Sample $BASE_NAME completed!"
done

echo "Illumina assembly and reporting complete."