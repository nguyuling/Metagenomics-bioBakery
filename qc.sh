#!/bin/bash
set -e

# bash qc.sh all        run all unfinished steps
# bash qc.sh allforce   run all steps regardlessly

# helper function script reference
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/lib.sh"

move_reads() {                          # to move all raw reads from outdir to subdir
    mkdir -p "$READS"
    if ls "$OUT_DIR"/*.fastq* >/dev/null 2>&1; then  # if read files exist in outdir
        mv "$OUT_DIR"/*.fastq* "$READS" || return 1  # move them to reads folder
        log "Moved reads to $READS"
    else
        log "No .fastq* files in $OUT_DIR, skipping move."   # to handle subsequent runs
    fi
}

run_fastqc() {                          # to run fastq on all read pairs
    need "$READS" 1 || return 1
    run mkdir -p "$FASTQC" || return 1
    for f in "$READS"/*.fastq*; do      # for each read file                  
        [ -f "$f" ] || continue         # if it's a file
        base=$(basename "$f")           # extract filename from full path
        name=${base%.fastq.gz}          # strip .fastq.gz extension
        name=${name%.fastq}             # strip .fastq extension
        skip_if_done "fastqc $base" "$FASTQC/${name}_fastqc.html" && continue
        run fastqc -o "$FASTQC" -t "$threads" "$f" || return 1
    done
}

run_multiqc() {                         # to run multiqc on the fastqc folder 
    need "$FASTQC" 1 || return 1
    skip_if_done "multiqc" "$MULTIQC/multiqc_report.pdf" && return 0
    mkdir -p "$MULTIQC"
    run multiqc \
        -n multiqc_report \
        -o "$MULTIQC" \
        "$FASTQC" || return 1
}

make_cutadapt_txt() {                   # to make default cutadapt specification txt allowing custom flag/value
    need "$READS" 1 || return 1
    skip_if_done "make_cutadapt_txt" "$cutadapt_txt" && return 0
    mkdir -p "$(dirname "$cutadapt_txt")"
    mkdir -p "$CUTADAPT"

    # mapping table: column_header     value
        mapping='sample_name    "$sample_name"
    r1    "$r1"
    r2    "$r2"
    -a    AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
    -A    AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
    -q    20
    -m    20'

    header=$(printf '%s\n' "$mapping" | while read -r header value; do  # read header and value columns from mapping
        printf '%s\t' "$header"         # append headers separated by tab
    done)
    printf '%s\n' "${header%\t}" > "$cutadapt_txt"  # write header to cutadapt.txt

    paired=0                            # init counters
    unpaired=0

    for f in "$READS"/*_R1*; do         # for each read file containing _R1
        [ -f "$f" ] || continue         # check if it is file
        r1=$(basename "$f")             # store r1 file name
        sample_name="${r1%%_R1*}"       # sample name = text before _R1
        r2=$(echo "$r1" | sed 's/_R1/_R2/') # generate expected r2 filename
        if [ ! -f "$READS/$r2" ]; then  # if expected r2 not found
            unpaired=$((unpaired + 1))  # increase unpaired counter
            continue                    # continue to find pair reads
        fi
        row=$(printf '%s\n' "$mapping" | while read -r header value; do # read header and value column from mapping
            eval "v=$value"             # expand value (to handle reference variables)
            printf '%s\t' "$v"          # append values separated by tab
        done)
        printf '%s\n' "${row%$\t}" >> "$cutadapt_txt"   # write values to cutadapt.txt
        paired=$((paired + 1))          # increase paired counter
    done

    for f in "$READS"/*_R2*; do         # for each read file containing _R2
        [ -f "$f" ] || continue         # check if it is file
        r2=$(basename "$f")             # store r2 file name
        r1=$(echo "$r2" | sed 's/_R2/_R1/')     # generate expected r1 filename
        [ -f "$READS/$r1" ] || unpaired=$((unpaired + 1))   # if expected r1 not found then increase unpaired counter
    done

    log "$paired pair(s) written to cutadapt.txt, edit the flag (-a, -A, -q, -m) if needed."
    if [ "$unpaired" -gt 0 ]; then
        log "WARNING: $unpaired file(s) have no matching pair and are NOT included in cutadapt.txt."
        log "Please add their sample_name, r1 and r2 to cutadapt.txt manually."
    fi
}

trim_adapters() {                       # to trim adapters of each read pair
    need "$cutadapt_txt" 1 || return 1

    header=$(head -n 1 "$cutadapt_txt")
    row=$(awk -F'\t' -v s="$sample_name" '$1 == s {print; exit}' "$cutadapt_txt")
    if [ -z "$row" ]; then
        log "ERROR: $sample_name not found in cutadapt.txt" && return 1
    fi

    IFS=$'\t' read -r -a cols <<< "$header"
    IFS=$'\t' read -r -a vals <<< "$row"
    r1_name=""
    r2_name=""
    cut_args=()
    for i in "${!cols[@]}"; do
        col="${cols[$i]}"
        val="${vals[$i]:-}"
        case "$col" in
            ""|sample_name) continue ;;
            r1)             r1_name="$val" ;;
            r2)             r2_name="$val" ;;
            -o|-p|-j)       log "WARN: column $col is set by the script, ignored" ;;
            -*)             [ -n "$val" ] && cut_args+=("$col" "$val") ;;
            *)              log "WARN: unknown column $col, ignored" ;;
        esac
    done

    if [ -z "$r1_name" ] || [ -z "$r2_name" ]; then
        log "ERROR: r1 and r2 columns must be set for $sample_name" && return 1
    fi

    r1="$READS/$r1_name"
    r2="$READS/$r2_name"
    need "$r1" 1 || return 1
    need "$r2" 1 || return 1

    skip_if_done "trim_adapters" "$trimmed_r1" "$trimmed_r2" && return 0
    mkdir -p "$sample_cutadapt"

    run cutadapt "${cut_args[@]}" \
        -j "$threads" \
        -o "$trimmed_r1" \
        -p "$trimmed_r2" \
        --report=full \
        "$r1" "$r2"
}

filter_rrna() {                         # to separate rRNA and non-rRNA reads
    need "$trimmed_r1" 2 || return 1
    need "$trimmed_r2" 2 || return 1
    need "$SORTMERNA_DB" || return 1
    skip_if_done "filter_rrna" "$unaligned" && return 0
    mkdir -p "$sample_sortmerna"
    run sortmerna \
        --ref "$SORTMERNA_DB/silva-arc-16s-id95.fasta" \
        --ref "$SORTMERNA_DB/silva-arc-23s-id98.fasta" \
        --ref "$SORTMERNA_DB/silva-bac-16s-id90.fasta" \
        --ref "$SORTMERNA_DB/silva-bac-23s-id98.fasta" \
        --ref "$SORTMERNA_DB/silva-euk-18s-id95.fasta" \
        --ref "$SORTMERNA_DB/silva-euk-28s-id98.fasta" \
        --ref "$SORTMERNA_DB/rfam-5s-database-id98.fasta" \
        --ref "$SORTMERNA_DB/rfam-5.8s-database-id98.fasta" \
        --workdir "$sample_sortmerna" \
        --aligned "$sample_sortmerna/${sample_name}_aligned" \
        --other "$sample_sortmerna/${sample_name}_unaligned" \
        --fastx \
        --out2 \
        --paired_out \
        --threads "$threads" \
        --reads "$trimmed_r1" \
        --reads "$trimmed_r2" || return 1
    cat \
        "$sample_sortmerna/${sample_name}_unaligned_fwd.fq" \
        "$sample_sortmerna/${sample_name}_unaligned_rev.fq" \
        > "$unaligned"
}

sample_done() {                         # usage: sample_done <sample_name>
    [ -e "$CUTADAPT/$1/${1}_trimmed_R1.fastq" ] &&
    [ -e "$CUTADAPT/$1/${1}_trimmed_R2.fastq" ] &&
    [ -e "$SORTMERNA/$1/${1}.fq" ]
}

trim_filter() {                         # to ask user before looping over pairs to trim and filter
    samples=$(awk -F'\t' 'NR > 1 && $1 != "" {print $1}' "$cutadapt_txt")

    pending=0
    for s in $samples; do
        [ -e "$CUTADAPT/$s/${s}_trimmed_R1.fastq" ] &&
        [ -e "$CUTADAPT/$s/${s}_trimmed_R2.fastq" ] &&
        [ -e "$SORTMERNA/$s/${s}.fq" ] || pending=1
    done
    if [ "$pending" = 0 ] && [ "$force" != 1 ]; then
        log "Trimming and rRNA filtering already done for all samples"
        return 0
    fi

    printf 'Continue with trim_adapters and filter_rrna using the current cutadapt.txt? [y/n] '
    read -r answer || answer=n          # no input counts as no
    case "$answer" in
        y|Y) ;;                         # y input runs trim & filter
        *)   log "Stopped before trimming. Edit cutadapt.txt, then run again."; return 0 ;;
    esac

    failed=""
    for sample_name in $samples; do
        log "=== Sample $sample_name ==="
        sample_cutadapt="$CUTADAPT/$sample_name"
        sample_sortmerna="$SORTMERNA/$sample_name"
        trimmed_r1="$sample_cutadapt/${sample_name}_trimmed_R1.fastq"
        trimmed_r2="$sample_cutadapt/${sample_name}_trimmed_R2.fastq"
        unaligned="$sample_sortmerna/${sample_name}.fq"

        { trim_adapters && filter_rrna; } || failed="$failed $sample_name"
    done
    if [ -n "$failed" ]; then
        log "Failed samples:$failed"
        return 1
    fi
}

case "${1:-all}" in
    all)      run_stages 0 move_reads run_fastqc run_multiqc make_cutadapt_txt trim_filter ;;
    allforce) run_stages 1 move_reads run_fastqc run_multiqc make_cutadapt_txt trim_filter ;;
    *)        echo "Usage: $0 [all|allforce]"; exit 1 ;;
esac