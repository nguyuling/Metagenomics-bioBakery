#!/bin/bash
set -e

# bash biobakery.sh all        run all unfinished steps
# bash biobakery.sh allforce   run all steps regardlessly

# helper function script reference
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/lib.sh"

run_pair() {                                # to run metaphlan & humann on a read pair
    sample_name="$1"

    # per-sample output directories
    Sortmerna="$SORTMERNA/$sample_name"
    Metaphlan="$METAPHLAN/$sample_name"
    Humann="$HUMANN/$sample_name"

    # per-sample output files
    trimmed="$CUTADAPT/${sample_name}_trimmed.fastq"
    taxonomic_profile="$Metaphlan/${sample_name}_taxonomic_profile.tsv"
    genefamilies="$Humann/${sample_name}_genefamilies.tsv"
    pathabundance="$Humann/${sample_name}_pathabundance.tsv"
    pathcoverage="$Humann/${sample_name}_pathcoverage.tsv"
    genes_relab_uniref90="$Genes_Relab/${sample_name}_uniref90.tsv"
    pathabundance_renormed="$Unpacked_Pathabundance/${sample_name}.tsv.tmp"
    pathabundance_unpacked="$Unpacked_Pathabundance/${sample_name}.tsv"
    pathcoverage_unpacked="$Unpacked_Pathcoverage/${sample_name}.tsv"

    log "=== Sample $sample_name ==="
    for stage in taxonomic_profiling functional_profiling genes_mapping unpack_pathways; do
        "$stage" || { log "Sample $sample_name stopped at $stage"; return 1; }
    done
}

run_all_pairs() {                           # to loop run_pair() on all samples
    force=$1
    before=${ran_count:-0}
    failed=""

    for s in $(awk -F'\t' 'NR > 1 && $1 != "" {print $1}' "$cutadapt_txt"); do
        run_pair "$s" || failed="$failed $s"
    done

    if [ -n "$failed" ]; then
        log "Step 8 skipped: failed samples:$failed" && return 1
    fi

    if [ "${ran_count:-0}" -gt "$before" ] || [ "$force" = 1 ] || [ ! -e "$CROSS_ANALYSIS/unpacked_pathcoverage.tsv" ]; then
        cross_analysis || { log "Cross-analysis failed"; return 1; }
    else
        log "Step 8 skipped: no upstream changes"
    fi
}

taxonomic_profiling() {                     # to run metaphlan and generate krona pie chart 
    need "$trimmed" 3 || return 1
    need "$METAPHLAN_DB" || return 1
    skip_if_done "taxonomic_profiling" "$taxonomic_profile" "$Metaphlan/$sample_name.krona.html" && return 0
    mkdir -p "$Metaphlan"
    run metaphlan \
        --input_type fastq \
        --bowtie2db "$METAPHLAN_DB" \
        -x mpa_vJun23_CHOCOPhlAnSGB_202307 \
        --bowtie2out "$Metaphlan/${sample_name}_metagenome_vJun23.bowtie2.bz2" \
        -o "$taxonomic_profile" \
        --nproc "$threads" \
        --offline \
        "$trimmed"  || return 1
    run metaphlan2krona.py \
        -p "$taxonomic_profile" \
        -k "$Metaphlan/${sample_name}_krona.txt"  || return 1
    run ktImportText \
        -o "$Metaphlan/$sample_name.krona.html" \
        -n 'all' \
        "$Metaphlan/${sample_name}_krona.txt"
}

functional_profiling() {                    # to run humann 
    need "$trimmed" 3 || return 1
    need "$taxonomic_profile" 4 || return 1
    need "$HUMANN_DB/chocophlan" || return 1
    need "$HUMANN_DB/uniref" || return 1
    skip_if_done "functional_profiling" "$genefamilies" "$pathabundance" "$pathcoverage" && return 0
    mkdir -p "$Humann"
    run humann \
        -i "$trimmed" \
        -o "$Humann" \
        --taxonomic-profile "$taxonomic_profile" \
        --nucleotide-database "$HUMANN_DB/chocophlan" \
        --protein-database "$HUMANN_DB/uniref" \
        --memory-use maximum \
        --threads "$threads"
}

genes_mapping() {                           # to regroup, sort, rename and renorm gene families abundance
    need "$genefamilies" 5 || return 1
    need "$HUMANN_DB/utility_mapping" || return 1
    output=()
    for c in $CATEGORIES; do
        output+=("$Genes_Abs/${sample_name}_${c}.tsv" "$Genes_Relab/${sample_name}_${c}.tsv")
    done
    skip_if_done "genes_mapping" "${output[@]}" && return 0
    mkdir -p "$Genes_Abs" "$Genes_Relab"

    while read -r category groupmap namemap; do
        regrouped="$GENES_MAPPING/${sample_name}_${category}.tsv.tmp"
        renamed="$Genes_Abs/${sample_name}_${category}.tsv"
        renormed="$Genes_Relab/${sample_name}_${category}.tsv"
        if [ "$groupmap" = "none" ]; then
            cp "$genefamilies" "$regrouped"  || return 1
        else
            run humann_regroup_table \
                -i "$genefamilies" \
                -c "$HUMANN_DB/utility_mapping/$groupmap" \
                -o "$regrouped" || return 1
        fi
        run humann_rename_table \
            -i "$regrouped" \
            -c "$HUMANN_DB/utility_mapping/$namemap" \
            -o "$renamed" || return 1
        run humann_renorm_table \
            -i "$regrouped" \
            -u relab \
            -o "$renormed" || return 1
        rm -f "$regrouped"
    done <<EOF
uniref90    none                            map_uniref90_name.txt.bz2
uniref50    map_uniref50_uniref90.txt.gz    map_uniref50_name.txt.bz2
eggnog      map_eggnog_uniref90.txt.gz      map_eggnog_name.txt.gz
go          map_go_uniref90.txt.gz          map_go_name.txt.gz
ko          map_ko_uniref90.txt.gz          map_ko_name.txt.gz
level4ec    map_level4ec_uniref90.txt.gz    map_ec_name.txt.gz
pfam        map_pfam_uniref90.txt.gz        map_pfam_name.txt.gz
EOF
}

unpack_pathways() {                         # to assign gene families to pathways
    need "$genes_relab_uniref90" 6 || return 1
    need "$pathabundance" 5 || return 1
    need "$pathcoverage" 5 || return 1
    need "$HUMANN_DB/pathways" || return 1
    skip_if_done "unpack_pathways" "$pathabundance_unpacked" "$pathcoverage_unpacked" && return 0
    rm -f "$pathabundance_renormed"
    run humann_renorm_table \
        -i "$pathabundance" \
        -u relab \
        -o "$pathabundance_renormed" || return 1
    run humann_unpack_pathways \
        --input-genes "$genes_relab_uniref90" \
        --input-pathways "$pathabundance_renormed" \
        --gene-mapping "$HUMANN_DB/pathways/metacyc_reactions_level4ec_only.uniref.bz2" \
        --pathway-mapping "$HUMANN_DB/pathways/metacyc_pathways_structured_filtered_v24" \
        --output "$pathabundance_unpacked" || return 1
    rm -f "$pathabundance_renormed"         # rm temp renormed pathways abundance after unpacking
    run humann_unpack_pathways \
        --input-genes "$genes_relab_uniref90" \
        --input-pathways "$pathcoverage" \
        --gene-mapping "$HUMANN_DB/pathways/metacyc_reactions_level4ec_only.uniref.bz2" \
        --pathway-mapping "$HUMANN_DB/pathways/metacyc_pathways_structured_filtered_v24" \
        --output "$pathcoverage_unpacked" || return 1
}

cross_analysis() {                          # to join gene families & pathways tables of all samples
    output=()
    for c in $CATEGORIES; do
        output+=("$CROSS_ANALYSIS/Genes_Abs_${c}.tsv" "$CROSS_ANALYSIS/Genes_Relab_${c}.tsv")
    done
    skip_if_done "cross_analysis" "${output[@]}" "$CROSS_ANALYSIS/Unpacked_Pathabundance.tsv" "$CROSS_ANALYSIS/Unpacked_Pathcoverage.tsv" && return 0
    mkdir -p "$CROSS_ANALYSIS"

    for c in $CATEGORIES; do
        genes_abs="$CROSS_ANALYSIS/Genes_Abs_${c}.tsv.tmp"
        genes_relab="$CROSS_ANALYSIS/Genes_Relab_${c}.tsv.tmp"
        run humann_join_tables \
            -i "$Genes_Abs" \
            --file_name "_${c}.tsv" \
            -o "$genes_abs" || return 1
        run humann_join_tables \
            -i "$Genes_Relab" \
            --file_name "_${c}.tsv" \
            -o "$genes_relab" || return 1
        run python "$SCRIPT_DIR/utils/merge_genes.py" \
            "$CROSS_ANALYSIS/genes_${c}.tsv" "$genes_abs" "$genes_relab" "$c" || return 1
        rm -f "$genes_abs" "$genes_relab"
    done

    find "$Unpacked_Pathabundance" -name "*.tsv.tmp" -delete
    run humann_join_tables \
        -i "$Unpacked_Pathabundance" \
        -o "$CROSS_ANALYSIS/unpacked_pathabundance.tsv" || return 1
    run humann_join_tables \
        -i "$Unpacked_Pathcoverage" \
        -o "$CROSS_ANALYSIS/unpacked_pathcoverage.tsv" || return 1
    run python "$SCRIPT_DIR/utils/make_excel.py" "genes_pathways" "$OUT_DIR" "$CROSS_ANALYSIS" || return 1
}

# entry point
case "${1:-}" in
    all)        run_all_pairs 0 ;;
    allforce)   run_all_pairs 1 ;;
    *)          echo "Usage: bash $0 [ all  |  allforce ]"; exit 1 ;;
esac