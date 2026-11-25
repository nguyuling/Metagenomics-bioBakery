#!/bin/bash
set -e

# SET output and database directories
OUT_DIR=/Volumes/T7/metagenomics            # output directory
METAPHLAN_DB=/Users/nguyuling/Documents/metaphlan_databases # metaphlan database
HUMANN_DB=/Volumes/T7/conda_envs/biobakery/lib/python3.10/site-packages/humann/data
LOG="$OUT_DIR/logs"                         # log directories
threads=8                                   # number of core

# output directories
READS="$OUT_DIR/0_Reads"
FASTQC="$OUT_DIR/1_Fastqc"
MULTIQC="$OUT_DIR/2_Multiqc"
CUTADAPT="$OUT_DIR/3_Cutadapt"
METAPHLAN="$OUT_DIR/4_MetaPhlAn"
HUMANN="$OUT_DIR/5_HUMAnN"
GENES_MAPPING="$OUT_DIR/6_Genes_Mapping"
PATHWAYS_ANALYSIS="$OUT_DIR/7_Pathways_Analysis"
CROSS_ANALYSIS="$OUT_DIR/8_Cross_Analysis"

# subdirectories
Genes_Abs="$GENES_MAPPING/Genes_Abs"
Genes_Relab="$GENES_MAPPING/Genes_Relab"
Unpacked_Pathabundance="$PATHWAYS_ANALYSIS/Unpacked_Pathabundance"
Unpacked_Pathcoverage="$PATHWAYS_ANALYSIS/Unpacked_Pathcoverage"

# output files
cutadapt_txt="$CUTADAPT/cutadapt.txt"

# mapping variable
CATEGORIES="uniref90 uniref50 eggnog go ko level4ec pfam"

mkdir -p "$LOG"
log_file="$LOG/pipeline_$(date '+%Y%m%d_%H%M%S').log"

log() {                                     # log datetime and message to terminal and log file
    msg="[$(date '+%Y-%m-%d %H:%M:%S')] $*"
    echo "$msg"
    echo "$msg" >> "$log_file"
}

run() {                                     # usage: run <command> <parameter>
    log "START: $*"                         # log the starting step's command
    start=$(date +%s)                       # store the start time
    if "$@" >> "$log_file" 2>&1; then       # run the command, write stdout / stderr to log file
        ran_count=$(( ${ran_count:-0} + 1 ))    # add 1 for each successful step run
        log "OK ($(( $(date +%s) - start ))s)"  # exit 0 : log run time (current time - start time)
    else
        rc=$?                               # exit != 0 : hold exit code of the last command
        log "FAILED (exit $rc): $*"         # log the exit code and the failed command
        return "$rc"                        # return the exit code
    fi
}

need() {                                    # usage: need <file or dir> <step number that creates it>
    if [ ! -e "$1" ]; then                  # if the prerequisite file is absent
        log "ERROR: missing $1. Run step $2 first."  && return 1   # log reminder message about running the upstream step first
    fi
}

skip_if_done() {                            # usage: skip_if_done <step name> <output>...
    name=$1; shift                          # remove first arguments (step name), leaving only the output files
    [ "${force:-0}" = 1 ] && { log "RUN $name: force=1"; return 1; }
    for f in "$@"; do                       # for each file
        if [ ! -e "$f" ]; then              # if the file is missing
            log "RUN $name: missing $f"     # log message to run the step
            return 1                        # stop the function
        fi
    done
    log "SKIPPED $name: outputs already exist"  # if file(s) exist, log skipping step message
    return 0
}

run_steps() {                               # usage: run_stages <force 0|1> <function>...
    force=$1; shift                         
    start=$(date +%s)                       # store the start time 
    for step in "$@"; do                    # for each step
        if ! "$step"; then                  # run the step, if it fails
            log "STOPPED at $step. Fix the error and run again."    # log failure message
            return 1                        # stop this function
        fi
    done
    total=$(( $(date +%s) - start ))        # if successful, count the total run time
    log "DONE. Total time: $(( total / 60 )) min $(( total % 60 )) sec" # log total run time
}