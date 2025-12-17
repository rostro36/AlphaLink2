#!/bin/bash
set -euo pipefail

# Default values for optional parameters
readonly DEFAULT_RECYCLING_ITERATIONS=20
readonly DEFAULT_NUMBER_OF_SAMPLES=25
readonly DEFAULT_NEFF=-1
readonly DEFAULT_DROPOUT_CROSSLINKS=-1
readonly DEFAULT_MODE="default"

# Function to display usage
show_usage() {
    cat << EOF
Protein Structure Prediction Pipeline

Usage: $(basename "$0") <fasta_path> <crosslinks> <output_dir_base> <param_path> 
                        <database_dir> <max_template_date> [options]

Required Arguments:
  fasta_path         Path to FASTA file containing protein sequence
  crosslinks         Path to crosslinking data file (CSV/TSV)
  output_dir_base    Base directory for output files
  param_path         Path to model parameters/weights file
  database_dir       Directory containing search databases
  max_template_date  Cutoff date for template search (YYYY-MM-DD)

Optional Arguments:
  recycling_iterations  Number of recycling iterations [default: $DEFAULT_RECYCLING_ITERATIONS]
  number_of_samples    Number of samples to generate [default: $DEFAULT_NUMBER_OF_SAMPLES]
  neff                 Effective number of sequences for subsampling (-1 for auto) [default: $DEFAULT_NEFF]
  dropout_crosslinks   Crosslink dropout rate (-1 for disabled) [default: $DEFAULT_DROPOUT_CROSSLINKS]
  mode                 Operation mode ("mmseqs", "precomputed", "default") [default: $DEFAULT_MODE]

Examples:
  # Minimal usage with defaults
  $0 /path/to/the/input.fasta /path/to/crosslinks.pkl.gz /path/to/the/output/directory/ /path/to/model_parameters.pt /path/to/database/directory/ 2020-05-01      
EOF
    exit 1
}

# Function to validate arguments
validate_arguments() {
    local missing_args=()
    
    # Check required arguments
    [[ -z "$fasta_path" ]] && missing_args+=("fasta_path")
    [[ -z "$crosslinks" ]] && missing_args+=("crosslinks")
    [[ -z "$output_dir_base" ]] && missing_args+=("output_dir_base")
    [[ -z "$param_path" ]] && missing_args+=("param_path")
    [[ -z "$database_dir" ]] && missing_args+=("database_dir")
    [[ -z "$max_template_date" ]] && missing_args+=("max_template_date")
    
    if [[ ${#missing_args[@]} -gt 0 ]]; then
        echo "ERROR: Missing required arguments: ${missing_args[*]}"
        show_usage
    fi
    
    # Validate file/directory existence
    [[ ! -f "$fasta_path" ]] && { echo "ERROR: FASTA file not found: $fasta_path"; exit 1; }
    [[ ! -f "$crosslinks" ]] && { echo "ERROR: Crosslinks file not found: $crosslinks"; exit 1; }
    [[ ! -f "$param_path" ]] && { echo "ERROR: Parameters file not found: $param_path"; exit 1; }
    [[ ! -d "$database_dir" ]] && { echo "ERROR: Database directory not found: $database_dir"; exit 1; }
    
    # Validate max_template_date format (basic check)
    if ! [[ "$max_template_date" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]]; then
        echo "WARNING: max_template_date '$max_template_date' may not be in YYYY-MM-DD format"
    fi
    
    # Validate numeric parameters
    [[ ! "$recycling_iterations" =~ ^[0-9]+$ ]] && { echo "ERROR: recycling_iterations must be a positive integer"; exit 1; }
    [[ ! "$number_of_samples" =~ ^[0-9]+$ ]] && { echo "ERROR: number_of_samples must be a positive integer"; exit 1; }
    
    # Create output directory if it doesn't exist
    mkdir -p "$output_dir_base"
}

# Function to log configuration
log_configuration() {
    cat << EOF
============================================
Pipeline Configuration
============================================
FASTA Path:          $fasta_path
Crosslinks File:     $crosslinks
Output Directory:    $output_dir_base
Parameters Path:     $param_path
Database Directory:  $database_dir
Max Template Date:   $max_template_date
Recycling Iterations: $recycling_iterations
Number of Samples:   $number_of_samples
Neff Parameter:      $neff
Dropout Crosslinks:  $dropout_crosslinks
Mode:                $mode
============================================
EOF
}

run_pipeline() {
    # Determine target name from FASTA file
    local fasta_file=$(basename "$fasta_path")
    local target_name="${fasta_file%.fa*}"
    
    # Create necessary directories
    create_output_directories
    
    # Generate or prepare MSAs based on mode
    generate_msa "$mode"
    
    # Run structure prediction
    run_prediction "$target_name"
}

# ============================================
# Helper Functions
# ============================================

create_output_directories() {
    echo "Creating output directories in: $output_dir_base"
    
    local dirs=(
        "$output_dir_base/msas"
        "$output_dir_base/predictions"
        "$output_dir_base/logs"
        "$output_dir_base/tmp"
    )
    
    for dir in "${dirs[@]}"; do
        mkdir -p "$dir"
        if [[ $? -ne 0 ]]; then
            echo "ERROR: Failed to create directory: $dir"
            exit 1
        fi
    done
}

# ============================================
# MSA Generation Functions
# ============================================

generate_msa() {
    local mode="$1"
    
    case "$mode" in
        "default")
            generate_default_msa
            ;;
        "precomputed")
            prepare_precomputed_msa
            ;;
        "mmseqs")
            generate_mmseqs_msa
            ;;
        *)
            echo "ERROR: Unknown mode '$mode'. Valid modes: default, precomputed, mmseqs"
            exit 1
            ;;
    esac
}

generate_default_msa() {
    echo "Starting default MSA generation..."
    
    # Define database paths
    local databases=(
        "--uniref90_database_path=$database_dir/uniref90/uniref90.fasta"
        "--mgnify_database_path=$database_dir/mgnify/mgy_clusters_2018_12.fa"
        "--bfd_database_path=$database_dir/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
        "--uniclust30_database_path=$database_dir/uniclust30/uniclust30_2018_08/uniclust30_2018_08"
        "--uniprot_database_path=$database_dir/uniprot/uniprot.fasta"
        "--pdb_seqres_database_path=$database_dir/pdb_seqres/pdb_seqres.txt"
        "--template_mmcif_dir=$database_dir/pdb_mmcif/mmcif_files"
    )
    
    # Check if required databases exist
    for db_arg in "${databases[@]}"; do
        local db_path="${db_arg#*=}"
        if [[ ! -e "$db_path" ]]; then
            echo "WARNING: Database file not found: $db_path"
        fi
    done
    
    # Run MSA generation
    local cmd=(
        python unifold/homo_search.py
        --fasta_path="$fasta_path"
        --max_template_date="$max_template_date"
        --output_dir="$output_dir_base"
        --obsolete_pdbs_path="$database_dir/pdb_mmcif/obsolete.dat"
        --use_precomputed_msas=True
    )
    
    # Add database arguments
    for db_arg in "${databases[@]}"; do
        cmd+=("$db_arg")
    done
    
    echo "Running: ${cmd[*]}"
    if ! "${cmd[@]}"; then
        echo "ERROR: Default MSA generation failed"
        exit 1
    fi
}

prepare_precomputed_msa() {
    echo "Preparing precomputed MSA..."
    
    local cmd=(
        python scripts/prepare_precomputed.py
        --input_path="$fasta_path"
        --output_dir="$output_dir_base"
    )
    
    echo "Running: ${cmd[*]}"
    if ! "${cmd[@]}"; then
        echo "ERROR: MMseqs2 MSA generation failed"
        exit 1
    fi
}

generate_mmseqs_msa() {
    echo "Starting MMseqs2 server MSA generation..."
    
    local cmd=(
        python scripts/prepare_mmseqs.py
        --fasta_path="$fasta_path"
        --max_template_date="$max_template_date"
        --output_dir="$output_dir_base"
    )
    
    echo "Running: ${cmd[*]}"
    if ! "${cmd[@]}"; then
        echo "ERROR: MMseqs2 MSA generation failed"
        exit 1
    fi
}

# ============================================
# Structure Prediction Function
# ============================================

run_prediction() {
    local target_name="$1"
    
    echo "Starting structure prediction for: $target_name"
    
    # Check if data directory exists after MSA generation
    if [[ ! -d "$output_dir_base" ]]; then
        echo "ERROR: Output directory not found: $output_dir_base"
        exit 1
    fi
    
    # Check if crosslinks file exists (if provided and not -1)
    if [[ "$dropout_crosslinks" != "-1" ]] && [[ ! -f "$crosslinks" ]]; then
        echo "WARNING: Crosslinks file not found: $crosslinks"
        echo "Continuing without crosslink constraints..."
        local crosslinks_arg=""
    else
        if [[ "$crosslinks" == *.csv ]]; then
            python generate_crosslink_pickle.py \
            --csv "$crosslinks" \
            --output "${crosslinks%.csv}.pkl.gz"
        fi
        local crosslinks_arg="--crosslinks=$crosslinks"
    fi
    
    # Build prediction command
    local cmd=(
        python inference.py
        --model_name="multimer_af2_crop"
        --param_path="$param_path"
        --data_dir="$output_dir_base"
        --target_name="$target_name"
        --output_dir="$output_dir_base"
        --bf16
        --use_uniprot
        --max_recycling_iters="$recycling_iterations"
        --times="$number_of_samples"
        --neff="$neff"
        --save_raw_output
        --relax
    )
    
    # Add crosslinks argument if provided
    [[ -n "$crosslinks_arg" ]] && cmd+=("$crosslinks_arg")
    
    # Add dropout_crosslinks if not disabled
    if [[ "$dropout_crosslinks" != "-1" ]]; then
        cmd+=(--dropout_crosslinks="$dropout_crosslinks")
    fi
    
    echo "Running prediction command..."
    echo "Command: ${cmd[*]}"
    
    # Run prediction with logging
    local log_file="$output_dir_base/logs/prediction_$(date +%Y%m%d_%H%M%S).log"
    echo "Logging to: $log_file"
    
    if ! "${cmd[@]}" 2>&1 | tee "$log_file"; then
        echo "ERROR: Structure prediction failed"
        echo "Check log file: $log_file"
        exit 1
    fi
    
    echo "Prediction completed successfully"
    echo "Results saved in: $output_dir_base/predictions"
    echo "Log file: $log_file"
}

# Main script execution
main() {
    # Parse command line arguments
    fasta_path="${1:-}"
    crosslinks="${2:-}"
    output_dir_base="${3:-}"
    param_path="${4:-}"
    database_dir="${5:-}"
    max_template_date="${6:-}"
    
    # Optional arguments with defaults
    recycling_iterations="${7:-$DEFAULT_RECYCLING_ITERATIONS}"
    number_of_samples="${8:-$DEFAULT_NUMBER_OF_SAMPLES}"
    neff="${9:-$DEFAULT_NEFF}"
    dropout_crosslinks="${10:-$DEFAULT_DROPOUT_CROSSLINKS}"
    mode="${11:-$DEFAULT_MODE}"
    
    # Display help if requested
    if [[ "$#" -eq 1 ]] && [[ "$1" == "--help" || "$1" == "-h" ]]; then
        show_usage
    fi
    
    # Validate arguments
    validate_arguments
    
    # Log configuration
    log_configuration
    
    # Add your pipeline execution code here
    run_pipeline
}

# Execute main function with all arguments
main "$@"
