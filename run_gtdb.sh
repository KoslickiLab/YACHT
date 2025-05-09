NUM_THREADS=1 # Adjust based on your machine's capabilities
cd demorun

# Prepare results directories
# rm -rf ./results/*
mkdir -p ./results/results_orig
mkdir -p ./results/results_mod

# Output file for logging
# runtime_log="runtime_log_gtdb1000_mod.txt"
runtime_log="runtime_log_orig_sample6-20.txt"
> $runtime_log  # Clear previous log file

# Initialize total time and sample count for calculating running averages
total_time_run=0
sample_count=0

# Function to convert seconds to HH:MM:SS format
format_time() {
    local total_seconds=$1
    local hours=$((total_seconds / 3600))
    local minutes=$(((total_seconds % 3600) / 60))
    local seconds=$((total_seconds % 60))
    printf "%02d:%02d:%02d" $hours $minutes $seconds
}

# Function to calculate running average
calculate_avg() {
    local total_time=$1
    local count=$2
    local avg_time=$((total_time / count))
    echo $(format_time $avg_time)
}

# yacht train --ref_file /scratch/akn5655/YACHT/demorun/gtdb-rs214-reps.k31_0.995_pretrained/1000sig.zip --ksize 31 --num_threads ${NUM_THREADS} --ani_thresh 0.995 --prefix 'gtdb1000_ani_thresh_0.995' --outdir ./ --force

# Start the timer
start_total=$(date +%s)
# /scratch/akn5655/YACHT/demorun/gtdb-rs214-reps.k31_0.995_pretrained/gtdb-rs214-reps.k31_0.995_config.json
# /scratch/akn5655/566/YACHT/demorun/gtdb1000/gtdb1000_ani_thresh_0.995_config.json

# Loop through samples 0 to 20
for i in {6..20}; do
    sample_count=$((sample_count + 1))
    # Yacht run per-sample timing
    start_run=$(date +%s)
    yacht run \
        --json "/scratch/akn5655/YACHT/demorun/gtdb-rs214-reps.k31_0.995_pretrained/gtdb-rs214-reps.k31_0.995_config.json" \
        --sample_file "/scratch/akn5655/YACHT-reproducibles/benchmark/CAMI_data/rhizosphere_data/sketches/rhimgCAMI2_sample_${i}.sig.zip" \
        --significance 0.99 \
        --num_threads ${NUM_THREADS} \
        --min_coverage_list 0.1 \
        --out "rhimgCAMI2_sample_${i}.sig.zip_result.xlsx"
    end_run=$(date +%s)

    per_sample_runtime_run=$((end_run - start_run))
    total_time_run=$((total_time_run + per_sample_runtime_run))
    avg_run=$(calculate_avg $total_time_run $sample_count)
    formatted_run_time=$(format_time $per_sample_runtime_run)

    # Log per-sample time and running average for yacht run
    echo "Sample $i (YACHT run): $formatted_run_time, Running Avg: $avg_run" | tee -a $runtime_log

    # Move the result
    mv "rhimgCAMI2_sample_${i}.sig.zip_result.xlsx" ./results/results_orig/
done

end_total=$(date +%s)

# Calculate total runtimes
formatted_total_time_run=$(format_time $total_time_run)


# Print and log total runtimes
echo -e "\nTotal execution time for YACHT run: $formatted_total_time_run" | tee -a $runtime_log

