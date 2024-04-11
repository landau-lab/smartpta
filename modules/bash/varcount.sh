#!/bin/bash

# Check for correct number of arguments
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <VCF file>"
    exit 1
fi

# Create a temporary directory for intermediate files
temp_dir=$(mktemp -d -t varcount-XXXXXXXXXX)
trap 'rm -rf -- "$temp_dir"' EXIT

# Create header with sample names
bcftools query -l "$1" | awk '{print $0}' ORS='\t' | sed 's/\t$//' > "$temp_dir/header.tsv"
echo -e "CHROM_POS_REF_ALT\t$(cat "$temp_dir/header.tsv")" > "$temp_dir/header_formatted.tsv"

# Prepare output files
vcf_basename=$(basename "$1" | cut -d. -f1)
total_depth_matrix="${vcf_basename}_NR.tsv"
variant_support_matrix="${vcf_basename}_NV.tsv"

# Initialize output files with headers
cat "$temp_dir/header_formatted.tsv" > "$total_depth_matrix"
cat "$temp_dir/header_formatted.tsv" > "$variant_support_matrix"

# Stream bcftools query instead of writing to tmp.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' "$1" | while IFS= read -r line; do
    # Split the line into an array
    IFS=$'\t' read -r -a arr <<< "$line"
    
    # Concatenate the first four columns with '_'
    pos_info="${arr[0]}_${arr[1]}_${arr[2]}_${arr[3]}"
    
    # Initialize strings for the current line's total depth and variant support
    total_depth_line="$pos_info"
    variant_support_line="$pos_info"
    
    # Loop through the allele depth (AD) fields (starting from the 5th element)
    for (( i=4; i<${#arr[@]}; i++ )); do
        # Replace '.' with '0' in the AD field
        ad_field="${arr[i]//./0}"
        
        # Split the AD field into reference and variant reads
        IFS=',' read -r ref_reads var_reads <<< "$ad_field"
        
        # Ensure ref_reads and var_reads are integers
        if ! [[ $ref_reads =~ ^[0-9]+$ ]] || ! [[ $var_reads =~ ^[0-9]+$ ]]; then
            echo "Error: Non-integer values found in AD field" >&2
            exit 2
        fi
        
        # Calculate total depth and append to the line string
        total_depth=$((ref_reads + var_reads))
        total_depth_line="$total_depth_line\t$total_depth"
        
        # Append variant reads to the line string
        variant_support_line="$variant_support_line\t$var_reads"
    done
    
    # Append the processed lines to their respective files
    echo -e "$total_depth_line" >> "$total_depth_matrix"
    echo -e "$variant_support_line" >> "$variant_support_matrix"

done
