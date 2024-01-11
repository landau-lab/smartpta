#!/bin/bash

# Check if the correct number of input arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_file> <output_file>"
    exit 1
fi

input_file="$1"
output_file="$2"

# Initialize the merged table using the IDs from the first file
read -r first_file < "$input_file"
awk '{print $1}' "$first_file" > "$output_file"

# Create a temporary file for the header
header_temp="header.txt"
echo -n "ID" > "$header_temp"

# Append the counts from each file and build the header simultaneously
while read -r file; do
    # Append the counts
    paste "$output_file" <(awk '{print $2}' "$file") > temp && mv temp "$output_file"
    # Add the trimmed filename to the header, removing ".counts"
    filename=$(basename "$file" .counts)
    echo -n -e "\t$filename" >> "$header_temp"
done < "$input_file"

# Add a newline at the end of the header
echo "" >> "$header_temp"

# Add the header to the merged table
cat "$header_temp" "$output_file" > temp && mv temp "$output_file"

# Clean up: Remove the temporary header file
rm "$header_temp"