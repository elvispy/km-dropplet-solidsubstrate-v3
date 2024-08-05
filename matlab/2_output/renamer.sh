#!/bin/bash

# Traverse all subdirectories and find .mat files
find . -type f -name "*.mat" | while read file; do
    # Get the directory and filename
    dir=$(dirname "$file")
    base=$(basename "$file")
    
    # Replace ":" with "-" in the filename
    newbase=$(echo "$base" | tr ':' '-')
    
    # Rename the file
    mv "$file" "$dir/$newbase"
done

echo "Renaming completed."
