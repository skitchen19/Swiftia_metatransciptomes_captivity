#!bin/bash

for f in *.*; do
    filename="${f%.*}"
    extension="fasta"
    new_filename=$(echo "$filename" | cut -d'.' -f1)
    mv "$f" "${new_filename}.${extension}"
done
