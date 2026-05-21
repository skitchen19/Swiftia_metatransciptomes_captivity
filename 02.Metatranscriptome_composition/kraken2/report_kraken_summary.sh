#!/bin/bash

# Output file
echo -e "Sample\tBacteria\tEukaryota\tViruses\tArchaea\tUnclassified" > kraken_domain_summary.tsv

for REPORT in *.report; do
    SAMPLE=$(basename "$REPORT" .report)

    # Strip leading/trailing whitespace, match exact taxon name at end of line
    BACT=$(awk '$6 ~ /Bacteria$/ {print $1}' "$REPORT" | head -n 1)
    EUK=$(awk '$6 ~ /Eukaryota$/ {print $1}' "$REPORT" | head -n 1)
    VIR=$(awk '$6 ~ /Viruses$/ {print $1}' "$REPORT" | head -n 1)
    ARCH=$(awk '$6 ~ /Archaea$/ {print $1}' "$REPORT" | head -n 1)
    UNCLASS=$(awk '$6 ~ /unclassified$/ {print $1}' "$REPORT" | head -n 1)

    # Replace empty fields with 0
    BACT=${BACT:-0}
    EUK=${EUK:-0}
    VIR=${VIR:-0}
    ARCH=${ARCH:-0}
    UNCLASS=${UNCLASS:-0}

    echo -e "${SAMPLE}\t${BACT}\t${EUK}\t${VIR}\t${ARCH}\t${UNCLASS}" >> kraken_domain_summary.tsv
done

