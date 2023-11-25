#!/bin/bash

for file in *_R1.fq.gz; do
    sample_name=$(echo "$file" | cut -d '_' -f 1)
    mkdir -p "$sample_name"
    mv "$sample_name"_R1.fq.gz "$sample_name"/
    mv "$sample_name"_R2.fq.gz "$sample_name"/
done
