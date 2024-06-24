import argparse
import os
from Bio import SeqIO

def extract_qualifier(genbank_file, qualifier_key):
    features = {}
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                location = (int(feature.location.start), int(feature.location.end))
                qualifiers = feature.qualifiers
                qualifier_value = qualifiers.get(qualifier_key, [""])[0]
                features[location] = qualifier_value
    return features

def is_within_range(location1, location2, error_margin):
    return abs(location1[0] - location2[0]) <= error_margin and abs(location1[1] - location2[1]) <= error_margin

def create_mapping(old_features, new_features, output_file, error_margin):
    with open(output_file, 'w') as out_f:
        for old_location, old_annotation in old_features.items():
            matched = False
            for new_location, new_annotation in new_features.items():
                if is_within_range(old_location, new_location, error_margin):
                    out_f.write(f"{old_annotation}\t{new_annotation}\n")
                    matched = True
                    break
            if not matched:
                out_f.write(f"{old_annotation}\t\n")

def main():
    parser = argparse.ArgumentParser(description="Link old annotation tags with new annotation tags in GenBank files, allowing for a specified error margin in location matching.")
    parser.add_argument('--gbk1', required=True, help="Path to the old GenBank file.")
    parser.add_argument('--gbk2', required=True, help="Path to the new GenBank file.")
    parser.add_argument('--output', required=True, help="Path to the output file prefix for mappings.")
    parser.add_argument('--error_margin', type=int, default=20, help="Allowable error margin in bp for location matching (default: 20).")
    parser.add_argument('--gbk1_tag', required=True, help="Qualifier tag to extract from the old GenBank file.")
    parser.add_argument('--gbk2_tag', required=True, help="Qualifier tag to extract from the new GenBank file.")

    args = parser.parse_args()

    old_features = extract_qualifier(args.gbk1, args.gbk1_tag)
    new_features = extract_qualifier(args.gbk2, args.gbk2_tag)

    output_file1 = args.output + "_1.txt"
    output_file2 = args.output + "_2.txt"

    # Mapping from gbk1 to gbk2
    create_mapping(old_features, new_features, output_file1, args.error_margin)

    # Mapping from gbk2 to gbk1
    create_mapping(new_features, old_features, output_file2, args.error_margin)

if __name__ == "__main__":
    main()
