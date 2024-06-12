import argparse
import os
from Bio import SeqIO

def extract_locus_tags(genbank_file):
    features = {}
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                location = (int(feature.location.start), int(feature.location.end))
                qualifiers = feature.qualifiers
                locus_tag = qualifiers.get("locus_tag", [""])[0]
                features[location] = locus_tag
    return features

def extract_db_xrefs(genbank_file):
    features = {}
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                location = (int(feature.location.start), int(feature.location.end))
                qualifiers = feature.qualifiers
                db_xref = qualifiers.get("db_xref", [""])[0]
                features[location] = db_xref
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
    parser = argparse.ArgumentParser(description="Link old annotation locus tags with new annotation db_xrefs in GenBank files, allowing for a specified error margin in location matching.")
    parser.add_argument('--gbk1', required=True, help="Path to the old GenBank file with locus tags.")
    parser.add_argument('--gbk2', required=True, help="Path to the new GenBank file with db_xrefs.")
    parser.add_argument('--output', required=True, help="Path to the output file prefix for mappings.")
    parser.add_argument('--error_margin', type=int, default=20, help="Allowable error margin in bp for location matching (default: 20).")

    args = parser.parse_args()

    old_features = extract_locus_tags(args.gbk1)
    new_features = extract_db_xrefs(args.gbk2)

    output_file1 = args.output + "_1.txt"
    output_file2 = args.output + "_2.txt"

    # Mapping from gbk1 to gbk2
    create_mapping(old_features, new_features, output_file1, args.error_margin)

    # Mapping from gbk2 to gbk1
    create_mapping(new_features, old_features, output_file2, args.error_margin)

if __name__ == "__main__":
    main()
