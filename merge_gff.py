import argparse

def parse_gff(file_path):
    regions = []
    with open(file_path, 'r') as file:
        for line in file:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            start = int(parts[3])
            end = int(parts[4])
            seq_tag = parts[1]
            regions.append((start, end, seq_tag, line.strip()))
    return regions

def merge_regions(regions, max_gap):
    merged_regions = []
    prev_start, prev_end, prev_seq_tag, prev_line = regions[0]
    merged_seq_tags = [prev_seq_tag]

    for start, end, seq_tag, line in regions[1:]:
        if start - prev_end <= max_gap:
            # Merge regions
            prev_end = max(prev_end, end)
            merged_seq_tags.append(seq_tag)
        else:
            # Save previous region
            merged_seq_tag = ','.join(merged_seq_tags)
            merged_regions.append((prev_start, prev_end, merged_seq_tag, prev_line))
            prev_start, prev_end, prev_seq_tag, prev_line = start, end, seq_tag, line
            merged_seq_tags = [seq_tag]
    
    # Add the last region
    merged_seq_tag = ','.join(merged_seq_tags)
    merged_regions.append((prev_start, prev_end, merged_seq_tag, prev_line))
    
    return merged_regions

def filter_regions(merged_regions, min_length):
    if min_length is None:
        return merged_regions
    else:
        filtered_regions = []
        for start, end, merged_seq_tag, line in merged_regions:
            if end - start >= min_length:
                filtered_regions.append((start, end, merged_seq_tag, line))
        return filtered_regions

def update_gff_lines(filtered_regions):
    updated_lines = []
    for start, end, merged_seq_tag, line in filtered_regions:
        parts = line.split('\t')
        parts[1] = merged_seq_tag
        parts[3] = str(start)
        parts[4] = str(end)
        updated_line = '\t'.join(parts)
        updated_lines.append(updated_line)
    return updated_lines

def save_gff(file_path, updated_lines):
    with open(file_path, 'w') as file:
        for line in updated_lines:
            file.write(line + '\n')

def main():
    parser = argparse.ArgumentParser(description='Merge GFF regions less than a specified gap and filter by length.')
    parser.add_argument('--input_file', type=str, required=True, help='Input GFF file')
    parser.add_argument('--output_file', type=str, required=True, help='Output GFF file')
    parser.add_argument('--max_gap', type=int, required=True, help='Maximum gap between regions to merge')
    parser.add_argument('--drop_len_below', type=int, help='Drop merged regions with length below this value')
    
    args = parser.parse_args()
    
    # Processing
    regions = parse_gff(args.input_file)
    merged_regions = merge_regions(regions, args.max_gap)
    filtered_regions = filter_regions(merged_regions, args.drop_len_below)
    updated_lines = update_gff_lines(filtered_regions)
    save_gff(args.output_file, updated_lines)
    
    print(f"Regions merged and saved to {args.output_file}")

if __name__ == '__main__':
    main()
