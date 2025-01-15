# Process multiple GTF files and create a combined data structure.
# Returns a dictionary: chr_intron_chain -> {attributes from each file}
# Save attributes as a TSV file


from collections import defaultdict
import re
import sys

# Process a single GTF file and return two dictionaries:
# 1. transcript_attributes: transcript_id -> attributes
# 2. transcript_coordinates: transcript_id -> list of exon coordinates
def process_gtf_file(file_path):
    transcript_attributes = {}
    transcript_chr = {}
    transcript_exon_coordinates = defaultdict(list)
    
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            fields = line.strip().split('\t')
            if len(fields) != 9:
                print(f"Warning: Line does not have 9 columns: {line}", file=sys.stderr)
                continue
                
            feature_type = fields[2]
            if feature_type != 'exon':
                continue
                
            attributes = parse_gtf_attributes(fields[8])
            if 'transcript_id' not in attributes:
                print(f"Warning: lines does not have transcript_id: {line}", file=sys.stderr)
                continue
            # remove undesired attributes
            attributes.pop('exon', None)        
                
            transcript_id = attributes['transcript_id']
            start = int(fields[3])
            end = int(fields[4])
            transcript_exon_coordinates[transcript_id].append((start, end))
            
            # Store attributes
            if transcript_id not in transcript_attributes:
                transcript_attributes[transcript_id] = attributes
            if transcript_id not in transcript_chr:
                transcript_chr[transcript_id] = str(fields[0])
    
    # Sort exon coordinates by start position
    for id in transcript_exon_coordinates:
        transcript_exon_coordinates[id] = sorted(transcript_exon_coordinates[id])

    return transcript_attributes, transcript_chr, transcript_exon_coordinates


# Parse GTF attribute string into a dictionary.
def parse_gtf_attributes(attribute_string):
    attributes = {}
    # Split by semicolon and strip whitespace
    pairs = [x.strip() for x in attribute_string.split(';') if x.strip()]
    for pair in pairs:
        key, value = pair.split('"', 1)
        value = value.strip('"')
        key = key.strip()
        attributes[key] = value
    return attributes

# Generate intron chain string from sorted exon coordinates.
# Returns: list-of-tuples as intron coordinates
# Returns: if single exon, returns [(-1, -1), [exonStart, exonEnd]]
def get_intron_chain(exon_coordinates):
    introns = []
    sorted_exons = sorted(exon_coordinates, key=lambda x: x[0])
    
    if len(sorted_exons) == 0:
        return [(-1, -1)]
    
    if len(sorted_exons) == 1:
        return [(-1, -1), (sorted_exons[0][0], sorted_exons[0][1])]

    for i in range(len(sorted_exons) - 1):
        intron_start = sorted_exons[i][1] + 1
        intron_end = sorted_exons[i + 1][0] - 1
        assert intron_start <= intron_end
        introns.append((intron_start, intron_end))    
    return introns


# Process multiple GTF files and create a combined data structure.
# Returns a dictionary: intron_chain -> {attributes from each file}
def process_multiple_gtf_files(file_paths):
    all_file_data = []
    
    # Process each file
    for file_path in file_paths:
        attributes, chromosomes, coordinates = process_gtf_file(file_path)
        all_file_data.append((attributes, chromosomes, coordinates))
    
    # Create combined data structure
    combined_data = defaultdict(lambda: {f'attr{i}': {} for i in range(len(file_paths))})
    
    # Process each file's data
    for file_idx, (attributes, chromosomes, coordinates) in enumerate(all_file_data):
        # Group transcripts by their intron chain
        chain_to_transcripts = defaultdict(list)
        for transcript_id, coord_list in coordinates.items():
            intron_chain = get_intron_chain(coord_list)
            intron_chain_str = chromosomes[transcript_id] + ":" + ','.join(f'{start}-{end}' for start, end in intron_chain)
            chain_to_transcripts[intron_chain_str].append(transcript_id)
        
        # For each intron chain, combine attributes of all transcripts sharing it
        for intron_chain_str, transcript_ids in chain_to_transcripts.items():
            combined_attrs = {}
            for transcript_id in transcript_ids:
                combined_attrs.update(attributes[transcript_id])
            combined_data[intron_chain_str][f'attr{file_idx}'] = combined_attrs
    
    return combined_data

# Save the combined data to a TSV file
def save_results(combined_data, output_file):
    # First, collect all unique attribute keys from all files
    all_keys = set()
    for intron_chain, file_attributes in combined_data.items():
        for file_attrs in file_attributes.values():
            all_keys.update(file_attrs.keys())
    
    # Create header
    header = ['chr_intron_chain']
    for file_idx in range(len(next(iter(combined_data.values())))):
        for key in sorted(all_keys):
            header.append(f'attr{file_idx}_{key}')
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for intron_chain, file_attributes in combined_data.items():
            row = [intron_chain]
            for file_idx in range(len(file_attributes)):
                file_key = f'attr{file_idx}'
                attrs = file_attributes[file_key]
                for key in sorted(all_keys):
                    row.append(attrs.get(key, ''))
            f.write('\t'.join(str(x) for x in row) + '\n')

def main(listOfFiles):
    combined_data = process_multiple_gtf_files(listOfFiles)
    save_results(combined_data, 'combined_results.tsv')
    # print example
    first_chain = next(iter(combined_data))
    print(f"Total unique intron chains found: {len(combined_data)}")
    print("\nExample of the first intron chain and its attributes:")
    print(f"Intron chain: {first_chain}")
    print("Attributes from each file:")
    for file_num, attrs in combined_data[first_chain].items():
        print(f"{file_num}:", attrs)

def parse(argv):
    parser = argparse.ArgumentParser(description='Process GTF files and compare intron chains.')
    parser.add_argument('-a1', metavar='FILE', help='Allele1 GTF file (ground truth)')
    parser.add_argument('-a2', metavar='FILE', help='Allele2 GTF file (ground truth)')
    parser.add_argument('otherFiles', nargs='*', help='Additional GTF files')
    args = parser.parse_args(argv)

    listOfFiles = []
    if args.a1:
        listOfFiles.append(args.a1)
    if args.a2:
        listOfFiles.append(args.a2)
    listOfFiles.extend(args.otherFiles)

    if not listOfFiles:
        parser.error("At least one GTF file must be provided")

    for f in listOfFiles:
        if not f.endswith(".gtf"):
            raise ValueError(f"File {f} is not a GTF file!")
    
    args.files = listOfFiles

    return args

if __name__ == "__main__":
    args = parse(sys.argv[1:])
    main(args.files)