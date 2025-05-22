#!/usr/bin/env python3
import os
import sys
import argparse
import gzip
import csv

# === Constants ===
REMOVE_CHARS = str.maketrans('', '', '^$[]')
CHR_IDX, POS_IDX, REF_NUC_IDX = 0, 1, 2
PREV_IDX, CUR_IDX, NEXT_IDX = 0, 1, 2

# === Argument Parser ===
def parse_args():
    parser = argparse.ArgumentParser(description="Extract positions with conserved flanks and differing middle bases across species.")
    parser.add_argument("--n-species", type=int, required=True, help="Number of species in the pileup file")
    parser.add_argument("--pileup-file", required=True, help="Path to .pileup.gz file")
    parser.add_argument("--output", required=True, help="Path to output .csv.gz file")

    # Optional hardcoded test mode (for debugging)
    if len(sys.argv) == 1:
        test_args = [
            
            "--pileup-file", "../Output/test_run_mutiple_species/Leptophobia_aripa__Pieris_brassicae__Pieris_mannii__Pieris_napi__Pieris_rapae.pileup.gz",
            "--n-species", "5",
            "--output", "../Output/test_run_mutiple_species/Leptophobia_aripa__Pieris_brassicae__Pieris_mannii__Pieris_napi__Pieris_rapae.csv.gz"
        ]
        return parser.parse_args(test_args)
    return parser.parse_args()

# === Utilities ===
def all_same(seq):
    return len(seq) > 0 and all(ch == seq[0] for ch in seq)

def quality_check(fields):
    sample_fields = fields[3:]
    return (
        sample_fields
        and all('*' not in field for field in sample_fields)
        and all(all_same(field.translate(REMOVE_CHARS)) for field in sample_fields)
    )

def parse_line(line, n_species):
    parts = line.strip().split('\t')
    if len(parts) < n_species * 3:
        return None

    # Extract chromosome, position, reference base
    chrom, pos, ref_base = parts[:3]

    # Extract base strings: column 4 + 3*i for each sample
    base_calls = parts[4::3]

    # Normalize each base: take the first character or use ref if '.' or ','
    normalized = [
        base[0] if base and base[0] not in {',', '.'} else ref_base
        for base in base_calls
    ]
    return [chrom, pos, ref_base] + normalized

def detect_mutations(triplet_buffer):
    # Extract triplets from 3-line buffer (prev, curr, next)
    triplets = [fields[3:] for fields in triplet_buffer]

    prev_bases = triplets[PREV_IDX]
    curr_bases = triplets[CUR_IDX]
    next_bases = triplets[NEXT_IDX]

    if all_same(prev_bases) and all_same(next_bases) and len(set(curr_bases)) > 1:
        chrom = triplet_buffer[CUR_IDX][CHR_IDX]
        pos = triplet_buffer[CUR_IDX][POS_IDX]
        ref_base = triplet_buffer[CUR_IDX][REF_NUC_IDX]
        return [chrom, pos, prev_bases[0].upper(), next_bases[0].upper(), ref_base.upper()] + [base.upper() for base in curr_bases]
    return None

# === Main Function ===
def main():
    args = parse_args()
    n_species = args.n_species
    pileup_file = args.pileup_file
    output_path = args.output

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with gzip.open(pileup_file, 'rt') as infile, gzip.open(output_path, 'wt', newline='') as outfile:
        writer = csv.writer(outfile)
        header = ["chromosome", "position", "left", "right", "reference_base"] + [f"taxa{i}" for i in range(1, n_species)]
        writer.writerow(header)

        # Initialize 3-line buffer for sliding window
        buffer = [None, parse_line(infile.readline(), n_species), parse_line(infile.readline(), n_species)]
        qc_flags = [False, quality_check(buffer[1]), quality_check(buffer[2])]

        for line in infile:
            buffer = [buffer[1], buffer[2], parse_line(line, n_species)]
            qc_flags = [qc_flags[1], qc_flags[2], quality_check(buffer[2])]

            if all(qc_flags):
                result = detect_mutations(buffer)
                if result:
                    writer.writerow(result)

    print(f"✅ Matching triplet mutation positions written to: {output_path}")

if __name__ == "__main__":
    main()
