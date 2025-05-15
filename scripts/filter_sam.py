#!/usr/bin/env python3
import sys
import argparse

# Constants
low_mapq_threshold = 5
default_mapq_threshold = 60
default_offset = 75

def parse_args():
    parser = argparse.ArgumentParser(description="Filter SAM file based on MAPQ and read continuity.")
    parser.add_argument("input_sam", help="Input SAM file or '-' for stdin")
    parser.add_argument("--output", "-o", help="Output SAM file (default: stdout)")
    parser.add_argument("--mapq", type=int, default=default_mapq_threshold,
                        help=f"MAPQ threshold (default: {default_mapq_threshold})")
    parser.add_argument("--offset", type=int, default=default_offset,
                        help=f"Expected distance between consecutive reads (default: {default_offset})")
    return parser.parse_args()

# Counters for reporting
total_reads = 0
kept_reads = 0
filtered_mapq = 0
filtered_disjoint = 0
filtered_chr = 0

def filter_read(prev_pos, cur_pos, next_pos, cur_mapq, cur_chr, cur_lines, outfile, mapq_threshold, offset):
    global total_reads, kept_reads, filtered_mapq, filtered_disjoint, filtered_chr
    for line, pos, mapq, chr in zip(cur_lines, cur_pos, cur_mapq, cur_chr):
        total_reads += 1
        if any(keyword in chr for keyword in ['Un', 'random', 'alt', 'fix', 'hap']):
            filtered_chr += 1
            continue
        if mapq < low_mapq_threshold:
            filtered_mapq += 1
            continue
        if mapq >= mapq_threshold or pos - offset in prev_pos + next_pos or pos + offset in prev_pos + next_pos:
            outfile.write(line)
            kept_reads += 1
        else:
            filtered_disjoint += 1

def main():
    args = parse_args()

    # Input
    infile = sys.stdin if args.input_sam == "-" else open(args.input_sam, 'r')
    outfile = open(args.output, 'w') if args.output else sys.stdout

    prev_read, cur_read, next_read = None, None, None
    prev_chr, cur_chr, next_chr = [], [], []
    prev_lines, cur_lines, next_lines = [], [], []
    prev_pos, cur_pos, next_pos = [], [], []
    prev_mapq, cur_mapq, next_mapq = [], [], []

    for i, line in enumerate(infile):
        if line.startswith('@'):
            outfile.write(line)
            continue

        fields = line.split('\t')
        try:
            read_name = fields[0]
            chr = fields[2]
            pos = int(fields[3])
            mapq = int(fields[4])
        except ValueError:
            print(f"⚠️ Invalid line {i+1}: POS={fields[3]} MAPQ={fields[4]}", file=sys.stderr)
            continue

        if next_read == read_name:
            next_lines.append(line)
            next_pos.append(pos)
            next_mapq.append(mapq)
            next_chr.append(chr)
        else:
            filter_read(prev_pos, cur_pos, next_pos, cur_mapq, cur_chr, cur_lines, outfile, args.mapq, args.offset)
            prev_lines, cur_lines, next_lines = cur_lines, next_lines, [line]
            prev_pos, cur_pos, next_pos = cur_pos, next_pos, [pos]
            prev_mapq, cur_mapq, next_mapq = cur_mapq, next_mapq, [mapq]
            prev_chr, cur_chr, next_chr = cur_chr, next_chr, [chr]
            prev_read, cur_read, next_read = cur_read, next_read, read_name

    # Final batch
    filter_read(prev_pos, cur_pos, next_pos, cur_mapq, cur_chr, cur_lines, outfile, args.mapq, args.offset)

    # Cleanup
    if args.output:
        outfile.close()
    if args.input_sam != "-":
        infile.close()

    # Summary report
    print(f"🧾 Filter Summary:", file=sys.stderr)
    print(f"  Total reads processed:   {total_reads}", file=sys.stderr)
    print(f"  Reads kept:              {kept_reads}", file=sys.stderr)
    print(f"  Filtered (low MAPQ):     {filtered_mapq}", file=sys.stderr)
    print(f"  Filtered (disjoint):     {filtered_disjoint}", file=sys.stderr)
    print(f"  Filtered (alt contigs):  {filtered_chr}", file=sys.stderr)

if __name__ == "__main__":
    main()
