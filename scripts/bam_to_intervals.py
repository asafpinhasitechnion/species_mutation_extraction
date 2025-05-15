#!/usr/bin/env python3
import argparse
import pysam
import os


def parse_args():
    parser = argparse.ArgumentParser(description="Extract genomic intervals from SAM/BAM reads.")
    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("output_path", help="Path to write output intervals")
    parser.add_argument("--sorted", action="store_true", help="Indicate if file is coordinate-sorted to merge on the fly")
    parser.add_argument("--merge", action="store_true", help="Disable merging of intervals")
    return parser.parse_args()


def write_intervals(intervals, out_path):
    with open(out_path, 'w') as out:
        for chrom, start, end in intervals:
            out.write(f"{chrom}\t{start}\t{end}\n")


def merge_intervals(intervals):
    merged = []
    for chrom, start, end in sorted(intervals):
        if merged and merged[-1][0] == chrom and merged[-1][2] >= start:
            merged[-1] = (merged[-1][0], merged[-1][1], max(merged[-1][2], end))
        else:
            merged.append((chrom, start, end))
    return merged


def extract_intervals_sorted(bamfile):
    merged = []
    for read in bamfile.fetch():
        if read.is_unmapped:
            continue
        chrom = bamfile.get_reference_name(read.reference_id)
        start = read.reference_start
        end = read.reference_end
        if merged and merged[-1][0] == chrom and merged[-1][2] >= start:
            merged[-1] = (merged[-1][0], merged[-1][1], max(merged[-1][2], end))
        else:
            merged.append((chrom, start, end))
    return merged


def extract_raw_intervals(bamfile):
    intervals = []
    for read in bamfile.fetch():
        if read.is_unmapped:
            continue
        chrom = bamfile.get_reference_name(read.reference_id)
        start = read.reference_start
        end = read.reference_end
        intervals.append((chrom, start, end))
    return intervals


def main():
    args = parse_args()

    out_path = args.output_path
    os.makedirs(out_path, exist_ok=True)
    input_bam_name = os.path.splitext(os.path.basename(args.input_bam))[0]
    output_file = os.path.join(out_path, f"{input_bam_name}_intervals.tsv")

    bamfile = pysam.AlignmentFile(args.input_bam, "rb")

    if args.sorted:
        intervals = extract_intervals_sorted(bamfile)
    else:
        intervals = extract_raw_intervals(bamfile)
        if args.merge:
            intervals = merge_intervals(intervals)

    write_intervals(intervals, output_file)
    print(f"✅ Intervals written to: {output_file}")


if __name__ == "__main__":
    main()
