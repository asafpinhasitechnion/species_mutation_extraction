from io import TextIOWrapper
import os
import subprocess
import sys
import multiprocessing
from utils import run_cmd, log  
from typing import Optional, TextIO
import matplotlib.pyplot as plt


# Constants
LOW_MAPQ_THRESHOLD = 1
DEFAULT_MAPQ_THRESHOLD = 60
DEFAULT_OFFSET = 75


def log_to_file(log_path: Optional[str], message: str):
    if log_path:
        with open(log_path, 'a') as f:
            f.write(message + '\n')

def filter_sam(
    input_stream: TextIO,
    output_stream: TextIO,
    low_mapq: int = 10,
    mapq_threshold: int = 60,
    offset: int = 75,
    mapq_hist_folder: Optional[str] = None,
    hist_name: str = None,
    continuity: bool = True,
    verbose: bool = True,
    log_path: Optional[str] = None,
):
    total_reads = 0
    kept_reads = 0
    filtered_mapq = 0
    filtered_disjoint = 0
    filtered_chr = 0
    mapq_values = []

    def write_summary():
        lines = [
            "Filter Summary:",
            f"  Total reads processed:   {total_reads}",
            f"  Reads kept:              {kept_reads}",
            f"  Filtered (low MAPQ):     {filtered_mapq}",
            f"  Filtered (disjoint):     {filtered_disjoint}",
            f"  Filtered (alt contigs):  {filtered_chr}"
        ]
        if verbose:
            for line in lines:
                log(line, verbose)
        for line in lines:
            log_to_file(log_path, line)

    def plot_histogram():
        if mapq_hist_folder:
            plt.figure(figsize=(8, 5))
            plt.hist(mapq_values, bins=50, color='steelblue', edgecolor='black', log=True)
            plt.title("MAPQ Score Distribution")
            plt.xlabel("MAPQ")
            plt.ylabel("Read Count (log scale)")
            plt.tight_layout()
            os.makedirs(mapq_hist_folder, exist_ok=True)
            out_path = os.path.join(mapq_hist_folder, hist_name)
            plt.savefig(out_path)
            msg = f"MAPQ histogram saved to {out_path}"
            log(msg, verbose)
            log_to_file(log_path, msg)

    def filter_read(prev_pos, cur_pos, next_pos, cur_mapq, cur_chr, cur_lines):
        nonlocal total_reads, kept_reads, filtered_mapq, filtered_disjoint, filtered_chr
        for line, pos, mapq, chr in zip(cur_lines, cur_pos, cur_mapq, cur_chr):
            total_reads += 1
            if any(keyword in chr for keyword in ['Un', 'random', 'alt', 'fix', 'hap']):
                filtered_chr += 1
                continue
            if continuity:
                if mapq < low_mapq:
                    filtered_mapq += 1
                    continue
                if (mapq >= mapq_threshold or 
                    (pos - offset in prev_pos + next_pos) or 
                    (pos + offset in prev_pos + next_pos)):
                    output_stream.write(line)
                    kept_reads += 1
                else:
                    filtered_disjoint += 1
            else:
                if mapq < mapq_threshold:
                    filtered_mapq += 1
                else:
                    output_stream.write(line)
                    kept_reads += 1

    cur_read, next_read = None, None
    cur_chr, next_chr = [], []
    cur_lines, next_lines = [], []
    prev_pos, cur_pos, next_pos = [], [], []
    cur_mapq, next_mapq = [], []

    for i, line in enumerate(input_stream):
        if line.startswith('@'):
            output_stream.write(line)
            continue
        fields = line.split('\t')
        try:
            read_name = fields[0]
            chr = fields[2]
            pos = int(fields[3])
            mapq = int(fields[4])
            mapq_values.append(mapq)
        except Exception:
            msg = f"Invalid line {i+1}"
            print(msg, file=sys.stderr)
            log_to_file(log_path, msg)
            continue
        if next_read == read_name:
            next_lines.append(line)
            next_pos.append(pos)
            next_mapq.append(mapq)
            next_chr.append(chr)
        else:
            filter_read(prev_pos, cur_pos, next_pos, cur_mapq, cur_chr, cur_lines)
            prev_pos, cur_pos, next_pos = cur_pos, next_pos, [pos]
            cur_mapq, next_mapq = next_mapq, [mapq]
            cur_chr, next_chr = next_chr, [chr]
            cur_lines, next_lines = next_lines, [line]
            cur_read, next_read = next_read, read_name

    filter_read(prev_pos, cur_pos, next_pos, cur_mapq, cur_chr, cur_lines)
    write_summary()
    plot_histogram()



class Aligner:
    def __init__(
        self,
        species_genome,
        reference_genome,
        base_output_dir,
        aligner_cmd=None,
        aligner_name=None,
        no_cache=False,
        cores = None,
        verbose=True
    ):
        self.species = species_genome.name
        self.reference = reference_genome.name
        self.output_dir = base_output_dir
        self.no_cache = no_cache
        self.verbose = verbose
        self.cores = cores if cores else multiprocessing.cpu_count()

        self.species_fasta = species_genome.fasta_path
        self.reference_fasta = reference_genome.fasta_path
        self.fastq = species_genome.fastq_path
        self.bam_dir = f"{self.output_dir}/BAMs"
        self.plots_dir = f"{self.output_dir}/Plots"
        self.raw_bam = f"{self.bam_dir}/{self.species}_to_{self.reference}_raw.bam"
        self.final_bam = f"{self.bam_dir}/{self.species}_to_{self.reference}.bam"
        self.hist_name = f"{self.species}_to_{self.reference}.png"
        self.log_path = self.final_bam.replace(".bam", ".log")

        os.makedirs(self.bam_dir, exist_ok=True)
        os.makedirs(self.plots_dir, exist_ok=True)

        if aligner_cmd:
            self.aligner_cmd_template = aligner_cmd
        elif aligner_name:
            self.aligner_cmd_template = self.get_aligner_cmd_from_name(aligner_name)
        else:
            log("No aligner specified. Defaulting to 'bwa'", self.verbose)
            self.aligner_cmd_template = self.get_aligner_cmd_from_name("bwa")

        self.validate_aligner_cmd()

    def get_aligner_cmd_from_name(self, name):
        commands = {
            "bwa": "bwa mem -t {cores} {ref} {fq}",
            "minimap2": "minimap2 -t {cores} -ax sr {ref} {fq}",
            "bbmap": "bbmap.sh ref={ref} threads={cores} in={fq} out=stdout.sam"
        }
        if name not in commands:
            raise ValueError(f"Unsupported aligner: {name}")
        return commands[name]

    def validate_aligner_cmd(self):
        required = ["{ref}", "{fq}", "{cores}"]
        for r in required:
            if r not in self.aligner_cmd_template:
                raise ValueError(f"--aligner-cmd must include placeholders: {', '.join(required)}")


    def align_streamed(self, mapq=60, offset=75, continuity=True, max_sort_mem=None):
        """
        Full streaming: aligner | filter | samtools sort -> BAM (low RAM, minimal I/O)
        """
        if os.path.exists(self.final_bam) and os.path.exists(self.final_bam + '.bai') and not self.no_cache:
            log(f"Streamed alignment already exists: {self.final_bam}", self.verbose)
            return self.final_bam

        cmd = self.aligner_cmd_template.replace("{ref}", self.reference_fasta) \
                                    .replace("{fq}", self.fastq) \
                                    .replace("{cores}", str(self.cores))

        log("Running full streaming alignment + filtering + sorting...", self.verbose)

        align_proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

        sort_cmd = ["samtools", "sort", "-@", str(self.cores), "-o", self.final_bam]
        if max_sort_mem:
            sort_cmd += ["-m", str(max_sort_mem)]

        sort_proc = subprocess.Popen(sort_cmd, stdin=subprocess.PIPE)

        assert align_proc.stdout is not None and sort_proc.stdin is not None

        with TextIOWrapper(align_proc.stdout) as reader, TextIOWrapper(sort_proc.stdin, write_through=True, buffering=1) as writer:
            filter_sam(
                input_stream=reader,
                output_stream=writer,
                mapq_threshold=mapq,
                offset=offset,
                continuity=continuity,
                mapq_hist_folder=self.plots_dir,
                hist_name=self.hist_name,
                verbose=self.verbose,
                log_path=str(self.log_path),
            )

        sort_proc.stdin.close()
        sort_proc.wait()

        run_cmd(["samtools", "index", self.final_bam])
        log(f"Finished (streamed): {self.final_bam}", self.verbose)
        return self.final_bam
    

    def align_disk_cached(self, mapq=60, offset=75, continuity=True):
        """
        Disk-cached: aligner -> raw.bam -> filter -> final.bam (safer on RAM, uses disk)
        """
        raw_bam = self.raw_bam
        final_bam = self.final_bam

        # Step 1: Align to raw.bam if missing or no_cache
        if not os.path.exists(raw_bam) or self.no_cache:
            cmd = self.aligner_cmd_template.replace("{ref}", self.reference_fasta) \
                                        .replace("{fq}", self.fastq) \
                                        .replace("{cores}", str(self.cores))
            log("Running alignment to raw BAM...", self.verbose)
            run_cmd(f"{cmd} | samtools view -bS - > {raw_bam}", shell=True)

        # Step 2: Filter to final.bam if missing or no_cache
        if not os.path.exists(final_bam) or self.no_cache:
            log("Filtering from raw BAM to final BAM...", self.verbose)
            view_proc = subprocess.Popen(["samtools", "view", "-h", raw_bam], stdout=subprocess.PIPE)
            sort_proc = subprocess.Popen(
                ["samtools", "sort", "-@", str(self.cores), "-o", final_bam],
                stdin=subprocess.PIPE
            )
            assert view_proc.stdout is not None and sort_proc.stdin is not None

            filter_sam(
                input_stream=TextIOWrapper(view_proc.stdout),
                output_stream=TextIOWrapper(sort_proc.stdin),
                mapq_threshold=mapq,
                offset=offset,
                continuity=continuity,
                mapq_hist_folder=self.plots_dir,
                hist_name=self.hist_name,
                verbose=self.verbose,
                log_path=str(self.log_path)
            )

            sort_proc.stdin.close()
            sort_proc.wait()

        # Step 3: Index final.bam if .bai missing or no_cache
        if not os.path.exists(final_bam + '.bai') or self.no_cache:
            run_cmd(["samtools", "index", final_bam])

        log(f"Finished (disk-cached): {final_bam}", self.verbose)
        return final_bam

