import os
import subprocess
from utils import log, run_cmd

class Pileup:
    def __init__(self, outgroup, aligners, base_output_dir, run_id = None, no_cache=False, verbose=True):
        self.reference = outgroup.name
        self.output_dir = base_output_dir
        self.outgroup = outgroup
        self.no_cache = no_cache
        self.verbose = verbose
        self.ref_fasta = outgroup.fasta_path
        self.bams = aligners
        self.taxon_names = [aligner.species for aligner in self.bams]
        self.run_id = run_id if run_id else f"{self.reference}__{'__'.join(self.taxon_names)}"

        self.pileup_path = f"{self.output_dir}/{self.run_id}.pileup.gz"

    def _check_file(self, path):
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Missing file: {path}")

    def generate(self):
        # Check files
        for path in [self.ref_fasta] + [aligner.final_bam for aligner in self.bams]:
            self._check_file(path)

        # Skip if exists and caching allowed
        if os.path.exists(self.pileup_path) and not self.no_cache:
            log(f"Pileup already exists: {self.pileup_path}", self.verbose)
            return self.pileup_path

        log(f"Generating pileup: {self.pileup_path}", self.verbose)
        cmd = ["samtools", "mpileup", "-f", self.ref_fasta, "-B", "-d", "100"] + \
            [bam.final_bam for bam in self.bams]

        with open(self.pileup_path, "wb") as out:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            gzip_proc = subprocess.Popen(["gzip"], stdin=proc.stdout, stdout=out)
            proc.stdout.close()
            gzip_proc.communicate()

        log(f"Pileup written to: {self.pileup_path}", self.verbose)
        return self.pileup_path
