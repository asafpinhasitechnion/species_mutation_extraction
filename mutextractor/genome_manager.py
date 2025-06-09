import os
import subprocess
import shutil
from utils import run_cmd, log
from Bio import SeqIO



class Genome:
    def __init__(self, name, accession, output_dir, no_cache=False, verbose = True):
        self.name = name
        self.accession = accession
        self.output_dir = os.path.join(output_dir, name)
        self.no_cache = no_cache
        self.fasta_path = os.path.join(self.output_dir, f"{name}.fasta")
        self.fastq_path = None
        self.verbose = verbose

    def download(self):
        if os.path.exists(self.fasta_path) and not self.no_cache:
            log(f"Genome for {self.name} found in cache. Skipping download.", self.verbose)
            return

        os.makedirs(self.output_dir, exist_ok=True)
        temp_dir = os.path.join(os.path.dirname(self.output_dir), f"temp_{self.name}")
        os.makedirs(temp_dir, exist_ok=True)

        log(f"Downloading genome for {self.name} ({self.accession})", self.verbose)
        zip_path = os.path.join(temp_dir, f"{self.name}.zip")
        run_cmd([
            "datasets", "download", "genome", "accession", self.accession,
            "--filename", zip_path
        ], verbose=self.verbose)

        log("Unzipping genome...", self.verbose)
        run_cmd(["unzip", "-q", "-o", zip_path, "-d", temp_dir], verbose=self.verbose)


        fna_source = None
        for root, _, files in os.walk(temp_dir):
            for file in files:
                if file.endswith(".fna"):
                    fna_source = os.path.join(root, file)
                    break

        if not fna_source:
            log("Error: No FASTA (.fna) file found in the downloaded data.", self.verbose)
            shutil.rmtree(temp_dir)
            raise FileNotFoundError("Missing FASTA")

        # Move the FASTA file to final location (efficient alternative to copy)
        shutil.move(fna_source, self.fasta_path)

        # Strip header and clean up
        log("Stripping FASTA headers (keeping only the first word)...", self.verbose)
        subprocess.run(["sed", "-i", "s/ .*$//", self.fasta_path], check=True)

        shutil.rmtree(temp_dir)



    def index(self, aligner="bwa"):
        log(f"Indexing genome for {self.name} with {aligner}...", self.verbose)

        required_bwa_exts = ["amb", "ann", "bwt", "pac", "sa"]
        required_all = [f"{self.fasta_path}.{ext}" for ext in required_bwa_exts]
        required_all.append(f"{self.fasta_path}.fai")  # samtools index

        all_exist = all(os.path.exists(f) for f in required_all)

        if all_exist and not self.no_cache:
            log(f"Index files already exist for {self.name}. Skipping indexing.", self.verbose)
            return

        for ext in required_bwa_exts + ["fai"]:
            path = f"{self.fasta_path}.{ext}"
            if os.path.exists(path):
                os.remove(path)

        # if aligner == "bwa":
        run_cmd(["bwa", "index", self.fasta_path], verbose = self.verbose)
        # else:
        #     raise ValueError(f"Aligner '{aligner}' not supported for indexing.")

        run_cmd(["samtools", "faidx", self.fasta_path], verbose = self.verbose)

        log(f"Indexing complete for {self.name}", self.verbose)


    def generate_fragment_fastq(self, length=150, offset=75, force=False):
        output_fastq = os.path.join(self.output_dir, f"{self.name}.fastq")
        
        if os.path.exists(output_fastq) and not force:
            log(f"Fastq for {self.name} already exists. Skipping.", self.verbose)
            self.fastq_path = output_fastq
            return output_fastq


        def split_sequence(seq, offset, section_size):
            last_start = len(seq) - section_size
            sections = [(i, seq[i:i+section_size]) for i in range(0, last_start, offset)]
            if last_start > 0:
                sections.append((last_start, seq[last_start:]))
            return sections

        chromosomes = {
            record.id: record.seq
            for record in SeqIO.parse(self.fasta_path, "fasta")
        }

        with open(output_fastq, 'w') as out:
            for chrom_name, sequence in chromosomes.items():
                for i, (start, frag) in enumerate(split_sequence(sequence, offset, length)):
                    end = start + len(frag) - 1
                    out.write(f"@{chrom_name}_{start+1}_{end+1}\n")
                    out.write(f"{frag}\n+\n")
                    out.write(f"{'I' * len(frag)}\n")

        log(f"Wrote {output_fastq}", self.verbose)
        self.fastq_path = output_fastq
        return output_fastq
