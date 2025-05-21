
# 🧬 Multi-Species Mutation Extraction Pipeline

This repository contains a tool for scalable extraction, detection and analysis of mutations in species' evolutionary history.
The pipeline aligns multiple species to a reference genome, simulates reads, filters by mapping quality, extracts unambiguous triplet mutations, and visualizes mutation and coverage patterns.

---

## 🔧 Core Concepts

### 🧪 Mutation Context Extraction
The pipeline identifies positions in the genome where:
- All species share the **same flanking bases**
- One species has a **distinct middle base**

This allows the extraction of **unambiguous point mutations** in triplet format, enabling detailed mutation spectrum analysis. 
The distribution of mutation spectra, mutation rates and mutation types across the reference genome can be visualized and processed.

### 🔁 Reference-Based Alignment Strategy
Each non-outgroup species is aligned independently to a shared **reference genome (outgroup)**. This avoids full multiple alignment and makes the process **scalable**. 
The tool is increadibly lightweight compared to stadard multiple species aligner, making it **increadibly fit for parallel running** across different species groups.

---

## 🚀 Functional Workflow

### Step 1: Genome Preparation
- Download genomes by NCBI accession
- Index reference genome for alignment tools

### Step 2: Read Simulation and Alignment
- Simulate FASTQ reads by sliding a window across the genome
- Align simulated reads to the outgroup
- Filter alignments by MAPQ and coverage

### Step 3: Mutation Detection
- Generate a pileup of reference + aligned BAMs
- Extract unambiguous triplet mutation events
- Optionally: extract the full mutations list, including genomic position

### Step 4: Normalization & Analysis
- Normalize mutation counts by underlying triplet abundance
- Collapse complementary strands for canonical spectra
- Visualize results: mutation spectra, genomic distribution, coverage

---

## ⚙️ Example Usage

### 🧰 Run the full pipeline
```bash
python run_pipeline.py Pieris_rapae GCF_905147795.1 \
                       Pieris_napi GCF_905475465.1 \
                       Leptophobia_aripa GCA_951799465.1 \
                       --mapq 60 --aligner bwa --remove-temp --genomic-position-plots
```

### 🛠️ Use an individual tool
```bash
python extract_mutations.py Leptophobia_aripa Pieris_rapae Pieris_napi \
  --pileup-dir Output --output-dir Output
```

---

## 💡 Caching & Cleanup

All stages support optional caching:
- Use `--no-cache` to **force regeneration** of intermediate files
- Default behavior: **skip re-running** steps if outputs exist

To manage disk usage, pass:
```bash
--remove-temp
```
This deletes intermediate BAMs, pileups, FASTQs, and intervals once they're no longer needed.

A dedicated script `cleanup_output.sh` allows **manual cleanup** by component:
```bash
bash cleanup_output.sh Output/ --bams --pileup --intervals --genomes
```

---

## ⚡ Scalability & Parallelization

- Each **species vs reference** alignment is independent → can be parallelized
- Output's for each run are saved in separate folders, allowing parallelization without conflict
- No full multiple sequence alignment required
- Easily scaled to hundreds of genomes with a job scheduler (e.g. SLURM/Condor)

This makes the pipeline ideal for large-scale evolutionary or mutational studies across clades.

---

## 📦 Output Summary

- `*.pileup.gz` – multi-taxa pileup for reference + aligned BAMs
- `*_mutations.csv.gz` – full triplet mutation calls
- `*_mutations.json` – mutation context counts
- `Tables/*.tsv` – normalized, collapsed spectra
- `Plots/*.png` – spectra, genome-wide mutation distributions, MAPQ histograms

---

## 📊 Optional Visualizations

- **Mutation spectra** (raw & normalized)
- **Triplet distributions** per taxon
- **MAPQ score histogram** (for filtering diagnostics)
- **Mutation & coverage distribution** across chromosomes

---

## 🧪 Requirements

- Python ≥ 3.7
- Tools: `bwa`, `samtools`, `datasets`, and optionaly `minimap2`, `bbmap` or other aligners
- Python libs: `pandas`, `numpy`, `matplotlib`

Install with:
```bash
conda create -n mutextractor_env python=3.10 pandas numpy matplotlib
conda activate mutextractor_env
```
