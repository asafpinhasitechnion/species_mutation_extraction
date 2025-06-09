# species_mutation_extraction/mutextractor/pipeline.py

import os

import pandas as pd
from genome_manager import Genome
from alignment_manager import Aligner
from mutation_extractor_manager import FiveMerExtractor, MutationExtractor, MutationNormalizer, TripletExtractor
from pileup_manager import Pileup
from plot_utils import CoveragePlotter, MutationDensityPlotter, MutationSpectraPlotter
from utils import get_top_n_chromosomes, log
import pysam

class MutationExtractionPipeline:
    def __init__(self, 
                 species_list,
                 outgroup,
                 aligner="bwa", 
                 base_output_dir="../Output", 
                 no_cache = False,
                 verbose = True, 
                 run_id = None,
                 **kwargs):
        self.species_list = species_list  # list of (name, accession)
        self.outgroup = outgroup          # (name, accession)
        self.aligner = aligner
        self.run_id = run_id
        if run_id is None:
            self.run_id = '__'.join([species[0] for species in [outgroup] + species_list])
        suffix = ''
        if 'suffix' in kwargs:
            suffix = '_' + kwargs['suffix']
        self.output_dir = f"{base_output_dir}/{self.run_id}{suffix}"
        self.params = kwargs

        # Will hold references to internal data
        self.reference = None
        self.genomes = []
        self.alignments = []
        self.tree = None
        self.mutations = None
        self.verbose = verbose
        self.no_cache = no_cache

    def run(self):
        log("Starting mutation extraction pipeline...", self.verbose)

        self.download_index_and_fragment_genomes()
        self.align_species()
        self.generate_pileup()
        self.extract_mutations_and_triplets()
        self.extract_intervals()
        self.run_plots()
        
        # self.analyze_and_plot()

        print("Pipeline completed successfully.")


    def download_index_and_fragment_genomes(self):
        log("Downloading, indexing, and fragmenting genomes...", self.verbose)

        # Reference genome (outgroup)
        ref_name, ref_acc = self.outgroup
        self.reference = Genome(
            name=ref_name,
            accession=ref_acc,
            output_dir=self.output_dir,
            no_cache=self.no_cache,
            verbose=self.verbose
        )
        self.reference.download()
        self.reference.index(aligner='bwa')

        # Ingroup genomes
        for name, acc in self.species_list:
            genome = Genome(
                name=name,
                accession=acc,
                output_dir=self.output_dir,
                no_cache=self.no_cache,
                verbose=self.verbose
            )
            genome.download()
            genome.generate_fragment_fastq(length=self.params.get("fragment_length", 75), offset=self.params.get("fragment_offset", 75), force=self.no_cache)
            self.genomes.append(genome)

    def align_species(self):
        log("Aligning species to reference...", self.verbose)

        for genome in self.genomes:
            aligner = Aligner(
                species_genome=genome,
                reference_genome=self.reference,
                base_output_dir=self.output_dir,
                aligner_name=self.aligner,              
                cores=self.params.get("cores", None),
                verbose=self.params.get("verbose", True),
            )

            if self.params.get("streamed", False):
                aligner.align_streamed(
                    mapq=self.params.get("mapq", 60),
                    offset=self.params.get("offset", 75),
                    continuity=self.params.get("continuity", True),
                    max_sort_mem=self.params.get("max_samtools_mem", None)
                )
            else:
                aligner.align_disk_cached(
                    mapq=self.params.get("mapq", 60),
                    offset=self.params.get("offset", 75),
                    continuity=self.params.get("continuity", True)
                )
            self.alignments.append(aligner)


    def generate_pileup(self):
        print("Generating pileup from alignments...")

        # Ensure aligners were run and final BAMs are available
        for aligner in self.alignments:
            if not os.path.exists(aligner.final_bam):
                raise FileNotFoundError(f"BAM not found: {aligner.final_bam}")

        pileup_generator = Pileup(
            outgroup=self.reference,
            aligners=self.alignments,
            base_output_dir=self.output_dir,
            run_id=self.run_id,
            no_cache=self.no_cache,
            verbose=self.verbose
        )

        self.pileup_path = pileup_generator.generate()


    def extract_mutations_and_triplets(self):
        log("Extracting mutations from pileup...", self.verbose)
        mutation_extractor = MutationExtractor(reference=self.reference.name,
                              taxon1=self.genomes[0].name,
                              taxon2=self.genomes[1].name,
                              pileup_file=self.pileup_path,
                              output_dir=os.path.join(self.output_dir, 'Mutations'),
                              no_full_mutations=False,
                              no_cache=False,
                              verbose=self.verbose)
        mutation_extractor.extract()

        fivemer_extractor = FiveMerExtractor(reference=self.reference.name,
                              taxon1=self.genomes[0].name,
                              taxon2=self.genomes[1].name,
                              pileup_file=self.pileup_path,
                              output_dir=os.path.join(self.output_dir, 'Mutations'),
                              no_cache=False,
                              verbose=self.verbose)
        fivemer_extractor.extract()

        log("Extracting triplets from pileup...", self.verbose)
        triplet_extractor = TripletExtractor(reference=self.reference.name,
                              taxon1=self.genomes[0].name,
                              taxon2=self.genomes[1].name,
                              pileup_file=self.pileup_path,
                              output_dir=os.path.join(self.output_dir, 'Triplets'),
                              no_cache=False,
                              verbose=self.verbose)
        triplet_extractor.extract()

        normalizer = MutationNormalizer(
            input_dir=self.output_dir,
            output_dir= os.path.join(self.output_dir, "Tables"),
            verbose=True
        )
        normalizer.normalize()

    def _extract_bam_intervals(self, input_bam, output_dir, sorted=False, merge=False, no_cache=False):
            print(output_dir)
            os.makedirs(output_dir, exist_ok=True)

            base_name = os.path.basename(input_bam).rsplit(".", 1)[0]
            output_file = os.path.join(output_dir, f"{base_name}_intervals.tsv.gz")

            if os.path.exists(output_file) and not no_cache:
                print(f"Intervals already exist: {output_file}")
                return output_file

            bamfile = pysam.AlignmentFile(input_bam, "rb")

            def extract_raw_intervals(bamfile):
                intervals = []
                for read in bamfile.fetch():
                    if not read.is_unmapped:
                        chrom = bamfile.get_reference_name(read.reference_id)
                        intervals.append((chrom, read.reference_start, read.reference_end))
                return intervals

            def extract_intervals_sorted(bamfile):
                merged = []
                for read in bamfile.fetch():
                    if read.is_unmapped:
                        continue
                    chrom = bamfile.get_reference_name(read.reference_id)
                    start = read.reference_start
                    end = read.reference_end
                    if merged and merged[-1][0] == chrom and merged[-1][2] >= start:
                        merged[-1] = (chrom, merged[-1][1], max(merged[-1][2], end))
                    else:
                        merged.append((chrom, start, end))
                return merged

            def merge_intervals(intervals):
                merged = []
                for chrom, start, end in sorted(intervals):
                    if merged and merged[-1][0] == chrom and merged[-1][2] >= start:
                        merged[-1] = (chrom, merged[-1][1], max(merged[-1][2], end))
                    else:
                        merged.append((chrom, start, end))
                return merged

            intervals = (
                extract_intervals_sorted(bamfile)
                if sorted
                else merge_intervals(extract_raw_intervals(bamfile)) if merge
                else extract_raw_intervals(bamfile)
            )

            df = pd.DataFrame(intervals, columns=["chromosome", "start", "end"])
            df.to_csv(output_file, sep='\t', index=False, compression="gzip")

            print(f"Intervals written to: {output_file}")
            return output_file
    
    def extract_intervals(self):
        for bam in self.alignments:
            self._extract_bam_intervals(bam.final_bam, os.path.join(self.output_dir, 'Intervals'))    

    def run_plots(self):
        spectra_plotter = MutationSpectraPlotter()
        spectra_plotter.plot(tables_dir = os.path.join(self.output_dir, 'Tables'))
        fai_file = self.reference.fasta_path + '.fai'
        coverage_plotter = CoveragePlotter(fai_file=fai_file)
        mutation_density_plotter = MutationDensityPlotter(fai_file=fai_file)

        top_chroms = get_top_n_chromosomes(fai_file, n=3)
        print("Plotting coverage and mutation density for top chromosomes...")
        for chrom in top_chroms:
            print(f"Plotting for {chrom}...")

            coverage_plotter.plot(interval_dir=os.path.join(self.output_dir, 'Intervals'),
                                 chromosome=chrom,
                                 output_dir=os.path.join(self.output_dir, 'Plots', f"coverage_{chrom}.png"))

            mutation_density_plotter.plot(mutation_dir=os.path.join(self.output_dir, 'Mutations'),
                                 chromosome=chrom,
                                 output_dir=os.path.join(self.output_dir, 'Plots'))
            
            mutation_density_plotter.plot(mutation_dir=os.path.join(self.output_dir, 'Mutations'),
                                 chromosome=chrom,
                                 output_dir=os.path.join(self.output_dir, 'Plots'),
                                 mutation_category = r"[ACTG][C>T]G")
            

if __name__ == "__main__":
    species = [
        ("Drosophila_pseudoobscura", "GCF_009870125.1"),
        ("Drosophila_miranda", "GCF_003369915.1")
    ]
    outgroup = ("Drosophila_helvetica", "GCA_963969585.1")

    pipeline = MutationExtractionPipeline(
        species_list=species,
        outgroup=outgroup,
        aligner="bwa",
        base_output_dir="../Output_OO",
        mapq=30, 
        suffix= 'MAPQ30'
    )
    pipeline.run()

