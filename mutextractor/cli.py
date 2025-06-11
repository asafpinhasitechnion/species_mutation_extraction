import argparse
import json
import os
from pipeline import MutationExtractionPipeline, MultiSpeciesMutationPipeline
from run_phylip import run_phylip
from plot_utils import MutationSpectraPlotter

def main():
    parser = argparse.ArgumentParser(description="Species Mutation Extraction CLI")
    subparsers = parser.add_subparsers(dest="command")

    # === Single-Pipeline ===
    single = subparsers.add_parser("run_single", help="Run 3-species pipeline (outgroup + 2)")
    single.add_argument("--outgroup", nargs=2, metavar=("NAME", "ACCESSION"), required=True)
    single.add_argument("--species", nargs=4, metavar=("NAME1", "ACC1", "NAME2", "ACC2"), required=True)
    single.add_argument("--output", required=True)
    single.add_argument("--no-cache", action="store_true")
    single.add_argument("--verbose", action="store_true")
    single.add_argument("--suffix", default=None)
    single.add_argument("--aligner", default="bwa")
    single.add_argument("--aligner-cmd", default=None)
    single.add_argument("--streamed", action="store_true")
    single.add_argument("--mapq", type=int, default=60)
    single.add_argument("--low-mapq", type=int, default=1)
    single.add_argument("--cores", type=int, default=None)

    # === Multi-Species ===
    multi = subparsers.add_parser("run_multi", help="Run multi-species pipeline from Newick")
    multi.add_argument("--tree", required=True)
    multi.add_argument("--output", required=True)
    multi.add_argument("--no-cache", action="store_true")
    multi.add_argument("--verbose", action="store_true")
    multi.add_argument("--outgroup", default=None)
    multi.add_argument("--aligner", default="bwa")
    multi.add_argument("--aligner-cmd", default=None)
    multi.add_argument("--streamed", action="store_true")
    multi.add_argument("--mapq", type=int, default=60)
    multi.add_argument("--low-mapq", type=int, default=1)
    multi.add_argument("--cores", type=int, default=None)

    # === Plotting ===
    plot = subparsers.add_parser("plot", help="Generate plots from output Tables directory")
    plot.add_argument("--tables", required=True)

    # === Run PHYLIP ===
    phylip = subparsers.add_parser("run_phylip", help="Run PHYLIP on mutation matrix")
    phylip.add_argument("--df", required=True, help="Path to matching_bases.csv.gz")
    phylip.add_argument("--tree", required=True, help="Path to annotated_tree.nwk")
    phylip.add_argument("--mapping", required=True, help="Path to species_mapping.json")
    phylip.add_argument("--command", default="dnapars", help="PHYLIP command: dnapars, dnapenny, etc.")
    phylip.add_argument("--prefix", default="phylip_run")
    phylip.add_argument("--input-string", default="Y\n")

    args = parser.parse_args()

    if args.command == "run_single":
        pipeline = MutationExtractionPipeline(
            species_list=[(args.species[0], args.species[1]), (args.species[2], args.species[3])],
            outgroup=(args.outgroup[0], args.outgroup[1]),
            base_output_dir=args.output,
            no_cache=args.no_cache,
            verbose=args.verbose,
            suffix=args.suffix,
            aligner=args.aligner,
            streamed=args.streamed,
            mapq=args.mapq,
            low_mapq=args.low_mapq,
            cores=args.cores,
        )
        pipeline.run()

    elif args.command == "run_multi":
        newick_tree = open(args.tree).read()
        pipeline = MultiSpeciesMutationPipeline(
            newick_tree=newick_tree,
            base_output_dir=args.output,
            outgroup=args.outgroup,
            no_cache=args.no_cache,
            verbose=args.verbose,
            aligner=args.aligner,
            streamed=args.streamed,
            mapq=args.mapq,
            low_mapq=args.low_mapq,
            cores=args.cores,
        )
        pipeline.run()

    elif args.command == "plot":
        plotter = MutationSpectraPlotter()
        plotter.plot(tables_dir=args.tables)

    elif args.command == "run_phylip":
        with open(args.mapping) as f:
            mapping = json.load(f)
        run_phylip(
            command=args.command,
            df_path=args.df,
            tree_path=args.tree,
            output_dir=os.path.dirname(args.df),
            prefix=args.prefix,
            input_string=args.input_string,
            mapping=mapping
        )

    else:
        parser.print_help()


if __name__ == "__main__":
    main()
