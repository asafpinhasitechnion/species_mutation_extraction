import sys
import os
import subprocess
from pathlib import Path
import json
from datetime import datetime
from Bio import Phylo
from io import StringIO
from ete3 import Tree


def run_cmd(cmd, shell=False):
    print(f"➡️ Running: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
    result = subprocess.run(cmd, shell=shell)
    if result.returncode != 0:
        print(f"❌ Command failed: {cmd}")
        sys.exit(result.returncode)


def write_metadata(args_dict, output_dir):
    metadata = {
        "timestamp": datetime.now().isoformat(),
        "arguments": args_dict
    }

    metadata_path = Path(output_dir) / "metadata.json"
    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=2)

    print(f"📝 Metadata written to {metadata_path}")


def parse_args():
    # if len(sys.argv) < 3:
    #     print("Usage: run_pipeline.py <newick_species_tree> <run_id> [optional args...]")
    #     sys.exit(1)

    # newick_tree, run_id = sys.argv[1:3]
    # optional_args = sys.argv[3:]


    newick_tree = "(Leptophobia_aripa|GCA_951799465.1, (Pieris_brassicae|GCF_905147105.1, (Pieris_napi|GCF_905475465.1, (Pieris_rapae|GCF_905147795.1, Pieris_mannii|GCA_028984075.1))));"
    run_id = 'test_run_mutiple_species'
    optional_args = ['--mapq', '1', '--no-cache']

    # === Argument groups ===
    download_args, index_args = [], []
    align_filter_args, pileup_args, mutation_args, triplet_args, interval_args = [], [], [], [], []
    global GENOMIC_PLOTS
    GENOMIC_PLOTS = False
    reference = None

    i = 0
    while i < len(optional_args):
        arg = optional_args[i]
        if arg in ["--remove-temp"]:
            align_filter_args.append(arg)
            pileup_args.append(arg)
            i += 1
        elif arg == "--no-cache":
            # Add to all affected stages
            download_args.append(arg)
            index_args.append(arg)
            align_filter_args.append(arg)
            pileup_args.append(arg)
            mutation_args.append(arg)
            triplet_args.append(arg)
            interval_args.append(arg)
            i += 1
        elif arg == "--mapq":
            align_filter_args.extend([arg, optional_args[i + 1]])
            i += 2
        elif arg == "--no-full-mutations":
            mutation_args.append(arg)
            i += 1
        elif arg == "--aligner":
            align_filter_args.extend([arg, optional_args[i + 1]])
            i += 2
        elif arg == "--aligner-cmd":
            align_filter_args.extend([arg, optional_args[i + 1]])
            i += 2
        elif arg == '--reference':
            reference = optional_args[i + 1]
            i += 2
        elif arg == '--genomic-position-plots':
            GENOMIC_PLOTS = True
            i += 1
        else:
            print(f"❗ Unknown argument: {arg}")
            sys.exit(1)

    return {
        "outgroup": reference,
        "newick_tree": newick_tree,
        "run_id": run_id,
        "download_args": download_args,
        "index_args": index_args,
        "align_filter_args": align_filter_args,
        "pileup_args": pileup_args,
        "mutation_args": mutation_args,
        "triplet_args": triplet_args,
        "interval_args": interval_args,
    }

def get_top_n_chromosomes(fai_path, n=2):
    chroms = []
    with open(fai_path) as f:
        for line in f:
            fields = line.strip().split('\t')
            chroms.append((fields[0], int(fields[1])))
    chroms.sort(key=lambda x: -x[1])
    return [c[0] for c in chroms[:n]]

def parse_species_accession_from_newick(newick_str):
    tree = Phylo.read(StringIO(newick_str), "newick")
    species_accession_dict = {}
    for leaf in tree.get_terminals():
        # Expecting leaf names in the format: species|accession
        if "|" in leaf.name:
            species, accession = leaf.name.split("|", 1)
            species_accession_dict[species] = accession
        else:
            print(f"❗ Leaf name '{leaf.name}' does not contain a '|' separator.")
            sys.exit(1)

    # Extract outgroup from the tree, or raise an error if no single outgroup exists
    root_clades = tree.root.clades
    outgroup = None
    for clade in root_clades:
        terminals = clade.get_terminals()
        if len(terminals) == 1:
            outgroup = terminals[0].name.split("|", 1)[0]
            break

    if outgroup is None:
        print("❗ Could not determine a single outgroup from the Newick tree. Please ensure the tree is rooted and has a single outgroup.")
        sys.exit(1)

    return species_accession_dict, outgroup


def annotate_tree_with_indices(newick_str, outgroup_name, file_path=None):
    """
    Traverses the Newick tree using ete3, annotates each node with an index.
    Terminal nodes: 0 - outgroup, 1, 2, ... for the rest.
    Internal nodes: Node(idx).
    Returns:
      - tree (with .index and .custom_name on each node)
      - terminal_mapping: {species_name <-> idx} (bidirectional for terminals only)
    If file_path is provided, saves the annotated tree (as newick) and the mapping (as json)
    with suffixes '_annotated.nwk' and '_mapping.json'.
    """
    # Step 1: Load tree
    tree = Tree(newick_str, format=1)

    # Step 2: Clean names (strip suffix after "|")
    for leaf in tree.iter_leaves():
        if '|' in leaf.name:
            leaf.name = leaf.name.split('|', 1)[0]

    # Step 3: Sort terminals (outgroup first)
    terminals = tree.get_leaves()
    sorted_terminals = [t for t in terminals if t.name == outgroup_name] + \
                       sorted([t for t in terminals if t.name != outgroup_name], key=lambda x: x.name)

    # Step 4: Assign terminal indices and build mapping
    terminal_mapping = {}
    for idx, node in enumerate(sorted_terminals):
        node.add_feature("index", idx)
        node.add_feature("custom_name", node.name)
        species_name = node.name
        terminal_mapping[idx] = species_name
        terminal_mapping[species_name] = idx

    # Step 5: Assign internal node indices
    next_internal_idx = len(sorted_terminals)
    for node in tree.traverse("postorder"):
        if not node.is_leaf():
            node.add_feature("index", next_internal_idx)
            node.add_feature("custom_name", f"Node({next_internal_idx})")
            next_internal_idx += 1

    # Step 6: Output to files (if needed)
    if file_path is not None:
        # Temporarily rename nodes to custom_name for writing
        original_names = {}
        for node in tree.traverse():
            original_names[node] = node.name
            node.name = getattr(node, "custom_name", node.name)

        annotated_tree_path = f"{os.path.splitext(file_path)[0]}_annotated.nwk"
        tree.write(format=1, outfile=annotated_tree_path)

        # Restore original names
        for node in tree.traverse():
            node.name = original_names[node]

        # Save mapping as JSON
        mapping_path = f"{os.path.splitext(file_path)[0]}_mapping.json"
        with open(mapping_path, "w") as f:
            json.dump(terminal_mapping, f, indent=2)

        print(f"📝 Annotated tree saved to {annotated_tree_path}")
        print(f"📝 Terminal mapping saved to {mapping_path}")

    return tree, terminal_mapping


def main():
    args = parse_args()

    run_id = args['run_id']
    base_output_dir = Path("../Output") / run_id
    base_output_dir.mkdir(parents=True, exist_ok=True)
    newick_tree = args["newick_tree"]
    species_accession_dict, outgroup = parse_species_accession_from_newick(newick_tree)
    if args['outgroup']:
        outgroup = args['outgroup']

    print(f"📁 Run ID: {run_id}")
    print(f"📂 Base output directory: {base_output_dir}")
    write_metadata(args, base_output_dir)

    # === GENOME DOWNLOADS ===
    print("⬇️ Downloading genomes...")
    for species, accession in species_accession_dict.items():
        run_cmd(["bash", "download_genome.sh", species, accession, str(base_output_dir)] + args["download_args"])

    # REFERENCE INDEXING
    print(f"🧬 Indexing outgroup genome: {outgroup}")
    run_cmd(["bash", "index_reference_genome.sh", outgroup, str(base_output_dir)] + args["index_args"])

    # ALIGNMENTS
    for species, accession in species_accession_dict.items():
        print(f"🔗 Aligning {species} to {outgroup}")
        run_cmd(["bash", "customizable_align_and_filter.sh", species, outgroup, str(base_output_dir)] + args["align_filter_args"])

    print(f"✅ Alignment and filtering complete for {run_id}")

    tree, terminal_mapping = annotate_tree_with_indices(newick_tree, outgroup)
    # Extract terminal indices in order, skipping outgroup index 0
    ordered_taxa = [
        terminal_mapping[idx]
        for idx in sorted(k for k in terminal_mapping if isinstance(k, int) and k != 0)
    ]

    # PILEUP
    print("📊 Creating pileup...")
    run_cmd(["bash", "create_multiple_species_pileup.sh", outgroup, str(base_output_dir)] + ordered_taxa + args["pileup_args"])

    # # MUTATIONS
    # print("🧪 Extracting mutations...")
    # run_cmd(["python3", "extract_multiple_species_mutations.py", args["out_name"], args["t1_name"], args["t2_name"],
    #          "--pileup-dir", str(base_output_dir),
    #          "--output-dir", str(base_output_dir / "Mutations")] + args["mutation_args"])
    

if __name__ == "__main__":
    main()
