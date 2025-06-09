from collections import defaultdict
from ete3 import Tree
from Bio import Phylo
from io import StringIO
import sys
import os
import json

def parse_species_accession_from_newick(newick_str):
    tree = Phylo.read(StringIO(newick_str), "newick")
    species_accession_dict = {}
    for leaf in tree.get_terminals():
        if "|" in leaf.name:
            species, accession = leaf.name.split("|", 1)
            species_accession_dict[species] = accession
        else:
            print(f"Leaf name '{leaf.name}' does not contain a '|' separator.")
            sys.exit(1)

    # Guess outgroup: any subtree with a single terminal at root level
    root_clades = tree.root.clades
    outgroup = None
    for clade in root_clades:
        terminals = clade.get_terminals()
        if len(terminals) == 1:
            outgroup = terminals[0].name.split("|", 1)[0]
            break

    if outgroup is None:
        print("Could not determine a single outgroup from the Newick tree. Please ensure the tree is rooted and has a single outgroup.")
        sys.exit(1)

    return species_accession_dict, outgroup


def annotate_tree_with_indices(newick_str, outgroup_name, file_path=None):
    tree = Tree(newick_str, format=1)

    # Normalize leaf names
    for leaf in tree.iter_leaves():
        if '|' in leaf.name:
            leaf.name = leaf.name.split('|', 1)[0]

    # Sort: outgroup first, rest alphabetically
    terminals = tree.get_leaves()
    sorted_terminals = [t for t in terminals if t.name == outgroup_name] + \
                       sorted([t for t in terminals if t.name != outgroup_name], key=lambda x: x.name)

    terminal_mapping = {}
    for idx, node in enumerate(sorted_terminals):
        node.add_feature("index", idx)
        node.add_feature("custom_name", node.name)
        terminal_mapping[idx] = node.name
        terminal_mapping[node.name] = idx

    next_internal_idx = len(sorted_terminals)
    for node in tree.traverse("postorder"):
        if not node.is_leaf():
            node.add_feature("index", next_internal_idx)
            node.add_feature("custom_name", f"Node({next_internal_idx})")
            next_internal_idx += 1

    if file_path is not None:
        original_names = {}
        for node in tree.traverse():
            original_names[node] = node.name
            node.name = getattr(node, "custom_name", node.name)

        annotated_tree_path = f"{os.path.splitext(file_path)[0]}_annotated.nwk"
        tree.write(format=1, outfile=annotated_tree_path)

        for node in tree.traverse():
            node.name = original_names[node]

        mapping_path = f"{os.path.splitext(file_path)[0]}_mapping.json"
        with open(mapping_path, "w") as f:
            json.dump(terminal_mapping, f, indent=2)

        print(f"Annotated tree saved to {annotated_tree_path}")
        print(f"Terminal mapping saved to {mapping_path}")

    return tree, terminal_mapping


def save_annotated_tree(tree, path):
    original_names = {}
    for node in tree.traverse():
        original_names[node] = node.name
        node.name = getattr(node, "custom_name", node.name)

    tree.write(format=1, outfile=path)

    for node in tree.traverse():
        node.name = original_names[node]

complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

def get_complement(mutation):
    comp = [complement[nuc] if nuc in complement else nuc for nuc in mutation]
    comp[0], comp[-1] = comp[-1], comp[0]
    return ''.join(comp)

def collapse_mutations(mutation_dict):
    collapsed = defaultdict(int)
    for mutation, count in mutation_dict.items():
        if mutation[2] in {'A', 'G'}:
            collapsed[get_complement(mutation)] += int(count)
        else:
            collapsed[mutation] += int(count)
    return collapsed