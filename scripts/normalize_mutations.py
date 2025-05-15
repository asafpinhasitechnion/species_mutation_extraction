import os
import json
import sys
import pandas as pd
from collections import defaultdict
import re

# Global variables
complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
valid_bases = {'A', 'G', 'C', 'T'}
mutation_characters = {'A', 'G', 'C', 'T', '[',']','>'}

# Folders
mutation_folder = '../Mutation_Dicts'
triplets_folder = '../Triplets_Dicts'
tables_folder = '../Tables'

def collapse_triplets(triplet_dict):
    collapsed_triplets = defaultdict(int)
    for triplet, count in triplet_dict.items():
        if triplet[1] in ['G', 'A']:
            collapsed_triplets[''.join(reversed([complement[ch] for ch in triplet]))] += count
        else:
            collapsed_triplets[triplet] += count
    return collapsed_triplets

def get_complement(mutation):
    """Get the complement mutation."""
    comp_mutation = [complement[nuc] if nuc in complement else nuc for nuc in mutation]
    comp_mutation[0], comp_mutation[-1] = comp_mutation[-1], comp_mutation[0]
    return ''.join(comp_mutation)

def collapse_mutations(mutation_dict):
    """Collapse the mutations properly to reach 96 mutational categories."""
    collapsed_mutations = defaultdict(int)
    for mutation, count in mutation_dict.items():
        if mutation[2] in {'A', 'G'}:
            collapsed_mutations[get_complement(mutation)] += int(count)
        else:
            collapsed_mutations[mutation] += int(count)
    return collapsed_mutations


def filter_triplets_dict(species_dict):
    """Remove keys that contain 'N'."""
    return {key: value for key, value in species_dict.items() if all(nuc in valid_bases for nuc in key)}


# Define the regex pattern for the mutation format
mutation_pattern = re.compile(r"^[ACGT]\[[ACGT]>[ACGT]\][ACGT]$")

def filter_mutations_dict(species_dict):
    return {key: value for key, value in species_dict.items() if mutation_pattern.match(key)}

def normalize_mutation_count_dict(mutation_count, triplet_count):
    """Normalize mutation counts by triplet count."""
    return {key: mutation_count[key] / triplet_count[f"{key[0]}{key[2]}{key[-1]}"] if triplet_count[f"{key[0]}{key[2]}{key[-1]}"] > 0 else 0 
            for key, value in mutation_count.items()}

def denormalize_mutation_count_dict(mutation_count, triplet_count):
    """Denormalize mutation counts by triplet count."""
    return {key: mutation_count[key] * triplet_count[f"{key[0]}{key[2]}{key[-1]}"] 
            for key, value in mutation_count.items()}

def scale_mutations(mutation_count, scaling=10000):
    """Scale mutation counts."""
    total_sum = sum(mutation_count.values())
    return {k: int(round(v / total_sum * scaling)) if total_sum > 0 else 0 for k, v in mutation_count.items()}

def create_folder_if_not_exists(folder):
    """Create the folder if it doesn't exist."""
    if not os.path.exists(folder):
        os.mkdir(folder)

def process_file_paths(reference, species1, species2):
    """Generate file paths and ensure required folders exist."""
    file_name = f'{reference}_{species1}({species2}).json'
    raw_mutation_path = os.path.join(mutation_folder, 'Raw_Mutations', file_name)
    raw_triplet_path = os.path.join(triplets_folder, 'Raw_Triplets', file_name)
    
    return raw_mutation_path, raw_triplet_path, file_name

def save_json(data, folder, file_name):
    """Save the given data as a JSON file."""
    with open(os.path.join(folder, file_name), "w") as f:
        json.dump(data, f)

def load_json(file_path):
    """Load data from a JSON file."""
    if os.path.exists(file_path):
        with open(file_path, "r") as f:
            return json.load(f)
    else:
        print(f"File not found: {file_path}")
    return None

# Create output folders
create_folder_if_not_exists(os.path.join(mutation_folder, 'Collapsed_Mutations'))
create_folder_if_not_exists(os.path.join(triplets_folder, 'Collapsed_Triplets'))
create_folder_if_not_exists(os.path.join(mutation_folder, 'Scaled_Mutations'))
create_folder_if_not_exists(os.path.join(mutation_folder, 'Normalized_Mutations'))
create_folder_if_not_exists(os.path.join(mutation_folder, 'Human_Normalized_Mutations'))
create_folder_if_not_exists(os.path.join(mutation_folder, 'Scaled_Normalized_Mutations'))

species_raw_mutation_dict = {}
species_normalized_mutation_dict = {}
species_human_normalized_mutation_dict = {}
species_mutation_dict = {}
species_triplets_dict = {}
species_name_to_file = {}

_, human_raw_triplet_path, _ = process_file_paths('Gorilla_gorilla_gorilla', 'Homo_sapiens', 'Pan_troglodytes')
human_triplets = load_json(human_raw_triplet_path)
human_triplets = filter_triplets_dict(human_triplets)
human_triplets = collapse_triplets(human_triplets)



with open('multiz_normalize_mutations_arguments.txt', 'r') as file:
    for line in file:
        reference,species1,species2 = line[:-1].split(',')
        # Process file paths and create folders
        raw_mutation_path, raw_triplet_path, file_name = process_file_paths(reference, species1, species2)
        
        # Load mutations and triplets
        species_mutations = load_json(raw_mutation_path)
        species_triplets = load_json(raw_triplet_path)
        
        if not species_mutations or not species_triplets:
            continue
        
        # Filter the loaded data
        species_mutations = filter_mutations_dict(species_mutations)
        species_triplets = filter_triplets_dict(species_triplets)
        
        # Collapse and save mutations and triplets
        collapsed_species_mutations = collapse_mutations(species_mutations)
        save_json(collapsed_species_mutations, os.path.join(mutation_folder, 'Collapsed_Mutations'), file_name)
        
        collapsed_species_triplets = collapse_triplets(species_triplets)
        save_json(collapsed_species_triplets, os.path.join(triplets_folder, 'Collapsed_Triplets'), file_name)

        species_name = species1
        i = 1
        while species_name in species_normalized_mutation_dict:
            species_name = species1 + '_' + str(i)
            i += 1
        species_name = file_name
        species_name_to_file[species_name] = file_name
        # Scale the mutations and save
        scaled_mutations = scale_mutations(collapsed_species_mutations)
        save_json(scaled_mutations, os.path.join(mutation_folder, 'Scaled_Mutations'), file_name)
        species_mutation_dict[species_name]=scaled_mutations
        
        # Normalize the mutations and save
        normalized_species_mutations = normalize_mutation_count_dict(collapsed_species_mutations, collapsed_species_triplets)
        save_json(normalized_species_mutations, os.path.join(mutation_folder, 'Normalized_Mutations'), file_name)

        # Scale normalized mutations and save
        scaled_normalized_mutations = scale_mutations(normalized_species_mutations)
        save_json(scaled_normalized_mutations, os.path.join(mutation_folder, 'Scaled_Normalized_Mutations'), file_name)

        human_normalized_species_mutations = denormalize_mutation_count_dict(normalized_species_mutations, human_triplets)
        scaled_human_normalized_mutations = scale_mutations(human_normalized_species_mutations)
        save_json(scaled_human_normalized_mutations, os.path.join(mutation_folder, 'Human_Normalized_Mutations'), file_name)

        species_normalized_mutation_dict[species_name]=scaled_normalized_mutations
        species_triplets_dict[species_name] = collapsed_species_triplets
        species_raw_mutation_dict[species_name] = collapsed_species_mutations
        species_human_normalized_mutation_dict[species_name] = scaled_human_normalized_mutations
        # print(species_name)

if not os.path.exists(tables_folder):
    os.mkdir(tables_folder)


with open('Multiz_name_to_file.json', 'w') as f:
    json.dump(species_name_to_file, f)

raw_mutations_df = pd.DataFrame(species_raw_mutation_dict)
print(raw_mutations_df.shape)
raw_mutations_df.to_csv(os.path.join(tables_folder, f'{raw_mutations_df.shape[1]}_raw_mutations.tsv'), sep = '\t')
pd.DataFrame(species_normalized_mutation_dict).to_csv(os.path.join(tables_folder, f'{raw_mutations_df.shape[1]}_normalized_mutations.tsv'), sep = '\t')
pd.DataFrame(species_human_normalized_mutation_dict).to_csv(os.path.join(tables_folder, f'{raw_mutations_df.shape[1]}_human_normalized_mutations.tsv'), sep = '\t')
pd.DataFrame(species_mutation_dict).to_csv(os.path.join(tables_folder, f'{raw_mutations_df.shape[1]}_scaled_mutations.tsv'), sep = '\t')
pd.DataFrame(species_triplets_dict).to_csv(os.path.join(tables_folder, f'{raw_mutations_df.shape[1]}_triplets.tsv'), sep = '\t')