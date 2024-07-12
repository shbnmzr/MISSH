import json
import cplex
import sys
from src.SetCovering import *


def main():
    fasta_file_path = sys.argv[1]
    seed_file_path = sys.argv[2]
    output_file_path = sys.argv[3]

    sequences = read_fasta(fasta_file_path)
    seed_patterns = read_seed_patterns(seed_file_path)

    elements, subsets = generate_subsets(sequences, seed_patterns)
    selected_subsets = solve_set_cover(elements, subsets)

    # Convert sets to lists for JSON serialization
    subsets = [(name, list(subset)) for name, subset in subsets]
    elements = list(elements)

    output_data = {
        "sequences": sequences,
        "seed_patterns": seed_patterns,
        "elements": elements,
        "subsets": subsets,
        "selected_subsets": selected_subsets
    }

    with open(output_file_path, 'w') as output_file:
        json.dump(output_data, output_file, indent=4)


if __name__ == "__main__":
    main()
