import cplex
import sys
from src.SetCovering import *
from src.functions import log_output
from src.functions import is_executable


def main():
    is_executable()

    fasta_file_path = sys.argv[1]
    seed_file_path = sys.argv[2]
    # output_file_path = sys.argv[3]

    sequences = read_fasta(fasta_file_path)
    seed_patterns = read_seed_patterns(seed_file_path)
    elements, subsets = generate_subsets(sequences, seed_patterns)
    selected_subsets = solve_set_cover(elements, subsets)
    # subsets = [(name, list(subset)) for name, subset in subsets]
    # elements = list(elements)

    # data = {
    #     "sequences": sequences,
    #     "seed_patterns": seed_patterns,
    #     "elements": elements,
    #     "subsets": subsets,
    #     "selected_subsets": selected_subsets
    # }
    # log_output(output_file_path, data)


if __name__ == "__main__":
    main()
