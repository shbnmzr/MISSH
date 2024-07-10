import cplex
from src.SetCovering import *
from itertools import combinations


def main():

    print(f'CPLEX Version: {cplex.__version__}')

    # Define the k value for k-mers
    k = 3

    # File paths for the input files
    file_paths = ["./input/paired.fna.1",
                  "./input/paired.fna.2"]

    # Read sequences from files
    sequences = read_sequences_from_files(file_paths)

    # Generate elements and subsets
    elements, subsets = generate_subsets(sequences, k)

    # Solve the set cover problem
    selected_subsets = solve_set_cover(elements, subsets)

    print("Sequences read from files:")
    for sequence in sequences:
        print(sequence)

    print("\nGenerated elements (k-mers):")
    print(elements)

    print("\nGenerated subsets:")
    for subset in subsets:
        print(subset)

    print("Selected subsets:", selected_subsets)


if __name__ == "__main__":
    main()
