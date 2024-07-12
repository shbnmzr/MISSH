import cplex
from Bio import SeqIO


def read_fasta(file_path):
    sequences = {}
    with open(file_path, 'r') as file:
        sequence_name = ''
        sequence_data = ''
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if sequence_name:
                    sequences[sequence_name] = sequence_data
                sequence_name = line[1:]
                sequence_data = ''
            else:
                sequence_data += line
        if sequence_name:
            sequences[sequence_name] = sequence_data
    return sequences

# Function to read seed patterns from a file
def read_seed_patterns(file_path):
    seed_patterns = []
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line:
                seed_patterns.append(line)
    return seed_patterns

# Function to generate spaced k-mers from a sequence
def generate_spaced_kmers(sequence, seed_pattern):
    kmers = set()
    seed_length = len(seed_pattern)
    for i in range(len(sequence) - seed_length + 1):
        kmer = ''.join([sequence[i + j] if seed_pattern[j] == '1' else '-' for j in range(seed_length)])
        kmers.add(kmer)
    return kmers

# Function to generate elements and subsets from sequences using spaced seeds
def generate_subsets(sequences, seed_patterns):
    elements = set()
    subsets = []
    for seq_name, sequence in sequences.items():
        for seed_pattern in seed_patterns:
            kmers = generate_spaced_kmers(sequence, seed_pattern)
            elements.update(kmers)
            subsets.append((seq_name + '_' + seed_pattern, kmers))
    return list(elements), subsets

# Function to solve the Set Cover problem using CPLEX
def solve_set_cover(elements, subsets):
    cpx = cplex.Cplex()
    cpx.set_results_stream(None)
    cpx.set_log_stream(None)

    # Create binary variables for each subset
    num_subsets = len(subsets)
    variable_names = [f"x{i}" for i in range(num_subsets)]
    cpx.variables.add(names=variable_names, types=['B'] * num_subsets)

    # Objective: Minimize the number of subsets used
    cpx.objective.set_sense(cpx.objective.sense.minimize)
    cpx.objective.set_linear([(f"x{i}", 1) for i in range(num_subsets)])

    # Constraints: Each element must be covered by at least one subset
    for element in elements:
        indices = [i for i in range(num_subsets) if element in subsets[i][1]]
        cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(ind=[f"x{i}" for i in indices], val=[1] * len(indices))],
            senses=["G"],
            rhs=[1]
        )

    cpx.solve()

    selected_subsets = [subsets[i][0] for i in range(num_subsets) if cpx.solution.get_values(f"x{i}") > 0.5]
    return selected_subsets
