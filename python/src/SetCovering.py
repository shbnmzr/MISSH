import cplex
from Bio import SeqIO


# Function to read sequences from a FASTA file
# Function to read sequences from a FASTA file
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

# Function to generate k-mers from a sequence
def generate_kmers(sequence, k):
    kmers = set()
    for i in range(len(sequence) - k + 1):
        kmers.add(sequence[i:i+k])
    return kmers


# Function to generate elements and subsets from sequences
def generate_subsets(sequences, k):
    elements = set()
    subsets = []
    for seq_name, sequence in sequences.items():
        kmers = generate_kmers(sequence, k)
        elements.update(kmers)
        subsets.append((seq_name, kmers))
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
        print(f"Element {element} is covered by subsets: {indices}")
        cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(ind=[f"x{i}" for i in indices], val=[1] * len(indices))],
            senses=["G"],
            rhs=[1]
        )

    cpx.solve()

    selected_subsets = [subsets[i][0] for i in range(num_subsets) if cpx.solution.get_values(f"x{i}") > 0.5]
    print(f"Selected subsets: {selected_subsets}")
    return selected_subsets


# def solve_set_cover(elements, subsets):
#     cpx = cplex.Cplex()
#     cpx.set_log_stream(None)
#     cpx.set_results_stream(None)
#
#     cpx.variables.add(names=[f"x{i}" for i in range(len(subsets))], types=["B"] * len(subsets))
#
#     for element in elements:
#         indices = [i for i in range(len(subsets)) if element in subsets[i][1]]
#         cpx.linear_constraints.add(
#             lin_expr=[cplex.SparsePair(ind=[f"x{i}" for i in indices], val=[1]*len(indices))],
#             senses=["G"],
#             rhs=[1]
#         )
#
#     cpx.objective.set_sense(cpx.objective.sense.minimize)
#     cpx.solve()
#
#     solution = cpx.solution.get_values()
#     selected_subsets = [subsets[i][0] for i in range(len(solution)) if solution[i] > 0.5]
#
#     return selected_subsets
#
#
# def read_sequences_from_files(file_paths):
#     sequences = {}
#     for file_path in file_paths:
#         with open(file_path, "r") as file:
#             for record in SeqIO.parse(file, "fasta"):
#                 sequences[record.id] = str(record.seq)
#     return sequences
#
#
# def generate_kmers(sequence, k):
#     return {sequence[i:i+k] for i in range(len(sequence) - k + 1)}
#
#
# def generate_subsets(sequences, k):
#     elements = set()
#     subsets = []
#
#     for seq_name, sequence in sequences.items():
#         kmers = generate_kmers(sequence, k)
#         elements.update(kmers)
#         subsets.append((seq_name, kmers))
#
#     return list(elements), subsets
