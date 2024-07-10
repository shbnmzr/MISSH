import cplex
from Bio import SeqIO


def solve_set_cover(elements, subsets):
    cpx = cplex.Cplex()
    cpx.set_log_stream(None)
    cpx.set_results_stream(None)

    cpx.variables.add(names=[f"x{i}" for i in range(len(subsets))], types=["B"] * len(subsets))

    for element in elements:
        indices = [i for i in range(len(subsets)) if element in subsets[i][1]]
        cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(ind=[f"x{i}" for i in indices], val=[1]*len(indices))],
            senses=["G"],
            rhs=[1]
        )

    cpx.objective.set_sense(cpx.objective.sense.minimize)
    cpx.solve()

    solution = cpx.solution.get_values()
    selected_subsets = [subsets[i][0] for i in range(len(solution)) if solution[i] > 0.5]

    return selected_subsets


def read_sequences_from_files(file_paths):
    sequences = {}
    for file_path in file_paths:
        with open(file_path, "r") as file:
            for record in SeqIO.parse(file, "fasta"):
                sequences[record.id] = str(record.seq)
    return sequences


def generate_kmers(sequence, k):
    return {sequence[i:i+k] for i in range(len(sequence) - k + 1)}


def generate_subsets(sequences, k):
    elements = set()
    subsets = []

    for seq_name, sequence in sequences.items():
        kmers = generate_kmers(sequence, k)
        elements.update(kmers)
        subsets.append((seq_name, kmers))

    return list(elements), subsets
