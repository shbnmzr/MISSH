#include "SpacedQmerLP.h"
#include <algorithm>
#include <vector>
#include <iostream>
#include <limits>
#include <cmath>
#include "ortools/linear_solver/linear_solver.h"

SpacedQmerLP::SpacedQmerLP() : SpacedQmer() {}

SpacedQmerLP::SpacedQmerLP(std::string spaced_qmer, size_t numprev) : SpacedQmer(spaced_qmer, numprev) {}

std::vector<Position> SpacedQmerLP::lp_spaced_seed(const std::vector<int>& data, int seed_length) const {
    using namespace operations_research;

    int n = data.size();
    std::unique_ptr<MPSolver> solver(MPSolver::CreateSolver("GLOP"));
    if (!solver) {
        std::cerr << "Linear solver not created." << std::endl;
        return {};
    }

    // Variables: x[i] will be 1 if position i is included in the seed, 0 otherwise.
    std::vector<MPVariable*> vars(n);
    for (int i = 0; i < n; ++i) {
        vars[i] = solver->MakeBoolVar("x_" + std::to_string(i));
    }

    // Constraints: Ensure every position is covered by at least one seed of length 'seed_length'.
    for (int i = 0; i < n - seed_length + 1; ++i) {
        LinearExpr constraint_expr;
        for (int j = i; j < i + seed_length && j < n; ++j) {
            constraint_expr += vars[j];
        }
        solver->MakeRowConstraint(constraint_expr >= 1.0);
    }

    // Objective: Minimize the number of positions used in the seed.
    LinearExpr objective;
    for (int i = 0; i < n; ++i) {
        objective += vars[i];
    }
    solver->MutableObjective()->MinimizeLinearExpr(objective);

    // Solve the problem.
    const MPSolver::ResultStatus result_status = solver->Solve();
    if (result_status != MPSolver::OPTIMAL) {
        std::cerr << "The problem does not have an optimal solution!" << std::endl;
        return {};
    }

    // Extract the solution.
    std::vector<Position> seeds;
    for (int i = 0; i < n; ++i) {
        if (vars[i]->solution_value() > 0.5) {
            Position seed;
            for (int j = i; j < i + seed_length && j < n; ++j) {
                seed.push_back(j);
            }
            seeds.push_back(seed);
        }
    }

    return seeds;
}

// Function to compute hashes using the LP-based method
void GetHashes_with_LP(const std::string& sequence, const SpacedQmerLP& spaced, Hash_Err_V& vHash, hash_type (*CharToInt)(char)) {
    std::vector<int> data(sequence.size());
    for (size_t i = 0; i < sequence.size(); ++i) {
        data[i] = CharToInt(sequence[i]);
    }
    int seed_length = spaced.GetQ();
    std::vector<Position> positions = spaced.lp_spaced_seed(data, seed_length);
    vHash.clear();
    for (const auto& pos : positions) {
        Hash_Err hash_err;
        GetHashFromPosOne(sequence, 0, pos, hash_err, CharToInt);
        vHash.push_back(hash_err);
    }
}
