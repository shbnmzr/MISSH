#include "SetCoveringSolver.h"
#include <ortools/linear_solver/linear_solver.h>

SetCoveringSolver::SetCoveringSolver(const std::set<int>& universal_set, const std::vector<std::set<int>>& subsets)
    : universal_set_(universal_set), subsets_(subsets) {}

std::vector<int> SetCoveringSolver::solve() {
    using namespace operations_research;
    MPSolver solver("SetCoveringSolver", MPSolver::CBC_MIXED_INTEGER_PROGRAMMING);
    const int num_subsets = subsets_.size();

    // Create variables
    std::vector<const MPVariable*> x(num_subsets);
    for (int i = 0; i < num_subsets; ++i) {
        x[i] = solver.MakeBoolVar("x" + std::to_string(i));
    }

    // Create constraints
    for (int element : universal_set_) {
        LinearExpr constraint;
        for (int i = 0; i < num_subsets; ++i) {
            if (subsets_[i].count(element) > 0) {
                constraint += x[i];
            }
        }
        solver.MakeRowConstraint(constraint >= 1.0);
    }

    // Create objective function
    LinearExpr objective;
    for (int i = 0; i < num_subsets; ++i) {
        objective += x[i];
    }
    solver.Minimize(objective);

    // Solve the problem
    const MPSolver::ResultStatus result_status = solver.Solve();

    if (result_status != MPSolver::OPTIMAL) {
        throw std::runtime_error("The problem does not have an optimal solution.");
    }

    // Retrieve the solution
    std::vector<int> solution;
    for (int i = 0; i < num_subsets; ++i) {
        if (x[i]->solution_value() > 0.5) {
            solution.push_back(i);
        }
    }

    return solution;
}
