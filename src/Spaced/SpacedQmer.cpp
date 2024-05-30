#include "SpacedQmer.h"
#include <algorithm>
#include <coin/CbcModel.hpp>
#include <coin/OsiClpSolverInterface.hpp>

SpacedQmer::SpacedQmer(const std::vector<std::string>& seedPatterns) : seedPatterns(seedPatterns) {
    parseSeedPatterns();
}

void SpacedQmer::parseSeedPatterns() {
    for (const auto& pattern : seedPatterns) {
        positions.push_back(parseSinglePattern(pattern));
    }
}

std::vector<int> SpacedQmer::parseSinglePattern(const std::string& pattern) const {
    std::vector<int> pos;
    for (int i = 0; i < pattern.size(); ++i) {
        if (pattern[i] == '1') {
            pos.push_back(i);
        }
    }
    return pos;
}

std::set<int> SpacedQmer::linearProgrammingSetCovering(const std::string& sequence) const {
    int numSeeds = seedPatterns.size();
    int seqLen = sequence.size();

    OsiClpSolverInterface solver;
    CoinModel model;

    // Create variables for each seed
    for (int i = 0; i < numSeeds; ++i) {
        model.setColumnBounds(i, 0, 1);
        model.setObjectiveCoefficient(i, 1);
    }

    // Map positions to seeds
    std::map<int, std::vector<int>> positionToSeeds;
    for (int i = 0; i < numSeeds; ++i) {
        for (int pos : positions[i]) {
            if (pos < seqLen) {
                positionToSeeds[pos].push_back(i);
            }
        }
    }

    // Add constraints to ensure each position is covered
    for (const auto& [pos, seeds] : positionToSeeds) {
        CoinPackedVector constraint;
        for (int seed : seeds) {
            constraint.insert(seed);
        }
        model.addRow(constraint, 1.0, COIN_DBL_MAX);
    }

    // Load problem to solver
    solver.loadFromCoinModel(model);

    // Create and solve the model
    CbcModel cbcModel(solver);
    cbcModel.solve();

    // Extract the selected seeds
    std::set<int> selectedSeeds;
    const double* solution = cbcModel.bestSolution();
    for (int i = 0; i < numSeeds; ++i) {
        if (solution[i] > 0.5) {
            selectedSeeds.insert(i);
        }
    }

    return selectedSeeds;
}

std::vector<std::string> SpacedQmer::generateSpacedQmers(const std::string& sequence) const {
    std::vector<std::string> spacedQmers;
    std::set<int> selectedSeeds = linearProgrammingSetCovering(sequence);
    for (int i : selectedSeeds) {
        for (const auto& pos : positions) {
            std::string qmer;
            for (int p : pos) {
                if (i + p < sequence.size()) {
                    qmer += sequence[i + p];
                }
            }
            spacedQmers.push_back(qmer);
        }
    }

    return spacedQmers;
}
