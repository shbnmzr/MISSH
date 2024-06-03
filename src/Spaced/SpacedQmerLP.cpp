#include "SpacedQmerLP.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>

using namespace std;

// Constructor definitions
SpacedQmerLP::SpacedQmerLP() : SpacedQmer() {}

SpacedQmerLP::SpacedQmerLP(string spaced_qmer, size_t numprev) : SpacedQmer(spaced_qmer, numprev) {}

// LP-based spaced seed hashing algorithm
vector<Seed> SpacedQmerLP::lp_spaced_seed(const vector<int>& data, int seed_length) {
    int n = data.size();
    vector<vector<int>> A(n, vector<int>(n, 0));
    for (int j = 0; j < n; ++j) {
        for (int i = j; i < min(n, j + seed_length); ++i) {
            A[i][j] = 1;
        }
    }

    vector<int> c(n, 1);
    vector<int> x(n, 0);

    for (int i = 0; i < n; ++i) {
        bool can_place_seed = true;
        for (int j = i; j < i + seed_length && j < n; ++j) {
            if (x[j] == 1) {
                can_place_seed = false;
                break;
            }
        }
        if (can_place_seed) {
            for (int j = i; j < i + seed_length && j < n; ++j) {
                x[j] = 1;
            }
        }
    }

    vector<Seed> seeds;
    for (int i = 0; i < n; ++i) {
        if (x[i] == 1) {
            seeds.push_back({i, seed_length});
        }
    }

    return seeds;
}
