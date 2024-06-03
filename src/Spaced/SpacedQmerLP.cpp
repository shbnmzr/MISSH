#include "SpacedQmerLP.h"
#include <algorithm>
#include <vector>

SpacedQmerLP::SpacedQmerLP() : SpacedQmer() {}

SpacedQmerLP::SpacedQmerLP(std::string spaced_qmer, size_t numprev) : SpacedQmer(spaced_qmer, numprev) {}

std::vector<Position> SpacedQmerLP::lp_spaced_seed(const std::vector<int>& data, int seed_length) const {
    int n = data.size();
    std::vector<std::vector<int>> A(n, std::vector<int>(n, 0));
    std::vector<int> x(n, 0);
    std::vector<int> b(n, 1);

    for (int i = 0; i < n - seed_length + 1; ++i) {
        for (int j = i; j < i + seed_length; ++j) {
            A[i][j] = 1;
        }
    }

    while (true) {
        bool found = false;
        for (int j = 0; j < n; ++j) {
            int sum = 0;
            for (int i = 0; i < n; ++i) {
                sum += A[i][j] * x[i];
            }
            if (sum < b[j]) {
                x[j] = 1;
                found = true;
                break;
            }
        }
        if (!found) break;
    }

    std::vector<Position> seeds;
    for (int i = 0; i < n; ++i) {
        if (x[i] == 1) {
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
