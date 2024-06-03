#ifndef SRC_HASH_SPACEDQMERLP_H_
#define SRC_HASH_SPACEDQMERLP_H_

#include "SpacedQmer.h"
#include <vector>

using namespace std;

class SpacedQmerLP : public SpacedQmer {
public:
    SpacedQmerLP();
    SpacedQmerLP(string spaced_qmer, size_t numprev);

    // New method for LP-based spaced seed hashing
    vector<Seed> lp_spaced_seed(const vector<int>& data, int seed_length);
};

#endif /* SRC_HASH_SPACEDQMERLP_H_ */
