#ifndef SRC_HASH_SPACEDQMERLP_H_
#define SRC_HASH_SPACEDQMERLP_H_

#include "SpacedQmer.h"
#include "../Hash/HashType.h"
#include "../Hash/HashFunction.h"
#include <vector>

// Include OR-Tools headers
#include "ortools/linear_solver/linear_solver.h"

// Class for LP-based spaced seed hashing
class SpacedQmerLP : public SpacedQmer {
public:
    SpacedQmerLP();  // Default constructor
    SpacedQmerLP(std::string spaced_qmer, size_t numprev);  // Parameterized constructor

    // LP-based spaced seed hashing algorithm
    std::vector<Position> lp_spaced_seed(const std::vector<int>& data, int seed_length) const;
};

void GetHashes_with_LP(const std::string& sequence, const SpacedQmerLP& spaced, Hash_Err_V& vHash, hash_type (*CharToInt)(char));

#endif /* SRC_HASH_SPACEDQMERLP_H_ */
