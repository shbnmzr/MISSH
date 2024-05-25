/*
 * MultiSpacedQmer.h
 *
 *      Author: enrico
 */

#ifndef SRC_HASH_MULTISPACEDQMER_H_
#define SRC_HASH_MULTISPACEDQMER_H_

#include <string>
#include <vector>
#include <algorithm>
#include <bitset>
#include <limits>
#include "./SpacedQmer.h"

using namespace std;

// Typedef for a vector of Pos_Ones structures
typedef vector<Pos_Ones> V_Pos_Ones;

// Typedef for a vector of indices referring to the positions inside the seed
typedef vector<size_t> Position;

// PreviousShiftMulti contains all the information needed to retrieve
// some positions from a hash that has already been computed
struct PreviousShiftMulti {
    Position one_to_remove; // Ones of the current seed that do not find a correspondence with the previous hash
    Position one_to_keep; // Ones of the current seed that find a correspondence with the previous hash

    // one_exit represents how much the previous hash has to be shifted to recreate the required overlapping
    // It can be negative, indicating a shift in the opposite direction
    int one_exit = 0;

    // shift_min represents the offset between the current position in the sequence and the one relative 
    // to the previous hash
    int shift_min = 0;

    // seed_num is the seed index in the vector of seeds from which the positions have to be recovered
    int seed_num = 0;

    // mask needed to remove the unwanted values (one_to_remove) from the partial current hash that can be 
    // computed using the previous hash to merge multiple partial hashes
    uint64_t mask = 0;
};

// The current hash is computed starting from a group of previous hashes stored in a vector
typedef vector<PreviousShiftMulti> V_PreviusShiftMulti;

// In the struct groupPrevious, information about all the previous hashes used and which 
// positions must be recomputed are stored
struct groupPrevious {
    V_PreviusShiftMulti prev; // Information about all the previous hashes used

    // not_covered contains the positions that must be computed anew from the sequence 
    // because they cannot be recovered from the previous hashes
    Position not_covered;
};

// V_groupPrevious contains all the information to manage both the standard case (where the 
// "best" group found by our greedy approach is fully available) and the transient case 
// (where the availability of previous hashes is limited)
typedef vector<groupPrevious> V_groupPrevious;

// Struct to store seed information
struct SeedInfo {
    V_groupPrevious group_previous; // Information about all the previous hashes used and positions to be recomputed

    // contains all the positions in the seed where a '1' is present. It is independent of the shift relative to the sequence,
    // so it is stored only once per seed.
    Position pos_ones;
};

// MultiSeedInfo contains all the computed information for each seed, for both standard and transient steps
// and for all the previous hashes for each case [seedindex][standard/transient_step].prev[previous_shift_index]
typedef vector<SeedInfo> MultiSeedInfo;

// Struct passed to the hashing function in case of hashing row ordered
// In this case, additional information about the initial and final transient parts is needed
struct MultiSeedInfoRow {
    MultiSeedInfo info; // Information about all the seeds
    size_t transient1; // Length of the first transient part
    size_t transient2; // Length of the second transient part
};

// Class that computes the information required for efficiently hashing DNA sequences
// by considering multiple spaced seeds simultaneously, all having the same weight and length
// It also returns this information as a reference
class MultiSpacedQmer {
public:
    MultiSpacedQmer(const vector<SpacedQmer> &spaced_qmers); // Constructor

    // Get the number of spaced q-mers
    inline size_t GetLength() const {
        return this->spaced_qmers.size();
    }

    // Get the information required for computing the hash by columns
    inline const MultiSeedInfo& Get_multi_seed_info_col() const {
        return this->multi_seed_info_col;
    }

    // Get the information required for computing the hash by rows
    inline const MultiSeedInfoRow& Get_multi_seed_info_row() const {
        return this->multi_seed_info_row;
    }

    // Reset function to reinitialize the object with new spaced q-mers
    void reset(vector<SpacedQmer> spaced_qmers);

private:
    vector<SpacedQmer> spaced_qmers; // Structure for each of the single seeds
    vector<string> spaced_qs; // Vector of spaced q-mer strings
    vector<Position> v_pos_ones; // Vector of positions of ones
    vector<Position> v_pos_pos_ones; // Vector of positions relative to pos_ones

    MultiSeedInfo multi_seed_info_col; // Final vector with all the needed information for column computation
    MultiSeedInfoRow multi_seed_info_row; // Final vector with all the needed information for row computation

    // Compute information for each seed for every case (standard and transient steps)
    void SetMultiSeedInfoCol();
    // Compute information for position i for each seed
    void ProcessMultiSeedCol(size_t i);

    // Compute information for each seed for every case (standard and transient steps)
    void SetMultiSeedInfoRow();
    // Compute information for position i for each seed
    void ProcessMultiSeedRow(MultiSeedInfo& multi_seed_info_row, int i);
};

// Function to visualize the information about a previous shift
void print_shift_multi(PreviousShiftMulti s);

#endif /* SRC_HASH_MULTISPACEDQMER_H_ */
