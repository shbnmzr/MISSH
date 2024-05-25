/*
 * SpacedQmer.h
 *
 *  Created on: 20/lug/2016
 *      Author: samuele
 *      Modified by: enrico
 */

#ifndef SRC_HASH_SPACEDQMER_H_
#define SRC_HASH_SPACEDQMER_H_

#include <string>
#include <vector>
#include <algorithm>
#include <bitset>
#include <limits>

using namespace std;

// Structure to store information about positions of ones in a spaced seed
struct Pos_Ones {
    size_t index_one = numeric_limits<size_t>::max(); // Index of '1' in the spaced seed
    size_t n_one = 0; // Number of '1's
    size_t pos_start = 0; // Starting position
    size_t n_one_before = 0; // Number of '1's before the first '1' in the unit
};

typedef vector<Pos_Ones> V_Pos_Ones;

// For previous hashes, position is a vector of indices referring to positions inside the seed
typedef vector<size_t> Position;
typedef vector<Position> V_Position;

// PreviousShift structure contains all the information regarding a previous shift of the spaced seed
// It is used to calculate the speedup previous more efficiently
struct PreviousShift {
    Position one_to_change; // Positions to change
    Position one_to_remove; // Positions to remove
    Position one_to_keep; // Positions to keep

    // Index of the first overlapping '1' in the spaced seed
    size_t one_exit = 0;

    // Shift of the current PreviousShift relative to the "main" one
    size_t shift_min = 0;

    // Mask required for calculation using ISSH
    uint64_t mask = 0;

    // Get the total size for the current shift
    inline size_t GetSize() const {
        return this->one_to_change.size() + this->one_to_remove.size() + this->one_exit + (this->one_to_remove.empty() ? 0 : 2);
    }
};

// Vector of PreviousShift structures associated with a spaced q-mer
typedef vector<PreviousShift> V_PreviusShift;
typedef vector<V_PreviusShift> V_V_PreviusShift;

class SpacedQmer {
public:
    SpacedQmer(); // Default constructor
    SpacedQmer(string spaced_qmer, size_t numprev); // Parameterized constructor

    // Get the weight (number of '1's) of the spaced q-mer
    inline size_t GetWeight() const {
        return this->pos_one.size();
    }

    // Get the length of the spaced q-mer
    inline size_t GetQ() const {
        return this->spaced_q.length();
    }

    // Check if a character at a given index is '1'
    inline bool isOne(size_t index) const {
        return this->spaced_q[index] == '1';
    }

    // Get the positions of '1's in the spaced q-mer
    inline const Position& GetPosOne() const {
        return this->pos_one;
    }

    // Get the shifts for a single previous hash
    inline const V_PreviusShift& GetShiftMinChange() const {
        return this->shift_min_change;
    }

    // Get the shifts for multiple previous hashes
    inline const V_V_PreviusShift& GetMultipleShifts() const {
        return this->multiple_shifts;
    }

    // Get the pointer to the shifts for multiple previous hashes
    inline const V_V_PreviusShift* GetMultipleShiftsPointer() const {
        return &this->multiple_shifts;
    }

    // Convert the spaced q-mer to string
    inline const string& toString() const {
        return this->spaced_q;
    }

    // Reset function to reinitialize the object
    void reset(string spaced_qmer, size_t numprev);

private:
    string spaced_q; // The actual string of ones and zeros, the original spaced seed
    Position pos_one; // Vector of indices corresponding to ones in the seed
    Position pos_pos_one; // Vector of positions in the pos_one vector where correspondence is yet to be recovered
    V_PreviusShift shift_min_change; // Data of all possible overlapping shifted seeds
    V_V_PreviusShift multiple_shifts; // Vector of groups of previous hashes to reuse, for both steady-state and transient phases

    size_t num_prev = 0; // Number of previous hashes to use, set from terminal

    void SaveIndexOne(); // Save the indices of '1's in the spaced q-mer
    void GetShiftMax(V_PreviusShift& shift_max); // Calculate maximum shifts for previous hashes
    void SetMultipleShifts(size_t index); // Calculate shifts for a specific previous hash group
    void SetAllMultipleShift(); // Calculate shifts for multiple previous hashes
    void SetBitMasks(); // Create bit masks for removing unnecessary values from previous hashes
};

// Support functions to print some results
void printp(Position p); // Print the elements of a Position vector
void print_shift(PreviousShift s); // Print the details of a PreviousShift object

// Support functions to manage the vector of positions that have not yet been recovered
void deleteElement(Position& pos, size_t index); // Delete positions identified as recoverable from the vector
bool isContained(Position pointer, Position pos, size_t index); // Check if a position is already recovered

#endif /* SRC_HASH_SPACEDQMER_H_ */
