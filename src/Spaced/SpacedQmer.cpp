/*
 * SpacedQmer.cpp
 *
 */

#include "SpacedQmer.h"

#include <iostream>

// Default constructor
SpacedQmer::SpacedQmer() {
    // Initialize with an empty spaced q-mer and num_prev as 0
    this->reset("", 0);
}

// Parameterized constructor
SpacedQmer::SpacedQmer(string spaced_qmer, size_t numprev) {
    // Initialize with the given spaced q-mer and num_prev
    this->reset(spaced_qmer, numprev);
}

// Reset function to reinitialize the object
void SpacedQmer::reset(string spaced_qmer, size_t numprev) {
    this->num_prev = numprev; // Set num_prev
    this->spaced_q = spaced_qmer; // Set spaced q-mer string
    this->SaveIndexOne(); // Save the indices of '1's in the spaced q-mer
    // Calculate the shifts for a single previous hash
    this->GetShiftMax(this->shift_min_change);
    // Calculate the groups of previous shifts for hash computation
    this->SetAllMultipleShift();
}

// Function to save the indices of '1's in the spaced q-mer
void SpacedQmer::SaveIndexOne() {
    this->pos_one.clear();
    this->pos_one.shrink_to_fit();
    size_t k = 0;
    for (size_t i = 0; i < this->spaced_q.length(); ++i)
        if (this->isOne(i)) {
            this->pos_one.push_back(i);
            this->pos_pos_one.push_back(k);
            k++;
        }
}

// FSH preprocessing: Calculate maximum shifts for previous hashes
void SpacedQmer::GetShiftMax(V_PreviusShift& shift_max) {
    shift_max.clear();
    shift_max.resize(this->spaced_q.size());
    shift_max.shrink_to_fit();
    size_t init = 0;

    bool find;
    for (size_t i = 1; i < this->spaced_q.size(); ++i) { // For all possible shifts except the first one
        find = false;
        for (size_t j = init; j < this->pos_one.size(); ++j) {
            if (this->pos_one[j] >= i) {
                init = j;
                find = true;
                break;
            }
        }
        if (!find) {
            init = this->pos_one.size(); // Skip next cycle without further checks
        }
        for (size_t j = init; j < this->pos_one.size(); ++j) { // For all positions of the second shifted vector
            if (this->pos_one[j-init] != this->pos_one[j]-i) {
                shift_max[i].one_to_remove.push_back(j-init);
                shift_max[i].one_to_change.push_back(j-init);
            } else {
                shift_max[i].one_to_keep.push_back(j-init);
            }
        }
        shift_max[i].one_exit = init;
        shift_max[i].shift_min = i; // Guess

        const PreviousShift& prev_shift_min = shift_max[shift_max[i-1].shift_min];
        size_t size_previus = prev_shift_min.GetSize();
        size_t size_current = shift_max[i].GetSize();
        if (i > 1 && size_previus < size_current) {
            shift_max[i] = shift_max[i-1];
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

// ISSH preprocessing

// Function to set all multiple shifts
void SpacedQmer::SetAllMultipleShift() {
    // Calculate the best previous shift group for steady-state
    this->multiple_shifts.resize(1);
    this->SetMultipleShifts(0);

    // Calculate the length of the transient phase: find the furthest previous shift
    size_t furthest_pos = this->multiple_shifts[0][0].shift_min;
    for (size_t i = 1; i < this->multiple_shifts[0].size(); i++) {
        if (furthest_pos < this->multiple_shifts[0][i].shift_min)
            furthest_pos = this->multiple_shifts[0][i].shift_min;
    }

    // Extend the vector length to include both transient and steady-state shifts
    this->multiple_shifts.resize(furthest_pos);

    // Calculate all shift groups for each step of the transient phase
    for (size_t k = 1; k < furthest_pos; k++) {
        this->SetMultipleShifts(k);
    }

    // Calculate the bit masks to remove unnecessary values from previous hashes
    this->SetBitMasks();

    // Debugging prints to understand hash calculation behavior
    if (false) {
        cout << "Seed: " << this->spaced_q << endl << "ha un transitorio lungo " << multiple_shifts.size() << endl;
        for (size_t i = 1; i < multiple_shifts.size(); i++) {
            if (multiple_shifts[i].size() == 0)
                cout << "\nGruppo di shift " << i << " ricalcola tutte le posizioni" << endl;
            else {
                cout << "\nGruppo di shift " << i << " utilizza " << multiple_shifts[i].size() << " hash\n";
                cout << "Numero posizioni da calcolare: " << multiple_shifts[i][0].one_to_change.size() << endl;
                // Uncomment to see details of individual hashes
                //for(size_t j=0; j<multiple_shifts[i].size(); j++)
                //{
                //    print_shift(multiple_shifts[i][j]);
                //}
            }
        }
        cout << "\nGruppo di shift a regime utilizza " << multiple_shifts[0].size() << " hash\n";
        cout << "Numero posizioni da calcolare: " << multiple_shifts[0][0].one_to_change.size() << endl << endl << endl;
        // Uncomment to see details of individual hashes
        //for(size_t j=0; j<multiple_shifts[0].size(); j++)
        //{
        //    print_shift(multiple_shifts[0][j]);
        //}
    }
}

// Function to set multiple shifts for a given index
void SpacedQmer::SetMultipleShifts(size_t index) {
    // Prepare the vector of multiple_shifts[index] containing information about previous hashes to reuse
    this->multiple_shifts[index].clear();
    this->multiple_shifts[index].shrink_to_fit();

    // Vector containing positions in the current group that couldn't be recovered
    Position pos_not_covered_yet = pos_pos_one;

    // Determine the maximum number of previous shifts to use
    size_t num_max_previous;
    if (num_prev == 0)
        num_max_previous = (index > this->spaced_q.size() || index == 0) ? this->spaced_q.size() : index;
    else
        num_max_previous = (index > this->spaced_q.size() || index == 0) ? num_prev : (index > num_prev ? num_prev : index);

    // Continue to recover positions from previous hashes until reaching the maximum or no positions remain
    for (size_t k = 0; k < num_max_previous && pos_not_covered_yet.size() > 1; k++) {
        PreviousShift curr_best;
        size_t max_num_shifts = index != 0 ? index : this->spaced_q.size();

        // Identify the previous shift that recovers the most positions
        for (size_t i = 1; i <= max_num_shifts; ++i) {
            size_t num_one_before_shift = 0;
            for (; pos_one[num_one_before_shift] < i; num_one_before_shift++);

            for (size_t init = 1; init <= num_one_before_shift; ++init) {
                PreviousShift temp;
                for (size_t j = init; j < this->pos_one.size(); ++j) {
                    if ((j - init) < this->pos_one.size()) {
                        if (pos_one[j] < i || this->pos_one[j-init] != this->pos_one[j]-i || !isContained(pos_not_covered_yet, pos_one, pos_one[j-init])) {
                            temp.one_to_remove.push_back(j-init);
                            temp.one_to_change.push_back(j-init);
                        } else {
                            temp.one_to_keep.push_back(j-init);
                        }
                    }
                }
                temp.one_exit = init;
                temp.shift_min = i;
                if (temp.one_to_keep.size() > curr_best.one_to_keep.size()) {
                    curr_best = temp;
                }
            }
        }

        for (size_t j = 0; j < curr_best.one_to_keep.size(); j++) {
            deleteElement(pos_not_covered_yet, curr_best.one_to_keep[j]);
        }

        if (curr_best.one_to_keep.size() != 0)
            multiple_shifts[index].push_back(curr_best);
    }

    if (this->multiple_shifts[index].size() > 0)
        this->multiple_shifts[index][0].one_to_change = pos_not_covered_yet;
}

// Function to create bit masks for removing unnecessary values from previous hashes
void SpacedQmer::SetBitMasks() {
    for (size_t k = 0; k < this->multiple_shifts.size(); k++) {
        for (size_t i = 0; i < this->multiple_shifts[k].size(); i++) {
            for (size_t j = 0; j < this->multiple_shifts[k][i].one_to_remove.size(); ++j) {
                this->multiple_shifts[k][i].mask |= (uint64_t)3 << (this->multiple_shifts[k][i].one_to_remove[j] * 2);
            }
            this->multiple_shifts[k][i].mask = ~this->multiple_shifts[k][i].mask;
        }
    }
}

// Function to delete positions identified as recoverable from the vector of positions yet to be recovered
void deleteElement(Position& pos, size_t index) {
    size_t i = 0;
    for (i = 0; i < pos.size() && pos[i] < index; i++);
    if (pos.size() > 0 && pos[i] == index) {
        for (; i < pos.size(); i++) {
            pos[i] = pos[i + 1];
        }
        pos.pop_back();
    }
}

// Function to check if a position is already recovered
bool isContained(Position pointer, Position pos, size_t index) {
    for (size_t i = 0; i < pointer.size(); i++) {
        if (pos[pointer[i]] == index)
            return true;
    }
    return false;
}

// Function to print the details of a PreviousShift object
void print_shift(PreviousShift s) {
    cout << "\none_to change= ";
    printp(s.one_to_change);
    cout << "\none_to_remove= ";
    printp(s.one_to_remove);
    cout << "\none_to_keep= ";
    printp(s.one_to_keep);
    cout << "\none_exit= " << s.one_exit;
    cout << "\nshift_min= " << s.shift_min;
    cout << "\nsize totale= " << s.GetSize() << endl;
    if (s.mask != 0)
        cout << "Maschera calcolata: " << bitset<42>(s.mask);
    cout << endl << endl;
}

// Function to print the elements of a Position vector
void printp(Position p) {
    for (size_t i = 0; i < p.size(); ++i)
        cout << p[i] << " ";
}
