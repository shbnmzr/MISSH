#include "./MultiSpacedQmer.h"

#include <iostream>
#include <cmath>

/**
 * Constructor
 * 
 * It calls the functions to generate the structures containing data that will be used during the hashing computation ISSH column and row 
 * 
 * @param spaced_qmers Array containing the information about each single spaced seed
 */
MultiSpacedQmer::MultiSpacedQmer(const vector<SpacedQmer> &spaced_qmers) 
{
    // Initialize local information about all the seeds given in input, locally create copies of v_pos_ones and v_pos_pos_ones for each seed
    this->spaced_qmers = spaced_qmers;
    spaced_qs.resize(GetLength());
    v_pos_ones.resize(GetLength());
    v_pos_pos_ones.resize(GetLength());
    size_t length = spaced_qmers[0].GetQ();
    size_t weight = spaced_qmers[0].GetWeight();

    // Check if all seeds have the same length and weight, and store their information
    for(size_t i = 0; i < GetLength(); i++)
    {
        if(length != spaced_qmers[i].GetQ() || weight != spaced_qmers[i].GetWeight())
        {
            cerr << endl << "The seed in the group does not have the same length or weight\n" << flush;
            exit(1);
        }

        spaced_qs[i] = spaced_qmers[i].toString();
        size_t k = 0;
        for(size_t j = 0; j < this->spaced_qs[i].length(); ++j)
            if(spaced_qmers[i].isOne(j))
            {
                v_pos_ones[i].push_back(j);
                v_pos_pos_ones[i].push_back(k);
                k++;
            }    
    }
    // Populate structure for row computation
    this->SetMultiSeedInfoRow();    
    
    // Locally create new copies of v_pos_pos_ones for each seed
    this->v_pos_pos_ones.clear(); this->v_pos_pos_ones.shrink_to_fit();
    v_pos_pos_ones.resize(GetLength());
    for(size_t i = 0; i < GetLength(); i++)
    {
        size_t k = 0;
        for(size_t j = 0; j < this->spaced_qs[i].length(); ++j)
            if(spaced_qmers[i].isOne(j))
            {
                v_pos_pos_ones[i].push_back(k);
                k++;
            }    
    }
    // Populate structure for column computation
    this->SetMultiSeedInfoCol();
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
// ISSH multi seed preprocessing
//////////////////////////////////////////////////////////////////////////////////

/**
 * Compute information needed to perform hashing by columns for both transient and standard computation
 * It populates a structure having as many items as the number of seeds, for each seed it is needed a "recipe" for the standard computation,
 * that is when all the best group of hashes have already been computed, and for each of the transient steps.
 * Each of these "recipes" is composed of a group of information about hashes that have already been computed and how to use them 
 * (performing shifts and masking unwanted bits), in order to recover from them positions that are useful for the current computation.
 */
void MultiSpacedQmer::SetMultiSeedInfoCol()
{    
    // Vector as long as the number of seeds
    this->multi_seed_info_col.resize(GetLength());
    
    // For each seed make room for information about standard computation
    for(size_t j = 0; j < GetLength(); j++)
    {
        multi_seed_info_col[j].group_previous.resize(1);
        multi_seed_info_col[j].pos_ones = v_pos_ones[j];
    }
    // Compute information about standard computation
    this->ProcessMultiSeedCol(0);
    
    // Once we know how far back is the furthest possible previous hash that is used 
    // in the standard computation, we know how long the transient will be
    int furthest_pos = 0;
    for(size_t j = 0; j < GetLength(); j++)
    {
        for(size_t i = 0; i < this->multi_seed_info_col[j].group_previous[0].prev.size(); i++)
        {
            furthest_pos = furthest_pos < this->multi_seed_info_col[j].group_previous[0].prev[i].shift_min ? this->multi_seed_info_col[j].group_previous[0].prev[i].shift_min : furthest_pos;
        }
    }
    furthest_pos++; 

    // Resize the info for each seed in order to store all the transient steps
    for(size_t j = 0; j < GetLength(); j++)
    {
        multi_seed_info_col[j].group_previous.resize(furthest_pos);
    }
    // Compute info for the transient steps
    for(int k = 1; k < furthest_pos; k++)
    {
        this->ProcessMultiSeedCol(k);
    }
}

/**
 * Function that computes the information needed to perform the hashing of a sequence by columns
 * @param index indicates if the computation is about the standard case (0) or a certain step of the transient (1...n)
 */ 
void MultiSpacedQmer::ProcessMultiSeedCol(size_t index)
{
    // Cycle for all the seeds
    for(size_t y = 0; y < GetLength(); y++)
    {    
        // Initialize the vector that will contain all the shifts identified for this particular case and store it inside 
        // the main structure in correspondence with the right seed and right index (standard or transient step)
        V_PreviusShiftMulti current_shift_group;
        this->multi_seed_info_col[y].group_previous[index].prev = current_shift_group;
        
        // Positions that are not recovered yet by using previous hashes, initially contains all the positions
        Position pos_not_covered_yet = v_pos_pos_ones[y];

        // y -> current seed
        // n -> seed from which I want to recover positions
        
        // done is set to true when no more positions can be recovered from any of the available hashes
        bool done = false;
        while(!done)
        {    
            // Initialize a current shift that is empty at the beginning
            PreviousShiftMulti curr_best;
        
            // max_num_shifts express how far from the current hash the previous hashes can be considered. It is unlimited for
            // the standard computation and limited in case of a transient step
            size_t max_num_shifts = index != 0 ? index - 1 : this->spaced_qs[index].size() - 1;

            // Cycle for each seed from which we will recover positions
            for(size_t n = 0; n < GetLength(); n++)
            {
                // It is possible to recover positions from hash computed in the same point of the sequence if those 
                // have already been computed
                size_t startingPoint = n < y ? 0 : 1;

                // Cycle for all the possible shifts between the two seeds
                for(size_t i = startingPoint; i <= max_num_shifts; ++i)
                { 
                    // init represents the number of positions that are discarded at the beginning of the hash that 
                    // is used to recover from. It can also be negative, meaning that the positions are discarded 
                    // at the beginning of the current hash. This represents the third degree of freedom. 
                    // First degree: seed from which we recover from
                    // Second degree: relative distance on the DNA sequence between the current hash and the one that we want to use 
                    // Third degree: offset relative to the overlapping between the two hashes.

                    // In order to be sure not to leave behind good possible overlappings, we, for now, start from -weight and
                    // arrive up to weight, knowing that the weight (number of "1"s inside the seed) corresponds to the q-gram length.
                    // Surely, a tighter bound can be proven, but for now it is good enough.
                    // In case that we are recovering from hashes computed thanks to the same seed there's no point in using negative init.
                    int startingInit = n == y ? 1 : -(this->spaced_qmers[y].GetWeight()-1);
                    int num_one_before_shift = (this->spaced_qmers[y].GetWeight()-1);
                    
                    for(int init = startingInit; init <= num_one_before_shift; ++init)
                    {
                        // Initialize temp shift about the current situation
                        PreviousShiftMulti temp;
                        // The current seed is scanned by the index (j - init) that is equal to 0 at the beginning and W-init at the end 
                        for(int j = 0; j < (int)this->v_pos_ones[n].size(); ++j) // j index on the old seed
                        {
                            // If the index scanning the current seed is still inside the hash, otherwise it is useless to start the comparison
                            if((j-init) < (int)this->v_pos_ones[y].size() && (j-init) >= 0) 
                            {
                                // If the index of the current '1' on the previous seed 'n' is useless because it will never overlap 
                                // || there is no overlapping between previous and current seed
                                // || current position is already recovered
                                if(v_pos_ones[n][j] < i || this->v_pos_ones[y][j-init] != this->v_pos_ones[n][j]-i || !isContained(pos_not_covered_yet, v_pos_ones[y], v_pos_ones[y][j-init]))
                                {
                                    // Info used to compute the masks in order to clear unwanted values
                                    temp.one_to_remove.push_back(j-init);
                                }
                                else
                                {
                                    // Info used to compare different previous shifts
                                    temp.one_to_keep.push_back(j-init);
                                }
                            }
                        }
                        // Update info about the computed situation inside temp
                        temp.one_exit = init;
                        temp.shift_min = i;
                        temp.seed_num = n;
                        
                        // Check in a greedy way if the temp is better than the one that was found before.
                        // If the number of recovered positions is the same, but it is closer to the current hash, it will be considered better to shorten the transient
                        if(temp.one_to_keep.size() > curr_best.one_to_keep.size() || (temp.one_to_keep.size() == curr_best.one_to_keep.size() && temp.shift_min < curr_best.shift_min))
                        {
                            curr_best = temp;
                        }
                    }
                }
            }
            // After checking all the possibilities, we find the best previous hash that allows recovering the largest number of positions.
            // These positions must be removed from the vector of positions that are not covered yet.
            for(size_t j = 0; j < curr_best.one_to_keep.size(); j++)
            {
                deleteElement(pos_not_covered_yet, curr_best.one_to_keep[j]);
            }
            // If the current best shift can recover at least one position, add it to the vector of shifts
            // that is used to compute, for this particular seed, for the hash relative to situation "index"
            // If not, the computation can no longer improve on the number of recovered positions.
            if(curr_best.one_to_keep.size() != 0)
            {
                multi_seed_info_col[y].group_previous[index].prev.push_back(curr_best);
            }
            else
            {
                done = true;
            }
        }
        // In order to know which positions are still not recovered and must be recomputed anew
        // store the remaining pos_not_covered_yet inside not_covered.
        this->multi_seed_info_col[y].group_previous[index].not_covered = pos_not_covered_yet;
    }    
    // Compute for all the used previous hashes the relative mask that has the ones on the positions that we want to keep and zeros elsewhere
    // For each seed
    for(size_t j = 0; j < this->multi_seed_info_col.size(); j++)
    {    
        // For each shift in the group
        for (size_t k = 0; k < this->multi_seed_info_col[j].group_previous[index].prev.size(); k++)
        {
            // For each position to keep set the relative position in the mask to 11 (3 in binary) 
            for (size_t i = 0; i < this->multi_seed_info_col[j].group_previous[index].prev[k].one_to_keep.size(); i++)
            {
                this->multi_seed_info_col[j].group_previous[index].prev[k].mask |= (uint64_t)3 << (this->multi_seed_info_col[j].group_previous[index].prev[k].one_to_keep[i] * 2);
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

// Compute information needed to perform hashing by rows
void MultiSpacedQmer::SetMultiSeedInfoRow()
{    
    // Vector as long as the number of seeds
    this->multi_seed_info_row.info.resize(GetLength());

    // Initialize the limitations for the standard case (actually no limitation)
    multi_seed_info_row.transient1 = this->spaced_qs[0].length()-1;
    multi_seed_info_row.transient2 = this->spaced_qs[0].length()-1;
    
    // Compute information about standard computation
    this->ProcessMultiSeedRow(multi_seed_info_row.info, 0);

    // Once we know the dependencies for the standard computation, we know how long
    // the two transients will be
    int furthest_pos_back = 0;
    int furthest_pos_front = 0;
    for(size_t j = 0; j < GetLength(); j++)
    {
        for(size_t i = 0; i < multi_seed_info_row.info[j].group_previous[0].prev.size(); i++)
        {
            furthest_pos_back = furthest_pos_back < multi_seed_info_row.info[j].group_previous[0].prev[i].shift_min ? multi_seed_info_row.info[j].group_previous[0].prev[i].shift_min : furthest_pos_back;
            furthest_pos_front = furthest_pos_front > multi_seed_info_row.info[j].group_previous[0].prev[i].shift_min ? multi_seed_info_row.info[j].group_previous[0].prev[i].shift_min : furthest_pos_front;
        }
        multi_seed_info_row.info[j].pos_ones = v_pos_ones[j];
    }
    
    multi_seed_info_row.transient1 = furthest_pos_back;
    
    multi_seed_info_row.transient2 = -furthest_pos_front;
    
    // Compute info for the transient steps
    for(int k = 1; k < furthest_pos_back+1; k++)
    {
        this->ProcessMultiSeedRow(multi_seed_info_row.info, k);
    }

    for(int k = furthest_pos_front; k < 0; k++)
    {
        ProcessMultiSeedRow(multi_seed_info_row.info, k);
    }
}

/**
 * Function that computes the information needed to perform the hashing of a sequence by rows.
 * At each call a step of the transient or the standard computation is addressed
 * @param step indicates if the computation is about the standard case (0), a certain step of the first transient (1...n) or the second transient (-1...-m)
 * @param multi_seed_info_row structure where to store the computed information
 */ 
void MultiSpacedQmer::ProcessMultiSeedRow(MultiSeedInfo& multi_seed_info_row, int step)
{
    // Cycle on all the seeds
    int index = 0;
    for(size_t y = 0; y < GetLength(); y++)
    {    
        groupPrevious current_shift_group;
        
        V_PreviusShiftMulti current_shift_prev;
        
        current_shift_group.prev = current_shift_prev;

        multi_seed_info_row[y].group_previous.push_back(current_shift_group);
        index = multi_seed_info_row[y].group_previous.size() - 1;
        Position pos_not_covered_yet = v_pos_pos_ones[y];
    
        // y -> current seed
        // n -> seed from which I want to recover positions
        bool done = false;
        while(!done)
        {    
            // Initialize a current shift that is empty at the beginning
            PreviousShiftMulti curr_best;
            // max_num_shifts express how far from the current hash the previous hashes can be considered. It is unlimited for
            // the standard computation and limited in case of a transient step
            int max_num_shifts = step > 0 ? step - 1 : this->multi_seed_info_row.transient1;
            
            // min_num_shifts represents the number of shifts that can be considered in the future 
            int min_num_shifts = step >= 0 ? -this->multi_seed_info_row.transient2 : step + 1; 
            
            // Consider only seeds that are ordered before the current one
            for(size_t n = 0; n <= y; n++)
            {
                // The starting point represents the offset from the current hash to the one that is considered for 
                // recovering the position. Since the previous row is always fully available, there are limitations only if
                // the row that is considered is the current one. 
                int startingPoint = n < y ? min_num_shifts : 1;
                
                for(int i = startingPoint; i <= max_num_shifts; ++i)
                { 
                    // At this point, the current hash and the previous hash are fixed, but we need to explore all the 
                    // possible overlapping between the two.

                    int num_one_before_shift = 0;
                    int num_one_after_shift = 0;
                    // The behavior is different if the shift is to the left or to the right; if it is positive, it is to the left
                    // If i >= 0, it is the same as the case for the column computation
                    if(i >= 0)
                    {
                        int startingInit = n == y ? 1 : -(this->spaced_qmers[y].GetWeight()-1);
                        int num_one_before_shift = (this->spaced_qmers[y].GetWeight()-1);

                        for(int init = startingInit; init <= num_one_before_shift; ++init)
                        {
                            // Initialize temp shift about the current situation
                            PreviousShiftMulti temp;
                            // The current seed is scanned by the index (j - init) that is equal to 0 at the beginning and W-init at the end 
                            for(size_t j = 0; j < this->v_pos_ones[n].size(); ++j) // j index on the old seed
                            {
                                if((j-init) < this->v_pos_ones[y].size() && (j-init) >= 0)
                                {
                                    if((int)v_pos_ones[n][j] < i || this->v_pos_ones[y][j-init] != this->v_pos_ones[n][j]-i || !isContained(pos_not_covered_yet, v_pos_ones[y], v_pos_ones[y][j-init]))
                                    {
                                        temp.one_to_remove.push_back(j-init);
                                    }
                                    else
                                    {
                                        temp.one_to_keep.push_back(j-init);
                                    }
                                }
                            }
                            
                            temp.one_exit = init;
                            temp.shift_min = i;
                            temp.seed_num = n;
                            if(temp.one_to_keep.size() > curr_best.one_to_keep.size() || (temp.one_to_keep.size() == curr_best.one_to_keep.size() && abs(temp.shift_min) < abs(curr_best.shift_min)))
                            {
                                curr_best = temp;
                            }
                        }                    
                    }
                    // Manage case in which the previously computed hash is on the right with respect to 
                    // the current hash
                    else
                    {
                        num_one_after_shift = 0;
                        for(; (int)this->v_pos_ones[y][num_one_after_shift] < -i; num_one_after_shift++);

                        num_one_before_shift = v_pos_ones[y].size()-1;
                        for(; (int)this->v_pos_ones[y][num_one_before_shift] > -i; num_one_before_shift--);
                        num_one_before_shift = v_pos_ones[y].size() - num_one_before_shift;
                        num_one_after_shift = -num_one_after_shift;

                        // Even the starting init can be modified to make fewer computations.
                        // init is an index in the vector of positions that are one. It represents, given the shift,
                        // how many ones at the beginning of the old hash we ignore
                        int startingInit = n == y ? 1 : num_one_after_shift;
                        for(int init = startingInit; init <= num_one_before_shift; ++init)
                        {
                            PreviousShiftMulti temp;
                            for(int j = 0; j < (int)this->v_pos_ones[n].size(); ++j)
                            {
                                if((j-init) < (int)this->v_pos_ones[y].size() && (j-init) >= 0)
                                {
                                    if(this->v_pos_ones[y][j-init] != this->v_pos_ones[n][j]-i || !isContained(pos_not_covered_yet, v_pos_ones[y], v_pos_ones[y][j-init]))
                                    {
                                        temp.one_to_remove.push_back(j-init);
                                    }
                                    else
                                    {
                                        temp.one_to_keep.push_back(j-init);
                                    }
                                }
                            }
                            temp.one_exit = init;
                            temp.shift_min = i;
                            temp.seed_num = n;
                            if(temp.one_to_keep.size() > curr_best.one_to_keep.size() || (temp.one_to_keep.size() == curr_best.one_to_keep.size() && abs(temp.shift_min) < abs(curr_best.shift_min)))
                            {
                                curr_best = temp;
                            }
                        }
                    }
                }
            }
            // Update list of positions that have not been recovered yet
            for(size_t j = 0; j < curr_best.one_to_keep.size(); j++)
            {
                deleteElement(pos_not_covered_yet, curr_best.one_to_keep[j]);
            }
            // Add computed information to the structure
            if(curr_best.one_to_keep.size() != 0)
            {
                multi_seed_info_row[y].group_previous[index].prev.push_back(curr_best);
            }
            else
            {
                done = true;
            }
        }
        multi_seed_info_row[y].group_previous[index].not_covered = pos_not_covered_yet;
    }
    // Compute for all the used previous hashes the relative mask that has the ones on the positions that we want to keep and zeros elsewhere
    // For each seed
    for(size_t j = 0; j < multi_seed_info_row.size(); j++)
    {    
        // For each shift in the group
        for (size_t k = 0; k < multi_seed_info_row[j].group_previous[index].prev.size(); k++)
        {
            // For each position to keep set the relative position in the mask to 11 (3 in binary) 
            for (size_t i = 0; i < multi_seed_info_row[j].group_previous[index].prev[k].one_to_keep.size(); i++)
            {
                multi_seed_info_row[j].group_previous[index].prev[k].mask |= (uint64_t)3 << (multi_seed_info_row[j].group_previous[index].prev[k].one_to_keep[i] * 2);
            }
        }
    }
}

/**
 * Print all the information about a single previous hash used for recovering positions. Used mostly for debugging purposes.
 * @param s structure that contains the information about a previously computed hash
 */ 
void print_shift_multi(PreviousShiftMulti s)
{
    cout << "\nC_" << s.shift_min - s.one_exit << "," << s.shift_min;
    //cout << "\none_to change= ";
    //printp(s.one_to_change);
    // Vector containing index at which value the hash has to be removed
    cout << "\none_to_remove= ";
    printp(s.one_to_remove);

    cout << "\none_to_keep= ";
    printp(s.one_to_keep);
    cout << "\none_exit= " << s.one_exit;
    cout << "\nshift_min= " << s.shift_min;
    cout << "\nseed_num= " << s.seed_num;
    if(s.mask != 0)
        cout << "\nMask: " << bitset<42>(s.mask);
    cout << endl << endl;
}