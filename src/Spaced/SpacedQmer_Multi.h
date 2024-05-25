#include "SpacedQmer.h"
#include "../Hash/HashType.h"

// Structure for multi-unit processing
struct MapUnit {
    vector<size_t> n_one; // Number of ones in the seeds
    vector<V_Pos_Ones> v_v_pos; // Positions of ones in the seeds
};

// Structure for storing extended previous shift information
struct PreviusShift_Ext : public PreviousShift {
    size_t prev_qmer = numeric_limits<size_t>::max(); // Previous q-mer index
    const SpacedQmer* current_sp_ptr; // Pointer to the current spaced q-mer
    const SpacedQmer* prev_sp_ptr; // Pointer to the previous spaced q-mer

    // Check if the previous spaced q-mer is correctly set
    bool isCorrectSpacedPrevious() const {
        return this->prev_qmer != numeric_limits<size_t>::max();
    }
};
typedef vector<PreviusShift_Ext> PreviusShift_Ext_V;

// Function to calculate shift changes between spaced q-mers
inline static void GetTableShiftFirstOnSecond(const vector<SpacedQmer>& v_spaced, vector<vector<PreviusShift_Ext_V>>& table_shift_change) {
    table_shift_change.clear(); 
    table_shift_change.resize(v_spaced.size());
    for(size_t i = 0; i < table_shift_change.size(); ++i)
        table_shift_change[i].resize(v_spaced.size());

    // Calculate positions to save with various shifts and then save
    for(size_t ss = 0; ss < v_spaced.size(); ++ss) {
        for(size_t s = 0; s < v_spaced.size(); ++s) { // Calculate shifts for all seeds
            table_shift_change[ss][s].resize(v_spaced[s].toString().size()); // Resize to match length of seed s
            for(size_t i = 0; i < table_shift_change[ss][s].size(); ++i)
                table_shift_change[ss][s][i].prev_qmer = s; // Reference q-mer for hash calculation of spaced ss

            size_t init = 0;
            for(size_t i = 0; i < v_spaced[s].toString().size(); ++i) { // For all possible shifts
                bool find = false;
                for(size_t j = init; j < v_spaced[s].GetPosOne().size(); ++j) {
                    if(v_spaced[s].GetPosOne()[j] >= i) {
                        init = j;
                        find = true;
                        break;
                    }
                }
                if(!find)
                    init = v_spaced[s].GetPosOne().size(); // Skip next cycle without further checks

                for(size_t j = init; j < v_spaced[s].GetPosOne().size() && (j-init) < v_spaced[ss].GetPosOne().size(); ++j) {
                    if(v_spaced[ss].GetPosOne()[j-init] != v_spaced[s].GetPosOne()[j]-i) {
                        table_shift_change[ss][s][i].one_to_remove.push_back(j-init);
                        table_shift_change[ss][s][i].one_to_change.push_back(j-init);
                    }
                }

                table_shift_change[ss][s][i].one_exit = init; // Shift of ones to perform on s
            }
        }
    }
}

// Function to get the minimum shift changes
inline static void GetShiftMinChange(const vector<SpacedQmer>& v_spaced, vector<PreviusShift_Ext_V>& v_shift_min) {
    // Initialize matrix to save all possible shifts between all pairs
    vector<vector<PreviusShift_Ext_V>> shift_firstid_on_secondid;
    GetTableShiftFirstOnSecond(v_spaced, shift_firstid_on_secondid);

    size_t max_size = 0;
    for(size_t ss = 0; ss < v_spaced.size(); ++ss)
        for(size_t s = 0; s < v_spaced.size(); ++s)
            if(max_size < shift_firstid_on_secondid[ss][s].size())
                max_size = shift_firstid_on_secondid[ss][s].size();

    // Choose which shifts to perform at each step, considering removals and insertions
    v_shift_min.clear(); 
    v_shift_min.resize(v_spaced.size());
    for(size_t ss = 0; ss < v_spaced.size(); ++ss) {
        for(size_t i = 0; i < max_size; ++i) {
            size_t limit_evaluation_previous = i == 0 ? ss : v_spaced.size();
            for(size_t s = 0; s < limit_evaluation_previous; ++s) {
                if(i < shift_firstid_on_secondid[ss][s].size()) {
                    PreviusShift_Ext current_shift_min = shift_firstid_on_secondid[ss][s][i];
                    size_t current_size = current_shift_min.GetSize();
                    current_shift_min.shift_min = i; // Guess

                    if(v_shift_min[ss].size() <= i) { // Nothing saved yet
                        if(i > 0) { // Evaluate first seed
                            const PreviusShift_Ext& saved_prev_shift_min = v_shift_min[ss][i-1];
                            size_t saved_size = saved_prev_shift_min.GetSize();

                            bool size_saved_better = saved_size < current_size;
                            if(size_saved_better && saved_prev_shift_min.isCorrectSpacedPrevious()) {
                                current_shift_min = saved_prev_shift_min;
                                current_size = current_shift_min.GetSize();
                            }
                        }
                        if(!(i == 0 && current_size > v_spaced[ss].GetWeight()))
                            v_shift_min[ss].push_back(current_shift_min);
                    } else {
                        const PreviusShift_Ext& saved_shift_min = v_shift_min[ss][i];
                        size_t saved_size = saved_shift_min.GetSize();

                        bool size_saved_better = saved_size < current_size; // Prefer new current shift
                        if(!size_saved_better || !saved_shift_min.isCorrectSpacedPrevious())
                            v_shift_min[ss][i] = current_shift_min;
                    }
                }
            }
            if(i == 0 && v_shift_min[ss].empty())
                v_shift_min[ss].push_back(PreviusShift_Ext());
        }
    }

    // Set pointers
    for(size_t ss = 0; ss < v_shift_min.size(); ++ss) {
        for(size_t i = 0; i < v_shift_min[ss].size(); ++i) {
            v_shift_min[ss][i].current_sp_ptr = &v_spaced[ss];
            if(v_shift_min[ss][i].isCorrectSpacedPrevious())
                v_shift_min[ss][i].prev_sp_ptr = &v_spaced[v_shift_min[ss][i].prev_qmer];
        }
    }
}

// Function to get rotated minimum shift changes
inline static void GetShiftMinChangeRotated(const vector<SpacedQmer>& v_spaced, vector<PreviusShift_Ext_V>& v_shift_min) {
    // Initialize matrix to save all possible shifts between all pairs
    vector<PreviusShift_Ext_V> shift_firstid_on_secondid;
    GetShiftMinChange(v_spaced, shift_firstid_on_secondid);
    size_t size_max = max_element(shift_firstid_on_secondid.begin(), shift_firstid_on_secondid.end(), [](PreviusShift_Ext_V& a, PreviusShift_Ext_V& b) { return a.size() < b.size(); })->size();
    v_shift_min.clear();
    v_shift_min.resize(size_max, PreviusShift_Ext_V(v_spaced.size()));
    for(size_t ss = 0; ss < shift_firstid_on_secondid.size(); ++ss) {
        for(size_t i = 0; i < size_max; ++i) {
            if(i < shift_firstid_on_secondid[ss].size()) // Element exists, copy it
                v_shift_min[i][ss] = shift_firstid_on_secondid[ss][i];
            else
                v_shift_min[i][ss] = shift_firstid_on_secondid[ss].back();
        }
    }
}

class SpacedQmer_Multi {
public:
    inline SpacedQmer& operator[](size_t i) { return this->v_spaced[i]; }
    inline const SpacedQmer& operator[](size_t i) const { return this->v_spaced[i]; }
    inline size_t size() const { return this->v_spaced.size(); }
    inline const MapUnit& getMapUnit() const { return map_unit; }
    inline const vector<PreviusShift_Ext_V>& getShiftMin() const { return v_shift_min; }
    inline const vector<PreviusShift_Ext_V>& getShiftMinRotated() const { return v_shift_min_rotated; }

    inline void init(const vector<SpacedQmer>& v_spaced) {
        this->v_spaced = v_spaced;
        GetShiftMinChange(this->v_spaced, this->v_shift_min);
        GetShiftMinChangeRotated(this->v_spaced, this->v_shift_min_rotated);
    }

private:
    vector<SpacedQmer> v_spaced;
    MapUnit map_unit;
    vector<PreviusShift_Ext_V> v_shift_min;
    vector<PreviusShift_Ext_V> v_shift_min_rotated;
};
