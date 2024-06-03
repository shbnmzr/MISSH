#ifndef SRC_HASH_SPACEDQMER_H_
#define SRC_HASH_SPACEDQMER_H_

#include <vector>
#include <string>

using namespace std;

typedef vector<size_t> Position;
typedef struct {
    vector<size_t> one_to_remove;
    vector<size_t> one_to_keep;
    vector<size_t> one_to_change;
} PreviousShift;
typedef vector<PreviousShift> V_PreviusShift;
typedef vector<V_PreviusShift> V_V_PreviusShift;

class SpacedQmer {
public:
    SpacedQmer();
    SpacedQmer(string spaced_qmer, size_t numprev);
    void reset(string spaced_qmer, size_t numprev);
    void SaveIndexOne();
    void GetShiftMax(V_PreviusShift& shift_max);
    void SetAllMultipleShift();
    void SetBitMasks();
    
    // New method for set covering
    void SetCoveringLP();

    inline bool isOne(size_t index) const { return this->spaced_q[index] == '1'; }
    inline const Position& GetPosOne() const { return this->pos_one; }
    inline const V_PreviusShift& GetShiftMinChange() const { return this->shift_min_change; }
    inline const V_V_PreviusShift& GetMultipleShifts() const { return this->multiple_shifts; }
    inline const V_V_PreviusShift* GetMultipleShiftsPointer() const { return &this->multiple_shifts; }
    inline const string& toString() const { return this->spaced_q; }

private:
    string spaced_q;
    Position pos_one;
    Position pos_pos_one;
    V_PreviusShift shift_min_change;
    V_V_PreviusShift multiple_shifts;
    size_t num_prev = 0;

    void SetMultipleShifts(size_t index);
};

// Support functions
void printp(Position p);
void print_shift(PreviousShift s);
void deleteElement(Position& pos, size_t index);
bool isContained(Position pointer, Position pos, size_t index);

#endif /* SRC_HASH_SPACEDQMER_H_ */