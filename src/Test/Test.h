#ifndef TEST_TEST_H_
#define TEST_TEST_H_

#include "../Hash/MultiHashFunction.h"
#include "../Input/FileScan.h"
#include "../Chrono/Chrono.h"
#include "../Parameter/Parameter.h"
#include "../Spaced/SpacedQmerLP.h" // Include the SpacedQmerLP header

class Test {
public:
    Test();
    virtual ~Test();
    bool load_sequences(const FileParameter& file);
    bool single_run(const SpacedQmer& spaced, bool test_equals);
    void single_save(const FileParameter& file, std::string dir_output);

    bool multi_run(const std::vector<SpacedQmer>& multi_spaced, bool test_equals);
    void multi_save(const FileParameter& file, const std::vector<SpacedQmer>& multi_spaced, std::string dir_output);

private:
    typedef std::pair<std::vector<SpacedQmer>, std::vector<PreviusShift_Ext_V>> MultiParameter;

    std::vector<std::string> to_hash;
    int64_t read_size_max;
    int64_t read_size_min;
    int64_t read_size_mean;

    Chrono::start_end_time times_naive;
    Chrono::start_end_time times_speedup_previous;
    Chrono::start_end_time times_ISSH;
    Chrono::start_end_time times_linear_programming_set_covering; // New timing for the new method

    Chrono::start_end_time times_multi_naive;
    Chrono::start_end_time times_multi_speedup_previous;
    Chrono::start_end_time times_multi_ISSH_single;    
    Chrono::start_end_time times_multi_speedup_multi_previous;
    Chrono::start_end_time times_multi_ISSH_multi_v1;
    Chrono::start_end_time times_multi_ISSH_multi_col;    
    Chrono::start_end_time times_multi_ISSH_multi_col_parallel;
    Chrono::start_end_time times_multi_ISSH_multi_row;    

    bool single_equals = true;

    SpacedQmer_Multi spaced_qmers;
    
    bool multi_equals = true;

    void single_test_hashes(const SpacedQmer& spaced);
    void single_test_equals(const SpacedQmer& spaced);
    void single_test_naive(const SpacedQmer& spaced);
    void single_test_speedup_previous(const SpacedQmer& spaced);
    void single_test_ISSH(const SpacedQmer& spaced);
    void single_test_linear_programming_set_covering(const SpacedQmer& spaced); // New method

    void multi_test_hashes(const std::vector<SpacedQmer>& multi_spaced);
    void multi_test_equals(const std::vector<SpacedQmer>& multi_spaced);
    void multi_test_naive(const std::vector<SpacedQmer>& multi_spaced);
    void multi_test_speedup_previous(const std::vector<SpacedQmer>& multi_spaced);
    void multi_test_ISSH_single(const std::vector<SpacedQmer>& multi_spaced);
    void multi_test_speedup_multi_previous(const std::vector<SpacedQmer>& multi_spaced);
    void multi_test_ISSH_multi_v1(const std::vector<SpacedQmer>& multi_spaced);
    void multi_test_ISSH_multi_col(const MultiSeedInfo& infoCol);
    void multi_test_ISSH_multi_col_parallel(const MultiSeedInfo& infoCol);
    void multi_test_ISSH_multi_row(const MultiSeedInfoRow& infoRow);
    
    void log_test_results(const SpacedQmer& spaced); // New method for logging results
};

#endif /* TEST_TEST_H_ */
