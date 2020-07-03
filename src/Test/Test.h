/*
 * Test.h
 */

#ifndef Test_Test_H_
#define Test_Test_H_

#include "../Hash/MultiHashFunction.h"
#include "../Input/FileScan.h"
#include "../Chrono/Chrono.h"
#include "../Parameter/Parameter.h"

class Test {
public:
	Test();
	virtual ~Test();
	bool load_sequences(const FileParameter& file);
	bool single_run(const SpacedQmer& spaced, bool test_equals);
	void single_save(const FileParameter& file, string dir_output);

	bool multi_run(const vector<SpacedQmer>& multi_spaced, bool test_equals);
	void multi_save(const FileParameter& file, const vector<SpacedQmer>& multi_spaced, string dir_output);

private:
	typedef pair<vector<SpacedQmer>, vector<PreviusShift_Ext_V>> MultiParameter;

	vector<string> to_hash;
	int64_t read_size_max;
	int64_t read_size_min;
	int64_t read_size_mean;

	Chrono::start_end_time times_naive;
	Chrono::start_end_time times_speedup_previous;
	Chrono::start_end_time times_ISSH;

	Chrono::start_end_time times_multi_naive;
	Chrono::start_end_time times_multi_speedup_previous;
	Chrono::start_end_time times_multi_ISSH_single;	
	Chrono::start_end_time times_multi_speedup_multi_previous;
	Chrono::start_end_time times_multi_ISSH_multi_v1;
	Chrono::start_end_time times_multi_ISSH_multi_col;	
	Chrono::start_end_time times_multi_ISSH_multi_col_parallel;
	Chrono::start_end_time times_multi_ISSH_multi_row;	

	bool single_equals = true;

	//Object defined inside SpacedQmer_Multi.h it contains all the needed information about FSH multiseed
	SpacedQmer_Multi spaced_qmers;
	
	bool multi_equals = true;

	void single_test_hashes(const SpacedQmer& spaced);
	void single_test_equals(const SpacedQmer& spaced);
	void single_test_naive(const SpacedQmer& spaced);
	void single_test_speedup_previous(const SpacedQmer& spaced);//used previus element
	void single_test_ISSH(const SpacedQmer& spaced);



	void multi_test_hashes(const vector<SpacedQmer>& multi_spaced);
	void multi_test_equals(const vector<SpacedQmer>& multi_spaced);
	void multi_test_naive(const vector<SpacedQmer>& multi_spaced);
	void multi_test_speedup_previous(const vector<SpacedQmer>& multi_spaced);//used previus element
	void multi_test_ISSH_single(const vector<SpacedQmer>& multi_spaced);
	void multi_test_speedup_multi_previous(const vector<SpacedQmer>& multi_spaced);//used previus element of different qmer
	void multi_test_ISSH_multi_v1(const vector<SpacedQmer>& multi_spaced);
	void multi_test_ISSH_multi_col(const MultiSeedInfo& infoCol);
	void multi_test_ISSH_multi_col_parallel(const MultiSeedInfo& infoCol);
	void multi_test_ISSH_multi_row(const MultiSeedInfoRow& infoRow);



};

#endif /* Test_Test_H_ */
