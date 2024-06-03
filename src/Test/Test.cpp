/*
 * Test.cpp
 *
 */

#include "Test.h"

Test::Test() {
    // TODO Auto-generated constructor stub
}

Test::~Test() {
    // TODO Auto-generated destructor stub
}

bool Test::load_sequences(const FileParameter& file) {
    this->to_hash.clear();
    this->to_hash.shrink_to_fit();
    cout << "\nLoad sequences... " << flush;
    bool ret = true;
    FilesScan scans;
    if (scans.init(file.getInputFiles().first.getPath(), file.getInputFiles().second.getPath())) {
        auto get_sequence = [](FileScan& filescan, vector<string>& to_hash) {
            size_t index_start = to_hash.size();
            to_hash.resize(to_hash.size() + filescan.getSequenceNumber());
            bool ret = true;
            for (size_t i = 0; i < filescan.getSequenceNumber() && ret; ++i) {
                Sequence seq;
                if (filescan.getSequenceWithIndex(i, seq))
                    to_hash[index_start + i] = seq.getSequence();
                else
                    ret = false;
            }
            return ret;
        };
        ret &= get_sequence(scans.first, this->to_hash);
        ret &= get_sequence(scans.second, this->to_hash);
    } else {
        ret = false;
    }

    // remove not valid
    auto remove_not_valid = [](string& seq) {
        for (size_t i = 0; i < seq.size(); ++i)
            if (CharToInt(seq[i]) == 4)
                return true;
        return false;
    };
    auto it = remove_if(this->to_hash.begin(), this->to_hash.end(), remove_not_valid);
    this->to_hash.resize(distance(this->to_hash.begin(), it));

    // count size mean read
    this->read_size_mean = 0;
    for (size_t i = 0; i < this->to_hash.size(); ++i)
        this->read_size_mean = this->read_size_mean + ((int64_t)this->to_hash[i].size() - (int64_t)this->read_size_mean) / (int64_t)(i + 1);

    auto min_max = minmax_element(this->to_hash.begin(), this->to_hash.end(), [](string& a, string& b) { return a.size() < b.size(); });
    this->read_size_min = min_max.first->size();
    this->read_size_max = min_max.second->size();

    cout << "Complete" << endl << flush;
    return ret;
}

bool Test::single_run(const SpacedQmer& spaced, bool test_equals) {
    this->single_equals = true;
    if (test_equals)
        this->single_test_equals(spaced);
    if (this->single_equals)
        this->single_test_hashes(spaced);
    return this->single_equals;
}

void Test::single_save(const FileParameter& file, string dir_output) {
    createDirAndSubDir(dir_output);

    auto naive = duration_cast<milliseconds>(this->times_naive.second - this->times_naive.first);
    auto speedup_previous = duration_cast<milliseconds>(this->times_speedup_previous.second - this->times_speedup_previous.first);
    auto ISSH = duration_cast<milliseconds>(this->times_ISSH.second - this->times_ISSH.first);
    auto lp_based = duration_cast<milliseconds>(this->times_linear_programming_set_covering.second - this->times_linear_programming_set_covering.first);

    string output_file = dir_output + "times.csv";
    bool exist = file_exist(output_file);
    ofstream output;
    if (exist)
        output.open(output_file, ofstream::out | ofstream::app);
    else {
        output.open(output_file, ofstream::out);
        output << "Files,n_seq,m_size,t_naive,t_previous,t_ISSH,t_lp_based," << endl << flush;
    }
    string str_equals = this->single_equals ? "yes" : "no";
    output << file.getInputFiles().getIdentify() << "," << this->to_hash.size() << "," << this->read_size_mean << "," << naive.count() << "," << speedup_previous.count() << "," << ISSH.count() << "," << lp_based.count() << "," << str_equals << "," << endl << flush;
    output.close();
}

void Test::single_test_hashes(const SpacedQmer& spaced) {
    Chrono chrono;

    cout << "Test naive... " << flush;
    chrono.exe(this, &Test::single_test_naive, spaced, this->times_naive);
    cout << "Complete" << endl << flush;
    cout << "Test speedup previous... " << flush;
    chrono.exe(this, &Test::single_test_speedup_previous, spaced, this->times_speedup_previous);
    cout << "Complete" << endl << flush;
    cout << "Test ISSH... " << flush;
    chrono.exe(this, &Test::single_test_ISSH, spaced, this->times_ISSH);
    cout << "Complete" << endl << flush;

    /**
     * Following line will be here for saftey and following potential use
     * please keep just in case.
    */
    // const SpacedQmerLP* spacedLP = dynamic_cast<const SpacedQmerLP*>(&spaced);

    cout << "Test linear programming set covering... " << flush;
    chrono.exe(this, &Test::single_test_linear_programming_set_covering, spaced, this->times_linear_programming_set_covering);
    cout << "Complete" << endl << flush;
}

void Test::single_test_equals(const SpacedQmer& spaced) {
    cout << "Checking if the hashes calculated with the different algorithms are equal..." << flush;
    Hash_Err_V vHash_naive;
    Hash_Err_V vHash_speedup_previous;
    Hash_Err_V vHash_ISSH;
    Hash_Err_V vHash_LP;
    // per ciascuna delle sequenza intervallate dai delimitatori calcola gli hash
    for (size_t i = 0; i < this->to_hash.size() && this->single_equals; ++i) {
        GetHashes_naive(this->to_hash[i], spaced, vHash_naive, CharToInt);
        GetHashes_speedup_previous(this->to_hash[i], spaced, vHash_speedup_previous, CharToInt);
        GetHashes_with_ISSH(this->to_hash[i], spaced, vHash_ISSH, CharToInt);

        const SpacedQmerLP* spacedLP = dynamic_cast<const SpacedQmerLP*>(&spaced);
        if (spacedLP) {
            GetHashes_with_LP(this->to_hash[i], *spacedLP, vHash_LP, CharToInt);
        }

        this->single_equals &= (vHash_naive.size() == vHash_speedup_previous.size()
            && vHash_naive.size() == vHash_ISSH.size()
            && (!spacedLP || vHash_naive.size() == vHash_LP.size()));
        if (this->single_equals) {
            for (size_t j = 0; j < vHash_naive.size(); ++j) {
                this->single_equals = (vHash_naive[j].isCorrect() == vHash_speedup_previous[j].isCorrect()
                    && vHash_naive[j].isCorrect() == vHash_ISSH[j].isCorrect()
                    && (!spacedLP || vHash_naive[j].isCorrect() == vHash_LP[j].isCorrect()));
                if (this->single_equals) {
                    this->single_equals &= (vHash_naive[j].hash == vHash_speedup_previous[j].hash
                        && vHash_naive[j].hash == vHash_ISSH[j].hash
                        && (!spacedLP || vHash_naive[j].hash == vHash_LP[j].hash));
                    if (!this->single_equals)
                        cout << "";
                } else
                    cout << "";
            }
        } else
            cout << "";
    }
    if (this->single_equals)
        cout << "Complete" << endl << flush;
    else
        cout << "Error" << endl << flush;
}

void Test::single_test_naive(const SpacedQmer& spaced) {
    for (size_t i = 0; i < this->to_hash.size(); ++i) {
        Hash_Err_V vHash;
        GetHashes_naive(this->to_hash[i], spaced, vHash, CharToInt);
    }
}

void Test::single_test_speedup_previous(const SpacedQmer& spaced) {
    for (size_t i = 0; i < this->to_hash.size(); ++i) {
        Hash_Err_V vHash;
        GetHashes_speedup_previous(this->to_hash[i], spaced, vHash, CharToInt);
    }
}

void Test::single_test_ISSH(const SpacedQmer& spaced) {
    for (size_t i = 0; i < this->to_hash.size(); ++i) {
        Hash_Err_V vHash;
        GetHashes_with_ISSH(this->to_hash[i], spaced, vHash, CharToInt);
    }
}

void Test::single_test_linear_programming_set_covering(const SpacedQmer& spaced) {
    const SpacedQmerLP* spacedLP = dynamic_cast<const SpacedQmerLP*>(&spaced);
    if (spacedLP) {
        for (size_t i = 0; i < this->to_hash.size(); ++i) {
            Hash_Err_V vHash;
            GetHashes_with_LP(this->to_hash[i], *spacedLP, vHash, CharToInt);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/* Da qui iniziano le funzioni che riguardano l'hashing di più spaced seed
contemporaneamente. Queste funioni non sono state implementate per l'algoritmo
ISSH*/

bool Test::multi_run(const vector<SpacedQmer>& multi_spaced, bool test_equals) {
    // initialize the complete structure for recovering hashes using all the seeds together
    this->spaced_qmers.init(multi_spaced);

    // if asked it performs the check between the results found
    this->multi_equals = true;
    if (test_equals)
        this->multi_test_equals(multi_spaced);

    if (this->multi_equals)
        this->multi_test_hashes(multi_spaced);
    return this->multi_equals;
}

void Test::multi_save(const FileParameter& file, const vector<SpacedQmer>& multi_spaced, string dir_output) {
    createDirAndSubDir(dir_output);

    auto naive = duration_cast<milliseconds>(this->times_multi_naive.second - this->times_multi_naive.first);
    auto speedup_previous = duration_cast<milliseconds>(this->times_multi_speedup_previous.second - this->times_multi_speedup_previous.first);
    auto ISSH_single = duration_cast<milliseconds>(this->times_multi_ISSH_single.second - this->times_multi_ISSH_single.first);
    auto speedup_multi_previous = duration_cast<milliseconds>(this->times_multi_speedup_multi_previous.second - this->times_multi_speedup_multi_previous.first);
    auto ISSH_multi_v1 = duration_cast<milliseconds>(this->times_multi_ISSH_multi_v1.second - this->times_multi_ISSH_multi_v1.first);
    auto ISSH_multi_col = duration_cast<milliseconds>(this->times_multi_ISSH_multi_col.second - this->times_multi_ISSH_multi_col.first);
    auto ISSH_multi_col_parallel = duration_cast<milliseconds>(this->times_multi_ISSH_multi_col_parallel.second - this->times_multi_ISSH_multi_col_parallel.first);
    auto ISSH_multi_row = duration_cast<milliseconds>(this->times_multi_ISSH_multi_row.second - this->times_multi_ISSH_multi_row.first);

    string output_file = dir_output + "multi_times.csv";
    bool exist = file_exist(output_file);
    ofstream output;
    if (exist)
        output.open(output_file, ofstream::out | ofstream::app);
    else {
        output.open(output_file, ofstream::out);
        output << "Spaced" << flush;
        for (size_t i = 0; i < multi_spaced.size(); ++i)
            output << "," << multi_spaced[i].toString() << "," << endl << flush;

        output << "Files,n_seq,m_size,t_naive,t_previous,ISSH_single,t_multi_previous,ISSH_multi_v1,ISSH_multi_col,ISSH_multi_col_parallel, ISSH_multi_row,equals" << endl << flush;
    }
    string str_equals = this->multi_equals ? "yes" : "no";
    output << file.getInputFiles().getIdentify() << "," << this->to_hash.size() << "," << this->read_size_mean << "," << naive.count() << "," << speedup_previous.count() << "," << ISSH_single.count() << "," << speedup_multi_previous.count() << "," << ISSH_multi_v1.count() << "," << ISSH_multi_col.count() << "," << ISSH_multi_col_parallel.count() << "," << ISSH_multi_row.count() << "," << str_equals << "," << endl << flush;
    output.close();

    output_file = dir_output + "changes.csv";
    output.open(output_file, ofstream::out);
    output << "changes_previous," << endl << flush;
    for (size_t s = 0; s < multi_spaced.size(); ++s) {
        const V_PreviusShift& curr_shift = multi_spaced[s].GetShiftMinChange();
        for (size_t i = 0; i < curr_shift.size(); ++i)
            output << curr_shift[i].GetSize() << "," << flush;
        output << endl << flush;
    }
    output << endl << flush;
    output << "changes_multi_previous," << endl << flush;
    const vector<PreviusShift_Ext_V>& v_shift_min = this->spaced_qmers.getShiftMin();
    for (size_t s = 0; s < v_shift_min.size(); ++s) {
        for (size_t i = 0; i < v_shift_min[s].size(); ++i)
            output << v_shift_min[s][i].GetSize() << "," << flush;
        output << endl << flush;
    }
    output << endl << flush;
    output << "changes_multi_previous_rotated," << endl << flush;
    const vector<PreviusShift_Ext_V>& v_shift_min_rotated = this->spaced_qmers.getShiftMinRotated();
    for (size_t i = 0; i < v_shift_min_rotated.size(); ++i) {
        for (size_t s = 0; s < v_shift_min_rotated[i].size(); ++s)
            output << v_shift_min_rotated[i][s].GetSize() << "," << flush;
        output << endl << flush;
    }
    output.close();
}

void Test::multi_test_hashes(const vector<SpacedQmer>& multi_spaced) {
    Chrono chrono;

    cout << "Test naive... " << flush;
    chrono.exe(this, &Test::multi_test_naive, multi_spaced, this->times_multi_naive);
    cout << "Complete" << endl << flush;

    cout << "Test previous... " << flush;
    chrono.exe(this, &Test::multi_test_speedup_previous, multi_spaced, this->times_multi_speedup_previous);
    cout << "Complete" << endl << flush;

    cout << "Test ISSH single... " << flush;
    chrono.exe(this, &Test::multi_test_ISSH_single, multi_spaced, this->times_multi_ISSH_single);
    cout << "Complete" << endl << flush;

    cout << "Test multi previous... " << flush;
    chrono.exe(this, &Test::multi_test_speedup_multi_previous, multi_spaced, this->times_multi_speedup_multi_previous);
    cout << "Complete" << endl << flush;

    cout << "Test ISSH multi v1... " << flush;
    chrono.exe(this, &Test::multi_test_ISSH_multi_v1, multi_spaced, this->times_multi_ISSH_multi_v1);
    cout << "Complete" << endl << flush;

    MultiSpacedQmer MultiSeed(multi_spaced);
    MultiSeedInfo infoCol = MultiSeed.Get_multi_seed_info_col();

    cout << "Test ISSH multi col... " << flush;
    chrono.exe(this, &Test::multi_test_ISSH_multi_col, infoCol, this->times_multi_ISSH_multi_col);
    cout << "Complete" << endl << flush;

    cout << "Test ISSH multi col parallel... " << flush;
    chrono.exe(this, &Test::multi_test_ISSH_multi_col_parallel, infoCol, this->times_multi_ISSH_multi_col_parallel);
    cout << "Complete" << endl << flush;

    MultiSeedInfoRow infoRow = MultiSeed.Get_multi_seed_info_row();

    cout << "Test ISSH multi row... " << flush;
    chrono.exe(this, &Test::multi_test_ISSH_multi_row, infoRow, this->times_multi_ISSH_multi_row);
    cout << "Complete" << endl << flush;

    // Add LP-based hashing test
    bool hasLP = false;
    for (const auto& spaced : multi_spaced) {
        if (dynamic_cast<const SpacedQmerLP*>(&spaced)) {
            hasLP = true;
            break;
        }
    }

    if (hasLP) {
        cout << "Test LP-based hashing... " << flush;
        chrono.exe(this, &Test::multi_test_LP_based, multi_spaced, this->times_multi_LP_based);
        cout << "Complete" << endl << flush;
    }
}

void Test::multi_test_equals(const vector<SpacedQmer>& multi_spaced) {
    cout << "Multi test equals hashes... " << flush;
    vector<vector<Position>> vErr_multi_previous(multi_spaced.size());

    V_V_Hash_Err vHashes_multi_previous;
    Hash_Err_V_V vHashes_single(multi_spaced.size());
    Hash_Err_V_V vHashes_ISSH_single(multi_spaced.size());
    Hash_Err_V_V vHashes_ISSH_multi_v1(multi_spaced.size());
    Hash_Err_V_V vHashes_ISSH_multi_col(multi_spaced.size());
    Hash_Err_V_V vHashes_ISSH_multi_col_parallel(multi_spaced.size());
    Hash_Err_V_V vHashes_ISSH_multi_row(multi_spaced.size());

    vector<V_V_PreviusShift> VV_shifts;
    VV_shifts.resize(multi_spaced.size());
    vector<Position> v_pos_one;
    v_pos_one.resize(multi_spaced.size());
    size_t max_transient_length = 0;
    size_t curr_max;
    for (size_t j = 0; j < multi_spaced.size(); ++j) {
        VV_shifts[j] = multi_spaced[j].GetMultipleShifts();
        v_pos_one[j] = multi_spaced[j].GetPosOne();
        curr_max = VV_shifts[j].size();
        max_transient_length = curr_max > max_transient_length ? curr_max : max_transient_length;
    }

    MultiSpacedQmer MultiSeed(multi_spaced);
    MultiSeedInfo infoCol = MultiSeed.Get_multi_seed_info_col();
    MultiSeedInfoRow infoRow = MultiSeed.Get_multi_seed_info_row();

    for (size_t i = 0; i < this->to_hash.size() && this->multi_equals; ++i) {
        GetHashes_speedup_multi_previous_Rotated(this->to_hash[i], this->spaced_qmers, vHashes_multi_previous, CharToInt);

        for (size_t j = 0; j < multi_spaced.size(); ++j)
            GetHashes_speedup_previous(this->to_hash[i], multi_spaced[j], vHashes_single[j], CharToInt);

        for (size_t j = 0; j < multi_spaced.size(); ++j)
            GetHashes_with_ISSH(this->to_hash[i], multi_spaced[j], vHashes_ISSH_single[j], CharToInt);

        GetHashes_with_ISSH_multi_v1(this->to_hash[i], multi_spaced, VV_shifts, v_pos_one, max_transient_length, vHashes_ISSH_multi_v1, CharToInt);

        GetHashes_with_ISSH_multi_col(this->to_hash[i], infoCol, vHashes_ISSH_multi_col, CharToInt);

        GetHashes_with_ISSH_multi_col_parallel(this->to_hash[i], infoCol, vHashes_ISSH_multi_col_parallel, CharToInt);

        GetHashes_with_ISSH_multi_row(this->to_hash[i], infoRow, vHashes_ISSH_multi_row, CharToInt);

        for (size_t j = 0; j < multi_spaced.size(); ++j) {
            vector<V_V_Hash_Err::reference_value_type> multi_previous;
            vHashes_multi_previous.get(j, multi_previous);

            this->multi_equals &= multi_previous.size() == vHashes_single[j].size();
            this->multi_equals &= vHashes_ISSH_single[j].size() == vHashes_single[j].size();
            this->multi_equals &= vHashes_ISSH_multi_v1[j].size() == vHashes_single[j].size();
            this->multi_equals &= vHashes_ISSH_multi_col[j].size() == vHashes_single[j].size();
            this->multi_equals &= vHashes_ISSH_multi_col_parallel[j].size() == vHashes_single[j].size();
            this->multi_equals &= vHashes_ISSH_multi_row[j].size() == vHashes_single[j].size();

            if (!this->multi_equals) {
                cout << "Error size j=" << j
                    << ", FSH multi =" << multi_previous.size()
                    << ", FSH single=" << vHashes_single[j].size()
                    << ", ISSH single=" << vHashes_ISSH_single[j].size()
                    << ", ISSH multi v1=" << vHashes_ISSH_multi_v1[j].size()
                    << ", ISSH multi col=" << vHashes_ISSH_multi_col[j].size()
                    << ", ISSH multi col parallel=" << vHashes_ISSH_multi_col_parallel[j].size()
                    << ", ISSH multi row=" << vHashes_ISSH_multi_row[j].size()
                    << endl
                    << flush;
            }
        }

        if (this->multi_equals) {
            for (size_t j = 0; j < multi_spaced.size(); ++j) {
                vector<V_V_Hash_Err::reference_value_type> multi_previous;
                vHashes_multi_previous.get(j, multi_previous);

                for (size_t h = 0; h < vHashes_single[j].size() && this->multi_equals; ++h) {
                    this->multi_equals &= multi_previous[h].get().isCorrect() == vHashes_single[j][h].isCorrect();
                    this->multi_equals &= vHashes_ISSH_single[j][h].isCorrect() == vHashes_single[j][h].isCorrect();
                    this->multi_equals &= vHashes_ISSH_multi_v1[j][h].isCorrect() == vHashes_single[j][h].isCorrect();
                    this->multi_equals &= vHashes_ISSH_multi_col[j][h].isCorrect() == vHashes_single[j][h].isCorrect();
                    this->multi_equals &= vHashes_ISSH_multi_col_parallel[j][h].isCorrect() == vHashes_single[j][h].isCorrect();
                    this->multi_equals &= vHashes_ISSH_multi_row[j][h].isCorrect() == vHashes_single[j][h].isCorrect();

                    if (this->multi_equals) {
                        this->multi_equals &= multi_previous[h].get().hash == vHashes_single[j][h].hash;
                        this->multi_equals &= vHashes_ISSH_single[j][h].hash == vHashes_single[j][h].hash;
                        this->multi_equals &= vHashes_ISSH_multi_v1[j][h].hash == vHashes_single[j][h].hash;
                        this->multi_equals &= vHashes_ISSH_multi_col[j][h].hash == vHashes_single[j][h].hash;
                        this->multi_equals &= vHashes_ISSH_multi_col_parallel[j][h].hash == vHashes_single[j][h].hash;
                        this->multi_equals &= vHashes_ISSH_multi_row[j][h].hash == vHashes_single[j][h].hash;

                        if (!this->multi_equals) {
                            cout << "Error hash j=" << j << ", h=" << h
                                << ", FSH multi=" << multi_previous[h].get().hash
                                << ", FSH single=" << vHashes_single[j][h].hash
                                << ", ISSH single=" << vHashes_ISSH_single[j][h].hash
                                << ", ISSH multi v1=" << vHashes_ISSH_multi_v1[j][h].hash
                                << ", ISSH multi col=" << vHashes_ISSH_multi_col[j][h].hash
                                << ", ISSH multi col parallel=" << vHashes_ISSH_multi_col_parallel[j][h].hash
                                << ", ISSH multi row=" << vHashes_ISSH_multi_row[j][h].hash
                                << endl
                                << flush;

                            cout << "i=" << i << "," << this->to_hash[i] << endl << flush;
                        }
                    } else {
                        cout << "Error charInt j=" << j << ", h=" << h
                            << ", FSH multi=" << multi_previous[h].get().hash
                            << ", FSH single=" << vHashes_single[j][h].hash
                            << ", ISSH single=" << vHashes_ISSH_single[j][h].hash
                            << ", ISSH multi v1=" << vHashes_ISSH_multi_v1[j][h].hash
                            << ", ISSH multi col=" << vHashes_ISSH_multi_col[j][h].hash
                            << ", ISSH multi row=" << vHashes_ISSH_multi_row[j][h].hash
                            << endl
                            << flush;
                        cout << "i=" << i << "," << this->to_hash[i] << endl << flush;
                    }
                }
            }
        } else {
            cout << "";
        }
    }
    if (this->multi_equals)
        cout << "Complete" << endl << flush;
    else
        cout << "Error" << endl << flush;
}

void Test::multi_test_LP_based(const vector<SpacedQmer>& multi_spaced) {
    for (size_t i = 0; i < this->to_hash.size(); ++i) {
        for (size_t j = 0; j < multi_spaced.size(); ++j) {
            const SpacedQmerLP* spacedLP = dynamic_cast<const SpacedQmerLP*>(&multi_spaced[j]);
            if (spacedLP) {
                Hash_Err_V vHash;
                GetHashes_with_LP(this->to_hash[i], *spacedLP, vHash, CharToInt);
            }
        }
    }
}

void Test::multi_test_naive(const vector<SpacedQmer>& multi_spaced) {
    Hash_Err_V_V hash_v(multi_spaced.size());
    for (size_t i = 0; i < this->to_hash.size(); ++i) {
        for (size_t j = 0; j < multi_spaced.size(); ++j) {
            GetHashes_naive(this->to_hash[i], multi_spaced[j], hash_v[j], CharToInt);
        }
    }
}

void Test::multi_test_speedup_previous(const vector<SpacedQmer>& multi_spaced) {
    Hash_Err_V_V hash_v(multi_spaced.size());
    for (size_t i = 0; i < this->to_hash.size(); ++i)
        for (size_t j = 0; j < multi_spaced.size(); ++j)
            GetHashes_speedup_previous(this->to_hash[i], multi_spaced[j], hash_v[j], CharToInt);
}

void Test::multi_test_ISSH_single(const vector<SpacedQmer>& multi_spaced) {
    Hash_Err_V_V hash_v(multi_spaced.size());
    for (size_t i = 0; i < this->to_hash.size(); ++i) {
        for (size_t j = 0; j < multi_spaced.size(); ++j) {
            GetHashes_with_ISSH(this->to_hash[i], multi_spaced[j], hash_v[j], CharToInt);
        }
    }
}

void Test::multi_test_speedup_multi_previous(const vector<SpacedQmer>& multi_spaced) {
    V_V_Hash_Err hashes_errs;
    hashes_errs.reserve(multi_spaced.size(), this->read_size_max);
    for (size_t i = 0; i < this->to_hash.size(); ++i)
        GetHashes_speedup_multi_previous_Rotated(this->to_hash[i], this->spaced_qmers, hashes_errs, CharToInt);
}

void Test::multi_test_ISSH_multi_v1(const vector<SpacedQmer>& multi_spaced) {
    Hash_Err_V_V hash_v(multi_spaced.size());
    vector<V_V_PreviusShift> VV_shifts;
    VV_shifts.resize(multi_spaced.size());
    vector<Position> v_pos_one;
    v_pos_one.resize(multi_spaced.size());
    size_t max_transient_length = 0;
    size_t curr_max;
    for (size_t j = 0; j < multi_spaced.size(); ++j) {
        VV_shifts[j] = multi_spaced[j].GetMultipleShifts();
        v_pos_one[j] = multi_spaced[j].GetPosOne();
        curr_max = VV_shifts[j].size();
        max_transient_length = curr_max > max_transient_length ? curr_max : max_transient_length;
    }
    for (size_t i = 0; i < this->to_hash.size(); ++i)
        GetHashes_with_ISSH_multi_v1(this->to_hash[i], multi_spaced, VV_shifts, v_pos_one, max_transient_length, hash_v, CharToInt);
}

void Test::multi_test_ISSH_multi_col(const MultiSeedInfo& infoCol) {
    Hash_Err_V_V hash_v(infoCol.size());
    for (size_t i = 0; i < this->to_hash.size(); ++i)
        GetHashes_with_ISSH_multi_col(this->to_hash[i], infoCol, hash_v, CharToInt);
}

void Test::multi_test_ISSH_multi_col_parallel(const MultiSeedInfo& infoCol) {
    Hash_Err_V_V hash_v(infoCol.size());
    for (size_t i = 0; i < this->to_hash.size(); ++i) {
        GetHashes_with_ISSH_multi_col_parallel(this->to_hash[i], infoCol, hash_v, CharToInt);
    }
}

void Test::multi_test_ISSH_multi_row(const MultiSeedInfoRow& infoRow) {
    Hash_Err_V_V hash_v(infoRow.info.size());
    for (size_t i = 0; i < this->to_hash.size(); ++i)
        GetHashes_with_ISSH_multi_row(this->to_hash[i], infoRow, hash_v, CharToInt);
}
