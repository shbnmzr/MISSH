//============================================================================
// Name        : MISSH.cpp
//============================================================================

#include "Test/Test.h"
#include "Spaced/SetCoveringSolver.h" // Include the set covering solver
#include "Spaced/SpacedQmer.h"
#include <set>
#include <vector>
#include <iostream>

int main(int argc, char* argv[]) {
    string dir_output = "../output/";
    bool sequence = false;
    FileParameter param;
    omp_set_num_threads(4);
    // default setting: performs the test multi
    int test_kind = 2;

    for(int i = 1; i < argc; i++) {
        if(strcmp(argv[i], "-si") == 0) {
            i++;
            if(!param.init(argv[i], "")) {
                cerr << endl << "Please enter an input filename single-end: -si <AbsPathFile>" << flush;
                return 0;
            }
            sequence = true;
        } else if(strcmp(argv[i], "-pi") == 0) {
            i++;
            if(!param.init(argv[i], argv[i+1])) {
                cerr << endl << "Please enter an input filenames paired-end: -pi <AbsPathFile1> <AbsPathFile2>\n" << flush;
                return 0;
            }
            sequence = true;
        } else if(strcmp(argv[i], "-q") == 0) {
            i++;
            vector<string> lines;
            string pathQmers(argv[i]);
            if(getLines(pathQmers, lines)) {
                vector<bool> correctQmer(lines.size(), false);
                regex rgx("^1(0|1)*1$");
                for(size_t i = 0; i < lines.size(); i++) {
                    correctQmer[i] = regex_match(lines[i], rgx);
                    if(!correctQmer[i]) {
                        cerr << endl << "Error on " << to_string(i+1) << "° spaced seed. Enter q-mer with 1 at begin and end of the string on input files. "
                             << "Ex. 1**1*11*1. 1 is the simbol considered, any others are not valid simbols.\n" << flush;
                        return 0;
                    } else {
                        param.addSpacedQmer(lines[i], lines[i]);
                    }
                }
            } else {
                cerr << endl << "Please enter a spaced seeds path as -q <AbsPathFile>. Every file's line must contain a spaced seeds.\n" << flush;
                return 0;
            }
        } else if(strcmp(argv[i], "-dirO") == 0) {
            i++;
            if(string(argv[i]) == "") {
                cerr << endl << "Enter valid path for output directory.\n" << flush;
                return 0;
            } else
                dir_output = argv[i];
        } else if(strcmp(argv[i], "-num") == 0) {
            i++;
            const char* s_num_prev(argv[i]);
            if(!isdigit(s_num_prev[0])) {
                cerr << endl << "Please enter valid number of previous hashes from which retrive positions.\n" << flush;
                return 0;
            } else
                param.setNumPrev((size_t)atoi(s_num_prev));
        } else if(strcmp(argv[i], "-test") == 0) {
            i++;
            if(strcmp(argv[i], "single") == 0) {
                test_kind = 0;
            } else {
                if(strcmp(argv[i], "multi") == 0) {
                    test_kind = 1;
                } else {
                    cerr << endl << "Please enter valid value (single or multi), both tests will be performed.\n" << flush;
                }
            }
        } else if(strcmp(argv[i], "-threads") == 0) {
            i++;
            const char* threads_num = argv[i];
            if(!isdigit(threads_num[0])) {
                cerr << endl << "Please enter valid number of previous hashes from which retrive positions.\n" << flush;
                return 0;
            } else
                omp_set_num_threads(atoi(threads_num));
        }
    }

    if(!sequence) {
        cerr << endl << "Please enter a DNA sequence with -si <AbsPathFile> or -pi <AbsPathFile1> <AbsPathFile2>\n" << flush;
        return 0;
    }

    if(param.getVSpaced().empty()) {
        // If no spaced seeds were provided, use the default ones
        param.addSpacedQmer("CLARK-S", "1111011101110010111001011011111");
        param.addSpacedQmer("CLARK-S", "1111101011100101101110011011111");
        param.addSpacedQmer("CLARK-S", "1111101001110101101100111011111");
        param.addSpacedQmer("rasbhari_minimizing_overlap_complexity", "1111010111010011001110111110111");
        param.addSpacedQmer("rasbhari_minimizing_overlap_complexity", "1110111011101111010010110011111");
        param.addSpacedQmer("rasbhari_minimizing_overlap_complexity", "1111101001011100111110101101111");
        param.addSpacedQmer("rasbhari_maximizing_sensitivity", "1111011110011010111110101011011");
        param.addSpacedQmer("rasbhari_maximizing_sensitivity", "1110101011101100110100111111111");
        param.addSpacedQmer("rasbhari_maximizing_sensitivity", "1111110101101011100111011001111");
        cout << endl << "Applied default spaced seed" << flush;
    }

    // Define the universal set and subsets
    std::set<int> universal_set = {0, 1, 2, 3, 4, 5}; // Example universal set
    std::vector<std::set<int>> subsets = {
        {0, 1, 2},
        {2, 3},
        {3, 4, 5},
        {0, 4, 5}
    };

    // Initialize SpacedQmer with set covering
    SpacedQmer spacedQmer;
    spacedQmer.initializeWithSetCovering(universal_set, subsets);

    // Create output directory if not present
    createDirAndSubDir(dir_output);

    cout << endl << "Applied the following spaced seed..." << flush;
    for(size_t i = 0; i < param.getVSpaced().size(); i++)
        cout << endl << "Type:" << param.getVSpaced()[i].first << ", Spaced seed: " << param.getVSpaced()[i].second.toString() << flush;

    // Test for a single spaced seed at a time
    if (test_kind == 0 || test_kind == 2) {
        cout << endl << "Performing test for a single spaced seed at a time" << flush;

        bool single_test_equals = false;
        string dir_output_1 = dir_output + "single/";
        Test test_single;
        if (test_single.load_sequences(param)) {
            for (size_t j = 0; j < param.getVSpaced().size(); ++j) {
                string dir_output_tmp = dir_output_1 + param.getVSpaced()[j].first + param.getVSpaced()[j].second.toString() + "/";

                test_single.single_run(param.getVSpaced()[j].second, single_test_equals);
                test_single.single_save(param, dir_output_tmp);
            }
        }
    }

    // Test for a group of spaced seeds
    if (test_kind == 1 || test_kind == 2) {
        cout << endl << "Performing test for multiple spaced seeds at a time" << flush;

        bool multi_test_equals = false;
        string dir_output_2 = dir_output + "multi/";
        vector<SpacedQmer> multi_spaced;
        for (size_t i = 0; i < param.getVSpaced().size(); ++i)
            multi_spaced.push_back(param.getVSpaced()[i].second);

        Test test_multi;
        if (test_multi.load_sequences(param)) {
            test_multi.multi_run(multi_spaced, multi_test_equals);
            test_multi.multi_save(param, multi_spaced, dir_output_2);
        }
    }
    cout << "End\n";
}
