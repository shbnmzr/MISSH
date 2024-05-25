//============================================================================
// Name        : ISSH.cpp
//============================================================================

#include "Test/Test.h"

int main(int argc, char* argv[]) {
    string dir_output = "../output/"; // Default output directory
    bool sequence = false; // Flag to check if sequence input is provided
    FileParameter param; // Object to handle file parameters
    omp_set_num_threads(4); // Set the number of threads for parallel processing

    // Default setting: performs the test multi
    int test_kind = 2;

    // Parse command-line arguments
    for(int i=1; i<argc; i++) {
        if(strcmp(argv[i], "-si") == 0) {
            // Single-end input file specified
            i++;
            if(!param.init(argv[i], "")) {
                cerr << endl << "Please enter an input filename single-end: -si <AbsPathFile>" << flush;
                return 0;
            }
            sequence = true;
        }
        else if(strcmp(argv[i], "-pi") == 0) {
            // Paired-end input files specified
            i++;
            if(!param.init(argv[i], argv[i+1])) {
                cerr << endl << "Please enter input filenames paired-end: -pi <AbsPathFile1> <AbsPathFile2>\n" << flush;
                return 0;
            }
            sequence = true;
        }
        else if(strcmp(argv[i], "-q") == 0) {
            // Spaced seeds input file specified
            i++;
            vector<string> lines;
            string pathQmers(argv[i]);
            if(getLines(pathQmers, lines)) {
                vector<bool> correctQmer(lines.size(), false);
                regex rgx("^1(0|1)*1$"); // Regex to check valid q-mers
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
        }
        else if(strcmp(argv[i], "-so") == 0) {
            // Output directory specified
            i++;
            dir_output = argv[i];
        }
        else if(strcmp(argv[i], "-sk") == 0) {
            // Test kind specified
            i++;
            test_kind = atoi(argv[i]);
        }
        else if(strcmp(argv[i], "-h") == 0) {
            // Help flag specified
            cout << endl << "Usage: " << argv[0] << " -si <AbsPathFile> or -pi <AbsPathFile1> <AbsPathFile2>" 
                 << " [-q <AbsPathFile>] [-so <AbsPathDir>] [-sk 0|1|2]" << endl;
            return 0;
        }
    }

    // If no spaced seeds are provided, use default seeds
    if(param.getVSpaced().empty()) {
        param.addSpacedQmer("max_pattern_compression", "110101011101100110100111111111");
        param.addSpacedQmer("rasbhari_maximizing_sensitivity", "1111110101101011100111011001111");
        cout << endl << "Applied default spaced seed" << flush;
    }
    // Create output directory if not present
    createDirAndSubDir(dir_output);

    cout << endl << "Applied the following spaced seed..." << flush;
    for(size_t i = 0; i < param.getVSpaced().size(); i++)
        cout << endl << "Type:" << param.getVSpaced()[i].first << ", Spaced seed: " << param.getVSpaced()[i].second.toString() << flush;

    // Test for a single spaced seed at a time
    if (test_kind == 0 || test_kind == 2) {
        cout << endl << "Performing test for a single spaced seed at a time"  << flush;

        bool single_test_equals = false;
        string dir_output_1 = dir_output + "single/";
        Test test_single;
        if(test_single.load_sequences(param)) {
            for(size_t j = 0; j < param.getVSpaced().size(); ++j) {
                string dir_output_tmp = dir_output_1 + param.getVSpaced()[j].first + param.getVSpaced()[j].second.toString() + "/";
                test_single.single_run(param.getVSpaced()[j].second, single_test_equals);
                test_single.single_save(param, dir_output_tmp);
            }
        }
    }
    // Test for a group of spaced seeds
    if (test_kind == 1 || test_kind == 2) {
        cout << endl << "Performing test for multiple spaced seeds at a time"  << flush;

        bool multi_test_equals = false;
        string dir_output_2 = dir_output + "multi/";
        // multi_spaced vector which contains all the spaced seeds given in input.
        vector<SpacedQmer> multi_spaced;
        for(size_t i = 0; i < param.getVSpaced().size(); ++i)
            multi_spaced.push_back(param.getVSpaced()[i].second);

        Test test_multi;
        if(test_multi.load_sequences(param)) {
            test_multi.multi_run(multi_spaced, multi_test_equals);
            test_multi.multi_save(param, multi_spaced, dir_output_2);
        }
    }
    cout << "End\n";
}
