#ifndef SPACEDQMER_H
#define SPACEDQMER_H

#include <string>
#include <vector>
#include <set>
#include <map>
#include <coin/CbcModel.hpp>
#include <coin/OsiClpSolverInterface.hpp>

class SpacedQmer {
public:
    SpacedQmer(const std::vector<std::string>& seedPatterns);
    std::vector<std::string> generateSpacedQmers(const std::string& sequence) const;

private:
    std::vector<std::string> seedPatterns;
    std::vector<std::vector<int>> positions;

    void parseSeedPatterns();
    std::vector<int> parseSinglePattern(const std::string& pattern) const;
    std::set<int> linearProgrammingSetCovering(const std::string& sequence) const;
};

#endif // SPACEDQMER_H
