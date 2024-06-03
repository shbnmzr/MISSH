#include "SpacedQmer.h"
#include <algorithm>

SpacedQmer::SpacedQmer(const std::vector<std::string>& seedPatterns) : seedPatterns(seedPatterns) {
    parseSeedPatterns();
}

void SpacedQmer::parseSeedPatterns() {
    for (const auto& pattern : seedPatterns) {
        positions.push_back(parseSinglePattern(pattern));
    }
}

std::vector<int> SpacedQmer::parseSinglePattern(const std::string& pattern) const {
    std::vector<int> pos;
    for (int i = 0; i < pattern.size(); ++i) {
        if (pattern[i] == '1') {
            pos.push_back(i);
        }
    }
    return pos;
}

std::vector<std::string> SpacedQmer::generateSpacedQmers(const std::string& sequence) const {
    std::vector<std::string> spacedQmers;
    std::set<int> selectedSeeds = linearProgrammingSetCovering(sequence);
    for (int i : selectedSeeds) {
        for (const auto& pos : positions) {
            std::string qmer;
            for (int p : pos) {
                if (i + p < sequence.size()) {
                    qmer += sequence[i + p];
                }
            }
            spacedQmers.push_back(qmer);
        }
    }

    return spacedQmers;
}
