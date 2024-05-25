#ifndef SETCOVERINGSOLVER_H
#define SETCOVERINGSOLVER_H

#include <vector>
#include <set>

class SetCoveringSolver {
public:
    SetCoveringSolver(const std::set<int>& universal_set, const std::vector<std::set<int>>& subsets);
    std::vector<int> solve();

private:
    std::set<int> universal_set_;
    std::vector<std::set<int>> subsets_;
};

#endif // SETCOVERINGSOLVER_H
