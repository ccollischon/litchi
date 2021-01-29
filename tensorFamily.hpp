#ifndef tensorFamilyHead
#define tensorFamilyHead

#include <vector>
#include <algorithm>

struct tensorFamily {    
    const uint rankA, rankB;
    const uint curvIndex;
    virtual double accessElement(std::vector<uint_fast8_t> indices) const = 0;
    tensorFamily(uint rankA, uint rankB, uint curvIndex=0) : rankA(rankA), rankB(rankB), curvIndex(curvIndex) {}
};

#endif
