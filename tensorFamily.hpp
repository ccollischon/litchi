#ifndef tensorFamilyHead
#define tensorFamilyHead

#include <vector>
#include <algorithm>
#include "healpix_cxx/healpix_map.h"

struct tensorFamily {    
    const uint rankA, rankB;
    const uint curvIndex;
    pointing r;
    virtual double accessElement(std::vector<uint_fast8_t> indices) const = 0;
    tensorFamily(uint rankA, uint rankB, uint curvIndex=0) : rankA(rankA), rankB(rankB), curvIndex(curvIndex), r(0,0) {}
    tensorFamily(uint rankA, uint rankB, uint curvIndex, pointing r) : rankA(rankA), rankB(rankB), curvIndex(curvIndex), r(r) {}
};

#endif
