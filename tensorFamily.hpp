#ifndef tensorFamilyHead
#define tensorFamilyHead

#include <vector>
#include <algorithm>
#include "healpix_cxx/healpix_map.h"

struct tensorFamily {    
    const uint rankA{0}, rankB{0};
    const uint curvIndex{0};
    pointing r{1.5701963268,0};
    virtual ~tensorFamily() = default;
    virtual double accessElement(const std::vector<uint_fast8_t>& indices) const = 0;
    tensorFamily(uint rankA, uint rankB, uint curvIndex=0) : rankA(rankA), rankB(rankB), curvIndex(curvIndex) {}
    tensorFamily(uint rankA, uint rankB, uint curvIndex, const pointing& r) : rankA(rankA), rankB(rankB), curvIndex(curvIndex), r(r) {}
};

#endif
