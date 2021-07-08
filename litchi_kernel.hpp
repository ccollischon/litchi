#ifndef litchi_kernel
#define litchi_kernel

#include <vector>
//#include <memory>
#include "healpix_cxx/healpix_map.h"

//#include "tensor2D.hpp"
#include "tensorFamily.hpp"
#include "tensorOperations.hpp"


struct minkmapFamily {
    Healpix_Map<double>& originalMap;
    const uint rankA{}, rankB{};
    const uint curvIndex{};
    virtual ~minkmapFamily() = default;
    virtual minkTensorStack at(int pixnum) const = 0;
    explicit minkmapFamily(Healpix_Map<double>& map) : originalMap(map) {}
    minkmapFamily(Healpix_Map<double>& map, uint rank1, uint rank2, uint curvind) : originalMap(map), rankA(rank1), rankB(rank2), curvIndex(curvind) {}
    minkmapFamily(const minkmapFamily& otherMap) : originalMap(otherMap.originalMap), rankA(otherMap.rankA), rankB(otherMap.rankB), curvIndex(otherMap.curvIndex) {}
};


template <typename ltype, typename rtype> //For sums of minkmaps
struct minkmapSum : minkmapFamily
{
    const rtype rhs;
    const ltype lhs;
    
    minkmapSum(const ltype* left,const rtype* right) :  minkmapFamily(left->originalMap, left->rankA, left->rankB, left->curvIndex), rhs(right), lhs(left)
    {
        if( (right->rankA != left->rankA) || (right->rankB != left->rankB) || (right->curvIndex != left->curvIndex))
        {
            std::cerr << "Error: trying to add maps of tensors of different type. Left has (rankA, rankB, curvIndex) =  (" << left->rankA << "," << left->rankB << "," << left->curvIndex << "), Right has ("  << right->rankA << "," << right->rankB << "," << right->curvIndex << "), this makes no sense" << std::endl;
            throw std::invalid_argument( "minkmapSum: Different parameters in addition" );
        }
    }
    
    minkTensorStack at(int pixnum) const override
    {
        return minkTensorStack(rhs.at(pixnum), lhs.at(pixnum));
    }
    
};

template <typename left, typename right, std::enable_if_t< std::is_base_of<minkmapFamily,left>::value && std::is_base_of<minkmapFamily,right>::value>* = nullptr   > 
minkmapSum<left, right> operator +(const left& a,const right& b)
{
    return minkmapSum<left,right> (&a, &b);
}

template <typename maptype, typename scalar>
struct minkmapTimes : minkmapFamily
{
    const maptype mymap;
    const scalar myscalar;
    
    template <typename map, typename scal>
    minkmapTimes(const map& te, const scal& sc)  : minkmapFamily(te),
        mymap(te),
        myscalar(sc)
    {
    }
    
    minkTensorStack at(int pixel) const override
    {
        return mymap.at(pixel) * myscalar;
    }
};

//Map times scalar, once with scalar on right and once with scalar on left
template<typename left, typename right , typename std::enable_if_t<std::is_arithmetic<right>::value && std::is_base_of<minkmapFamily,left>::value>* = nullptr>
minkmapTimes<left,right> operator* (const left& lhs, const right& rhs)
{
    return minkmapTimes<left, right> (lhs, rhs);
}

template<typename left, typename right , typename std::enable_if_t<std::is_arithmetic<left>::value && std::is_base_of<minkmapFamily,right>::value>* = nullptr>
minkmapTimes<right,left> operator* (const left& lhs, const right& rhs)
{
    return minkmapTimes<right, left> (rhs, lhs);
}

template<typename maptype, typename std::enable_if_t<std::is_base_of<minkmapFamily,maptype>::value>* = nullptr>
struct minkmapStack : minkmapFamily
{
    std::vector<maptype> mapstack;
    
    minkmapStack(std::vector<maptype>& stack) : minkmapFamily(stack.at(0)), mapstack(stack) {}
    
    minkTensorStack at(int pixel) const override
    {
        minkTensorStack thistensor(mapstack.at(0).at(pixel));
        for(uint i=1; i<mapstack.size(); ++i)
        {
            thistensor += mapstack[i].at(pixel);
        }
        return thistensor;
    }
};


#endif
