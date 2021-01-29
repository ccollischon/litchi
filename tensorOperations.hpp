#ifndef tensorOperations
#define tensorOperations

#include "tensorFamily.hpp"
#include <iostream>
#include <type_traits>
#include <vector>
#include <stdexcept>

// Sum expression template, operator+ of minkTensorIntegrand should return this
template <typename ltype, typename rtype>
struct minkTensorSum : tensorFamily
{
    const rtype rhs;
    const ltype lhs;
    
    minkTensorSum(const ltype& left,const rtype& right)  : tensorFamily(right.rankA, right.rankB, right.curvIndex), 
        rhs(right), lhs(left)
    {
        if( (right.rankA != left.rankA) || (right.rankB != left.rankB) || (right.curvIndex != left.curvIndex))
        {
            std::cerr << "Error: trying to add tensors of different type. Left has (rankA, rankB, curvIndex) =  (" << left.rankA << "," << left.rankB << "," << left.curvIndex << "), Right has ("  << right.rankA << "," << right.rankB << "," << right.curvIndex << "), this makes no sense" << std::endl;
            throw std::invalid_argument( "minkTensorSum: Different parameters in addition" );
        }
    }
    
    double accessElement(std::vector<uint_fast8_t> indices) const
    {
        return (rhs.accessElement(indices)) + (lhs.accessElement(indices));
    }
};

//like minkTensorSum, but with * instead of + and scalar needs to come after tensor
template <typename tensortype, typename scalar>
struct minkTensorTimes : tensorFamily
{
    const tensortype mytensor;
    const scalar myscalar;
    
    template <typename tens, typename scal>
    minkTensorTimes(const tens& te, const scal& sc)  : tensorFamily(te),
        mytensor(te),
        myscalar(sc)
    {
    }
    
    double accessElement(std::vector<uint_fast8_t> indices) const
    {
        return mytensor.accessElement(indices) * myscalar;
    }
};

//Sum only for tensorFamily
template<typename left, typename right , typename std::enable_if_t<std::is_base_of<tensorFamily,left>::value && std::is_base_of<tensorFamily,left>::value>* = nullptr>
minkTensorSum<left,right> operator + (const left& lhs, const right& rhs)
{
    minkTensorSum<left,right> returnval (lhs, rhs);
    return returnval;
}

//Tensor times scalar, once with scalar on right and once with scalar on left
template<typename left, typename right , typename std::enable_if_t<std::is_arithmetic<right>::value && std::is_base_of<tensorFamily,left>::value>* = nullptr>
minkTensorTimes<left,right> operator* (const left& lhs, const right& rhs)
{
    return minkTensorTimes<left, right> (lhs, rhs);
}

template<typename left, typename right , typename std::enable_if_t<std::is_arithmetic<left>::value && std::is_base_of<tensorFamily,right>::value>* = nullptr>
minkTensorTimes<right,left> operator* (const left& lhs, const right& rhs)
{
    return minkTensorTimes<right, left> (rhs, lhs);
}




#endif
