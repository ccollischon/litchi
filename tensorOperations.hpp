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
    
    minkTensorSum(const ltype& left,const rtype& right)  : tensorFamily(right.rankA, right.rankB, right.curvIndex, right.r), 
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

/***** Functions that work on tensors ****/

double trace(tensorFamily& input) //sum of eigenvalues
{
    double sinT = sin(input.r.theta);
    if(input.rankA+input.rankB == 0) return input.accessElement({});
    
    std::vector<uint_fast8_t> indices(input.rankA+input.rankB,0);
    double summand = input.accessElement(indices); //zeroes
    
    for(uint i=0; i<indices.size(); i++) { indices.at(i) = 1; }
    summand += input.accessElement(indices)*pow(sinT, input.rankA+input.rankB); //ones with metric contribution
    return summand;
}

double eigenValueQuotient(tensorFamily& input) //TODO check
{
    uint ranksum = input.rankA+input.rankB;
    if (ranksum == 1)
    {
        std::cerr << "Error: Eigenvalue quotient not well-defined for rank 1! Trying to calculate rankA rankB = " << input.rankA <<" "<< input.rankB << std::endl;
        throw std::invalid_argument("eigenValueQuotient not defined for rank 1");
    } else if (ranksum == 0)
    {
        return input.accessElement({});
    } else if(ranksum == 2)
    { //eigenvalues of matrix (a b, c d)
        double sin2T = sin(input.r.theta)*sin(input.r.theta); //pull down one index = sin^2 (theta) factor wherever left index = 1 (arbitrary choice)
        double dplusa = input.accessElement({1,1})*sin2T+input.accessElement({0,0});
        double adminusbc = input.accessElement({0,0})*input.accessElement({1,1})*sin2T - pow(input.accessElement({0,1}),2)*sin2T; //tensors are symmetric here, one of them needs factor
        double twolambda1 = dplusa + sqrt(dplusa*dplusa - 4*adminusbc);
        double twolambda2 = dplusa - sqrt(dplusa*dplusa - 4*adminusbc);
        if(twolambda1>twolambda2) return twolambda1/twolambda2;
        else return twolambda2/twolambda1;
        
    } else
    {
        std::cerr << "Error: Eigenvalue quotient not implemented for rank higher than 2! Trying to calculate rankA rankB = " << input.rankA <<" "<< input.rankB << std::endl;
        throw std::invalid_argument("eigenValueQuotient not implemented for higher ranks");
    }
}





#endif
