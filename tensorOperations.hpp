#ifndef tensorOperations
#define tensorOperations

#include "tensorFamily.hpp"
#include "minkTensorIntegrand.hpp"
#include <iostream>
#include <type_traits>
#include <vector>
#include <stdexcept>

// Sum expression template, operator+ of minkTensorIntegrand should return this
/*
template <typename ltype, typename rtype>
struct minkTensorSum : tensorFamily
{
    rtype rhs;
    ltype lhs;
    
    minkTensorSum(const ltype& left, const rtype& right)  : tensorFamily(right.rankA, right.rankB, right.curvIndex, left.r), 
        rhs(right), lhs(left)
    {
        rhs.moveTo(left.r);
        if( (right.rankA != left.rankA) || (right.rankB != left.rankB) || (right.curvIndex != left.curvIndex))
        {
            std::cerr << "Error: trying to add tensors of different type. Left has (rankA, rankB, curvIndex) =  (" << left.rankA << "," << left.rankB << "," << left.curvIndex << "), Right has ("  << right.rankA << "," << right.rankB << "," << right.curvIndex << "), this makes no sense" << std::endl;
            throw std::invalid_argument( "minkTensorSum: Different parameters in addition" );
        }
        
    }
    
    double accessElement(const std::vector<uint_fast8_t>& indices) const override
    {
        return (rhs.accessElement(indices)) + (lhs.accessElement(indices));
    }
    
    void moveTo(const pointing& newR) override
    {
        rhs.moveTo(newR);
        lhs.moveTo(newR);
    }
    
};

//like minkTensorSum, but with * instead of + and scalar needs to come after tensor
template <typename tensortype, typename scalar>
struct minkTensorTimes : tensorFamily
{
    tensortype mytensor;
    const scalar myscalar;
    
    template <typename tens, typename scal>
    minkTensorTimes(const tens& te, const scal& sc)  : tensorFamily(te),
        mytensor(te),
        myscalar(sc)
    {
    }
    
    double accessElement(const std::vector<uint_fast8_t>& indices) const override
    {
        return mytensor.accessElement(indices) * myscalar;
    }
    
    void moveTo(const pointing& newR) override
    {
        mytensor.moveTo(newR);
    }
    
};

//Sum only for tensorFamily
template<typename left, typename right , typename std::enable_if_t<std::is_base_of<tensorFamily,left>::value && std::is_base_of<tensorFamily,right>::value>* = nullptr>
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
* 
*/


//Save all linear combinations of minkTensorIntegrands in one class
struct minkTensorStack : tensorFamily
{
    std::vector<pointing> ns{};
    std::vector<double> weights{};
    
    minkTensorStack(const minkTensorStack& left, const minkTensorStack& right) : tensorFamily(left.rankA, left.rankB, left.curvIndex, left.r), ns(left.ns), weights(left.weights)
    {
        appendStack(right);
    }
    
    minkTensorStack(uint rank1, uint rank2, uint curvInd, const pointing& rNew) : tensorFamily(rank1, rank2, curvInd, rNew)
    {}
    
    explicit minkTensorStack(const minkTensorIntegrand& inp, double weight=1) : tensorFamily(inp.rankA, inp.rankB, inp.curvIndex, inp.r), ns{inp.n}, weights{weight}
    {}
    
    //Move/copy constructors default
    minkTensorStack(minkTensorStack&& other) = default;
    minkTensorStack(const minkTensorStack& other) = default;
    
    //Move/copy assignments leave rank untouched
    minkTensorStack& operator= (const minkTensorStack& other)
    {
        assert(rankA==other.rankA && rankB==other.rankB && curvIndex==other.curvIndex && "Trying to copy assign minkTensorStacks of different rank!");
        ns = other.ns;
        weights = other.weights;
        r = other.r;
        return *this;
    }
    
    minkTensorStack& operator= (minkTensorStack&& other)
    {
        assert(rankA==other.rankA && rankB==other.rankB && curvIndex==other.curvIndex && "Trying to move assign minkTensorStacks of different rank!");
        ns = std::move(other.ns);
        weights = std::move(other.weights);
        r = other.r;
        return *this;
    }
    
    ~minkTensorStack() = default;
    
    double accessElement(const std::vector<uint_fast8_t>& indices) const override
    {
        assert(ns.size()==weights.size() && "Error: number of weights different from number of normal vectors!");
        
        double retval = 0.;
        for(uint i=0; i<ns.size(); ++i)
        {
            minkTensorIntegrand tensorHere(rankA, rankB, curvIndex, r, ns[i]);
            retval+= tensorHere.accessElement(indices)*weights[i];
        }
        return retval;
    }
    
    void moveTo(const pointing& newR) override
    {
        for(auto n : ns)
        {
            n = parallelTransport(r, newR, n);
        }
        r = newR;
    }
    
    void addTensor(pointing n, double weight)
    {
        ns.push_back(n);
        weights.push_back(weight);
    }
    
    void addMinkTensorIntegrand(const minkTensorIntegrand& tens, double weight=1)
    {
        pointing newn = tens.n;
        if(arclength(r,tens.r)>1e-12)   newn = parallelTransport(tens.r,r,newn);
        ns.push_back(newn);
        weights.push_back(weight);
    }
    
    void appendStack(minkTensorStack other)
    {
        if(arclength(r,other.r)>1e-12)   other.moveTo(r);
        
        ns.reserve(ns.size()+other.ns.size());
        weights.reserve(weights.size()+other.weights.size());
        
        ns.insert(ns.end(), other.ns.begin(), other.ns.end());
        weights.insert(weights.end(), other.weights.begin(), other.weights.end());
    }
    
    explicit operator double() const
    {
        if(rankA+rankB == 0)
        {
            return accessElement({});
        }
        else
        {
            std::cerr << "Error: trying to convert non-rank zero tensor to double. Rank is " << rankA+rankB << std::endl;
            throw std::invalid_argument( "minkTensorStack: Rank too high for double conversion" );
        }
    }
    
    minkTensorStack& operator+= (const minkTensorStack& other)
    {
        appendStack(other);
        return *this;
    }
    
    minkTensorStack& operator+= (const minkTensorIntegrand& other)
    {
        addMinkTensorIntegrand(other);
        return *this;
    }
    
    minkTensorStack& operator*= (double other)
    {
        std::for_each(weights.begin(), weights.end(), [&other](double& inp){inp*=other;});
        return *this;
    }
};

minkTensorStack operator + (const minkTensorStack& lhs, const minkTensorStack& rhs)
{
    minkTensorStack returnval (lhs, rhs);
    return returnval;
}

minkTensorStack operator + (const minkTensorIntegrand& lhs, const minkTensorStack& rhs)
{
    minkTensorStack returnval (rhs);
    returnval.addMinkTensorIntegrand(lhs);
    return returnval;
}

minkTensorStack operator + (const minkTensorStack& lhs, const minkTensorIntegrand& rhs)
{
    minkTensorStack returnval (lhs);
    returnval.addMinkTensorIntegrand(rhs);
    return returnval;
}

minkTensorStack operator + (const minkTensorIntegrand& lhs, const minkTensorIntegrand& rhs)
{
    minkTensorStack returnval (lhs);
    returnval.addMinkTensorIntegrand(rhs);
    return returnval;
}

template<typename left , typename std::enable_if_t<std::is_base_of<minkTensorStack,left>::value || std::is_base_of<minkTensorIntegrand,left>::value>* = nullptr>
minkTensorStack operator* (const left& lhs, double rhs)
{
    minkTensorStack retStack(lhs);
    retStack*=rhs;
    return retStack;
}

template<typename right , typename std::enable_if_t<std::is_base_of<minkTensorStack,right>::value || std::is_base_of<minkTensorIntegrand,right>::value>* = nullptr>
minkTensorStack operator* (double lhs, const right& rhs)
{
    minkTensorStack retStack(rhs);
    retStack*=lhs;
    return retStack;
}




/***** Functions that work on tensors ****/

double trace(const tensorFamily& input) //sum of eigenvalues
{
    double sinT = sin(input.r.theta);
    if(input.rankA+input.rankB == 0) return input.accessElement({});
    
    std::vector<uint_fast8_t> indices(input.rankA+input.rankB,0);
    double summand = input.accessElement(indices); //zeroes
    
    for(uint i=0; i<indices.size(); i++) { indices.at(i) = 1; }
    //summand += input.accessElement(indices); //EDIT
    summand += input.accessElement(indices)*pow(sinT, input.rankA+input.rankB); //ones with metric contribution
    
    return summand;
}

double eigenValueQuotient(const tensorFamily& input) //TODO check
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
