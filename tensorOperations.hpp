#ifndef tensorOperations
#define tensorOperations

//#include "tensorFamily.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <type_traits>
#include <vector>
#include <stdexcept>

#include "eigen/Eigen/Eigenvalues"
#include "minkTensorIntegrand.hpp"
/** \file tensorOperations.hpp
 * \brief minkTensorStack, trace, eigenvalue quotient: all operations involving tensors
 */

///Save all linear combinations of minkTensorIntegrands in one class
struct minkTensorStack
{
    const uint rankA{0}, rankB{0};
    const uint curvIndex{0};
    pointing r{1.5701963268,0};
    //std::vector<pointing> ns{}; ///< list of normal Vectors from which minkTensorIntegrands should be generated
    //std::vector<double> weights{}; ///< list of weights for minkTensorIntegrands
    std::vector<std::pair<pointing,double>> nweights{}; ///< list normal Vectors from which minkTensorIntegrands should be generated and their respective weights
    
    minkTensorStack(const minkTensorStack& left, const minkTensorStack& right) : rankA(left.rankA), rankB(left.rankB), curvIndex(left.curvIndex), r(left.r), nweights(left.nweights)
    {
        appendStack(right);
    }
    
    minkTensorStack(uint rank1, uint rank2, uint curvInd, const pointing& rNew) : rankA(rank1), rankB(rank2), curvIndex(curvInd), r(rNew)
    {}
    
    explicit minkTensorStack(const minkTensorIntegrand& inp, double weight=1) : rankA(inp.rankA), rankB(inp.rankB), curvIndex(inp.curvIndex), r(inp.r), nweights{std::make_pair(inp.n,weight)}
    {}
    
    //Move/copy constructors default
    minkTensorStack(minkTensorStack&& other) = default;
    minkTensorStack(const minkTensorStack& other) = default;
    ~minkTensorStack() = default;
    
    //Move/copy assignments leave rank untouched
    minkTensorStack& operator= (const minkTensorStack& other)
    {
        assert(rankA==other.rankA && rankB==other.rankB && curvIndex==other.curvIndex && "Trying to copy assign minkTensorStacks of different rank!");
        nweights = other.nweights;
        r = other.r;
        return *this;
    }
    
    minkTensorStack& operator= (minkTensorStack&& other)
    {
        assert(rankA==other.rankA && rankB==other.rankB && curvIndex==other.curvIndex && "Trying to move assign minkTensorStacks of different rank!");
        nweights = std::move(other.nweights);
        r = std::move(other.r);
        return *this;
    }
    
    /**
     * Access element of linear combination of tensors determined by weights
     * \param indices Indices of desired element
     * \return Value of desired element
     */
    double accessElement(const std::vector<uint_fast8_t>& indices) const
    {
        double retval = 0.;
        for(uint i=0; i<nweights.size(); ++i)
        {
            minkTensorIntegrand tensorHere(rankA, rankB, curvIndex, r, std::get<0>(nweights[i]));
            retval+= tensorHere.accessElement(indices)*std::get<1>(nweights[i]);
        }
        return retval;
    }
    
    /**
     * Access element of linear combination of tensors without multiplying by weights
     * \param indices Indices of desired element
     * \return Value of desired element
     */
    double accessElement_noweights(const std::vector<uint_fast8_t>& indices) const
    {
        double retval = 0.;
        for(uint i=0; i<nweights.size(); ++i)
        {
            minkTensorIntegrand tensorHere(rankA, rankB, curvIndex, r, std::get<0>(nweights[i]));
            retval+= tensorHere.accessElement(indices);
        }
        return retval;
    }
    
    /** Parallel transport all normal vectors in stack to newR along geodesic
     */
    void moveTo(const pointing& newR)
    {
        for(uint i=0; i<nweights.size(); ++i)
        {
            std::get<0>(nweights[i]) = parallelTransport(r, newR, std::get<0>(nweights[i]));
        }
        r = newR;
    }
    
    /** Add Tensor with given pointing and weight to stack
     */
    void addTensor(pointing n, double weight)
    {
        nweights.push_back(std::make_pair(n,weight));
    }
    
    /** Add normal vector of minkTensorIntegrand with same ranks to stack. Normal vector is parallel transported to position of stack if positions differ
     */
    void addMinkTensorIntegrand(const minkTensorIntegrand& tens, double weight=1)
    {
        assert(rankA==tens.rankA && rankB==tens.rankB && curvIndex==tens.curvIndex && "Trying to addMinkTensorIntegrand to minkTensorStack of different rank!");
        pointing newn = tens.n;
        if(arclength(r,tens.r)>1e-12)   newn = parallelTransport(tens.r,r,newn);
        nweights.push_back(std::make_pair(newn,weight));
    }
    
    /** Add normal vectors and weights of other stack with same ranks to this stack. Normal vectors are parallel transported to position of this stack if positions differ
     */
    void appendStack(minkTensorStack other)
    {
        assert(rankA==other.rankA && rankB==other.rankB && curvIndex==other.curvIndex && "Trying to append minkTensorStacks of different rank!");
        if(arclength(r,other.r)>1e-12)   other.moveTo(r);
        
        nweights.reserve(nweights.size()+other.nweights.size());
        
        nweights.insert(nweights.end(), other.nweights.begin(), other.nweights.end());
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
        std::for_each(nweights.begin(), nweights.end(), [&other](auto& inp){std::get<1>(inp)*=other;});
        return *this;
    }
};

minkTensorStack operator+ (const minkTensorStack& lhs, const minkTensorStack& rhs)
{
    minkTensorStack returnval (lhs, rhs);
    return returnval;
}

minkTensorStack operator+ (const minkTensorIntegrand& lhs, const minkTensorStack& rhs)
{
    minkTensorStack returnval (rhs);
    returnval.addMinkTensorIntegrand(lhs);
    return returnval;
}

minkTensorStack operator+ (const minkTensorStack& lhs, const minkTensorIntegrand& rhs)
{
    minkTensorStack returnval (lhs);
    returnval.addMinkTensorIntegrand(rhs);
    return returnval;
}

minkTensorStack operator+ (const minkTensorIntegrand& lhs, const minkTensorIntegrand& rhs)
{
    minkTensorStack returnval (lhs);
    returnval.addMinkTensorIntegrand(rhs);
    return returnval;
}

template<typename T>
concept minkTensor = std::is_base_of<minkTensorStack,T>::value || std::is_base_of<minkTensorIntegrand,T>::value;

template<minkTensor left>
minkTensorStack operator* (const left& lhs, double rhs)
{
    minkTensorStack retStack(lhs);
    retStack*=rhs;
    return retStack;
}

template<minkTensor right>
minkTensorStack operator* (double lhs, const right& rhs)
{
    minkTensorStack retStack(rhs);
    retStack*=lhs;
    return retStack;
}




/***** Functions that work on tensors ****/

template<typename tens>
double trace(const tens& input) //sum of eigenvalues
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

template<typename tens>
double eigenValueQuotient(const tens& input)
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
        double retval;
        if(twolambda1>twolambda2) retval = twolambda1/twolambda2;
        else retval = twolambda2/twolambda1;
        return std::isnan(retval) ? 0 : retval; //return 0 instead of nan in case of division by zero
        
    } else if(ranksum == 4)
    {
        double sin2T = sin(input.r.theta)*sin(input.r.theta);
        //Only need these 5 because of total symmetry:
        double W0000 = input.accessElement({0,0,0,0});
        double W1100 = input.accessElement({1,1,0,0});
        double W1111 = input.accessElement({1,1,1,1});
        double W0001 = input.accessElement({0,0,0,1});
        double W1101 = input.accessElement({1,1,0,1});
        //Create matrix, calculate eigenvalues
        Eigen::Matrix3d mehrabadimatrix{
            {W0000,2.*sin2T*W0001,sin2T*sin2T*W1100},
            {W0001,2.*sin2T*W1100,sin2T*sin2T*W1101},
            {W1100,2.*sin2T*W1101,sin2T*sin2T*W1111}
        };
        Eigen::EigenSolver<Eigen::Matrix3d> solver(mehrabadimatrix,false);
        auto EVvec = solver.eigenvalues();
        Eigen::Vector3d realVec(std::abs(EVvec(0)), std::abs(EVvec(1)),std::abs(EVvec(2)));
        double ratio = std::abs(realVec.maxCoeff()/realVec.minCoeff());
        
        return std::isnan(ratio) ? 0 : ratio;
    } else
    {
        std::cerr << "Error: Eigenvalue quotient not implemented for rank higher than 2! Trying to calculate rankA rankB = " << input.rankA <<" "<< input.rankB << std::endl;
        throw std::invalid_argument("eigenValueQuotient not implemented for higher ranks");
    }
}

//Should return direction of eigenvector with highest eigenvalue. Zero means south, pi/2 means east
double eigenVecDir(const auto& input)
{
    uint ranksum = input.rankA+input.rankB;
    if (ranksum == 0)
    {
        std::cerr << "Error: Eigenvector direction not defined for rank 0!\n";
        throw std::invalid_argument("eigenVecDir not defined for rank 0");
    }
    else if (ranksum==1)
    {
        return giveAngle(pointing(input.accessElement({0}),input.accessElement({1})), input.r);
    }
    else if (ranksum==2)
    {  //eigenvalues of matrix (a b, c d)
        double sin2T = sin(input.r.theta)*sin(input.r.theta); //pull down one index = sin^2 (theta) factor wherever left index = 1 (arbitrary choice)
        double dplusa = input.accessElement({1,1})*sin2T+input.accessElement({0,0});
        double adminusbc = input.accessElement({0,0})*input.accessElement({1,1})*sin2T - pow(input.accessElement({0,1}),2)*sin2T; //tensors are symmetric here, one of them needs factor
        double twolambda1 = dplusa + sqrt(dplusa*dplusa - 4*adminusbc);
        double twolambda2 = dplusa - sqrt(dplusa*dplusa - 4*adminusbc);
        
        pointing relevantVec(0,1); //set phi=1 as free parameter here
        if(twolambda1>twolambda2)
        {
            double aminusd = input.accessElement({0,0}) - input.accessElement({1,1})*sin2T;
            relevantVec.theta = (aminusd + sqrt(dplusa*dplusa-4*adminusbc)) / (2*input.accessElement({0,1})*sin2T);
        }
        else
        {
            double aminusd = input.accessElement({0,0}) - input.accessElement({1,1})*sin2T;
            relevantVec.theta = (aminusd - sqrt(dplusa*dplusa-4*adminusbc)) / (2*input.accessElement({0,1})*sin2T);
        }
        return giveAngle(relevantVec,input.r);
    }
    else
    {
        std::cerr << "Error: requesting eigenvector direction for ill-defined ranks: rankA rankB = " << input.rankA <<" "<< input.rankB << std::endl;
        throw std::invalid_argument("eigenVecDir not defined for this rank");
    }
}





#endif
