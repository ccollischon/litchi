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
#include "eigen/Eigen/Dense"
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
    uint numnan{0}; ///< tracking how many contributions from a masked pixel this tensor contains (in addition to content in nweights)
    uint numnull{0}; ///< tracking how many contributions from empty windows this tensor contains (in addition to content in nweights)
    
    std::vector<std::pair<pointing,double>> nweights{}; ///< list normal Vectors from which minkTensorIntegrands should be generated and their respective weights
    
    minkTensorStack(const minkTensorStack& left, const minkTensorStack& right) : rankA(left.rankA), rankB(left.rankB), curvIndex(left.curvIndex), r(left.r), numnan(left.numnan), numnull(left.numnull), nweights(left.nweights)
    {
        numnan += right.numnan;
        numnull += right.numnull;
        appendStack(right);
    }
    
    minkTensorStack(uint rank1, uint rank2, uint curvInd, const pointing& rNew) : rankA(rank1), rankB(rank2), curvIndex(curvInd), r(rNew)
    {}
    
    explicit minkTensorStack(const minkTensorIntegrand& inp, double weight=1) : rankA(inp.rankA), rankB(inp.rankB), curvIndex(inp.curvIndex), r(inp.r), nweights{}
    {
        if(std::isnan(weight)) 
        {
            numnan = 1;
        }
        else if (std::abs(weight)<1e-15)
        {
            numnull = 1;
        }
        else
        {
            nweights.push_back(std::make_pair(inp.n,weight));
        }
    }
    
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
        numnan = other.numnan;
        numnull = other.numnull;
        return *this;
    }
    
    minkTensorStack& operator= (minkTensorStack&& other)
    {
        assert(rankA==other.rankA && rankB==other.rankB && curvIndex==other.curvIndex && "Trying to move assign minkTensorStacks of different rank!");
        nweights = std::move(other.nweights);
        r = std::move(other.r);
        numnan = other.numnan;
        numnull = other.numnull;
        return *this;
    }
    
    /**
     * Checks whether the amont of masked/NAN contributions to this stack is too high.
     * Allows for a relative amount of masked pixels up to 1/16
     */
    bool isMasked() const
    {
        return 15*numnan > nweights.size()+numnull;
    }
    
    /**
     * Checks whether tensor contains contours/anything non-empty. Empty spots are usually set to NAN ind output map
     */
    bool isEmpty() const
    {
        return nweights.size() == 0;
    }
    
    /**
     * Access element of linear combination of tensors determined by weights
     * \param indices Indices of desired element
     * \return Value of desired element
     */
    double accessElement(const std::vector<uint_fast8_t>& indices) const
    {
        if(isMasked()) {return NAN;}
        
        double retval = 0.;
        for(uint i=0; i<nweights.size(); ++i)
        {
            minkTensorIntegrand tensorHere(rankA, rankB, curvIndex, r, std::get<0>(nweights[i]));
            retval+= tensorHere.accessElement(indices)*std::get<1>(nweights[i]);
        }
        return retval;
    }
    
    /**
     * Access element of linear combination of tensors and rescales with total/masked # pixels
     * \param indices Indices of desired element
     * \return Value of desired element
     */
    double accessElement_rescaled(const std::vector<uint_fast8_t>& indices) const
    {
        //if(isMasked()) {return NAN;}
        double factor = 1.*(nweights.size()+numnull+numnan)/(1.*(nweights.size()+numnull));
        
        double retval = 0.;
        for(uint i=0; i<nweights.size(); ++i)
        {
            minkTensorIntegrand tensorHere(rankA, rankB, curvIndex, r, std::get<0>(nweights[i]));
            retval+= tensorHere.accessElement(indices)*std::get<1>(nweights[i]);
        }
        return retval*factor;
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
        if(std::isnan(weight)) {
            numnan += 1;
        }
        else if (std::abs(weight)<1e-15) {
            numnull += 1;
        }
        else {
            nweights.push_back(std::make_pair(n,weight));
        }
    }
    
    /** Add normal vector of minkTensorIntegrand with same ranks to stack. Normal vector is parallel transported to position of stack if positions differ
     */
    void addMinkTensorIntegrand(const minkTensorIntegrand& tens, double weight=1)
    {
        assert(rankA==tens.rankA && rankB==tens.rankB && curvIndex==tens.curvIndex && "Trying to addMinkTensorIntegrand to minkTensorStack of different rank!");
        if(std::isnan(weight))
        {
            numnan += 1;
        }
        else if (std::abs(weight)<1e-15) 
        {
            numnull += 1;
        }
        else
        {
            pointing newn = tens.n;
            if(arclength(r,tens.r)>1e-12)   newn = parallelTransport(tens.r,r,newn);
            nweights.push_back(std::make_pair(newn,weight));
        }
    }
    
    /** Add normal vectors and weights of other stack with same ranks to this stack. Normal vectors are parallel transported to position of this stack if positions differ
     */
    void appendStack(minkTensorStack other)
    {
        assert(rankA==other.rankA && rankB==other.rankB && curvIndex==other.curvIndex && "Trying to append minkTensorStacks of different rank!");
        numnan += other.numnan;
        numnull += other.numnull;
        
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
        if(std::isnan(other))
        {
            numnan += nweights.size();
            nweights.clear();
        }
        else if(std::abs(other)<1e-15)
        {
            numnull += nweights.size();
            nweights.clear();
        }
        else
        {
            std::for_each(nweights.begin(), nweights.end(), [&other](auto& inp){std::get<1>(inp)*=other;});
        }
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

minkTensorStack nanTensor(uint rank1, uint rank2, uint curvInd, const pointing& rNew)
{
    minkTensorStack thistensor(rank1,rank2,curvInd,rNew);
    thistensor.numnan=1;
    return thistensor;
}

minkTensorStack nullTensor(uint rank1, uint rank2, uint curvInd, const pointing& rNew)
{
    minkTensorStack thistensor(rank1,rank2,curvInd,rNew);
    thistensor.numnull=1;
    return thistensor;
}


/***** Functions that work on tensors ****/

template<typename tens>
double trace(const tens& input) //sum of eigenvalues
{
    //if(input.isMasked()) {return NAN;} //setting pixels to nan is not necessary here
    
    double sinT = sin(input.r.theta);
    if(input.rankA+input.rankB == 0) return input.accessElement_rescaled({});
    
    std::vector<uint_fast8_t> indices(input.rankA+input.rankB,0);
    double summand = input.accessElement_rescaled(indices); //zeroes
    
    for(uint i=0; i<indices.size(); i++) { indices.at(i) = 1; }
    summand += input.accessElement_rescaled(indices)*pow(sinT, input.rankA+input.rankB); //ones with metric contribution
    
    return summand;
}

template<typename tens>
double eigenValueQuotient(const tens& input)
{
    if(input.isMasked()) {return NAN;}
    if(input.isEmpty()) {return NAN;}
    
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
        
        Eigen::Matrix2d zahlenfriedhof{
            {input.accessElement({0,0}),input.accessElement({0,1})*sin2T},
            {input.accessElement({0,1}),input.accessElement({1,1})*sin2T}
        };
        Eigen::EigenSolver<Eigen::Matrix2d> solver(zahlenfriedhof,false);
        auto EVvec = solver.eigenvalues().real();
        double ratio = std::abs(EVvec.maxCoeff()/EVvec.minCoeff());
        
        // Check if the above calculation does the same thing as manual calculation
#ifdef THISRUNSINATEST
        double dplusa = input.accessElement({1,1})*sin2T+input.accessElement({0,0});
        double adminusbc = input.accessElement({0,0})*input.accessElement({1,1})*sin2T - pow(input.accessElement({0,1}),2)*sin2T; //tensors are symmetric here, one of them needs factor
        double twolambda1 = dplusa + sqrt(dplusa*dplusa - 4*adminusbc);
        double twolambda2 = dplusa - sqrt(dplusa*dplusa - 4*adminusbc);
        double retval;
        if(twolambda1>twolambda2) retval = twolambda1/twolambda2;
        else retval = twolambda2/twolambda1;
        //return std::isnan(retval) ? 0 : retval; //return 0 instead of nan in case of division by zero
        std::cout << "selfmade evq "<< retval << " eigen evq " << ratio <<"\n";
#endif
        return ratio;
        
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
        Eigen::Vector3d reduced(std::abs(EVvec(0))-0., std::abs(EVvec(1))-0.,std::abs(EVvec(2))-0.);
        double retval = std::abs(reduced.norm());
        return retval;
        
        //double ratio = std::abs(realVec.maxCoeff()/realVec.minCoeff());
        
        //return std::isnan(ratio) ? 0 : ratio;
    } else
    {
        std::cerr << "Error: Eigenvalue quotient not implemented for rank 3 and higher than 4! Trying to calculate rankA rankB = " << input.rankA <<" "<< input.rankB << std::endl;
        throw std::invalid_argument("eigenValueQuotient not implemented for higher ranks");
    }
}

//Should return direction of eigenvector with highest eigenvalue. Zero means south, pi/2 means east
double eigenVecDir(const auto& input)
{
    if(input.isMasked()) {return NAN;}
    if(input.isEmpty()) {return NAN;}
    
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

        pointing relevantVec(0,1); //to be filled
        
        Eigen::Matrix2d zahlenfriedhof{
            {input.accessElement({0,0}),input.accessElement({0,1})*sin2T},
            {input.accessElement({0,1}),input.accessElement({1,1})*sin2T}
        };
        Eigen::EigenSolver<Eigen::Matrix2d> solver(zahlenfriedhof,true);
        auto EVvec = solver.eigenvalues().real();
        auto EVec = solver.eigenvectors().real();
        if (solver.info() != Eigen::Success)
        {
            std::cerr << "Error: Could not find eigenvalues/vectors of the following matrix:\n";
            std::cerr << zahlenfriedhof << std::endl;
        }
        
        if(EVvec[0]>EVvec[1]) // grab eigenvector with largest eigenvalue
        {
            relevantVec.theta = EVec.col(0)[0];
            relevantVec.phi = EVec.col(0)[1];
        }
        else
        {
            relevantVec.theta = EVec.col(1)[0];
            relevantVec.phi = EVec.col(1)[1];
        }
        double retangle = giveAngle(relevantVec,input.r);
        
#ifdef THISRUNSINATEST
        double dplusa = input.accessElement({1,1})*sin2T+input.accessElement({0,0});
        double adminusbc = input.accessElement({0,0})*input.accessElement({1,1})*sin2T - pow(input.accessElement({0,1}),2)*sin2T; //tensors are symmetric here, one of them needs factor
        double twolambda1 = dplusa + sqrt(dplusa*dplusa - 4*adminusbc);
        double twolambda2 = dplusa - sqrt(dplusa*dplusa - 4*adminusbc);
        
        pointing oldVec(0,1);
        if(twolambda1>twolambda2) // grab eigenvector with largest eigenvalue
        {
            double aminusd = input.accessElement({0,0}) - input.accessElement({1,1})*sin2T;
            oldVec.theta = (aminusd + sqrt(dplusa*dplusa-4*adminusbc)) / (2*input.accessElement({0,1})*sin2T);
        }
        else
        {
            double aminusd = input.accessElement({0,0}) - input.accessElement({1,1})*sin2T;
            oldVec.theta = (aminusd - sqrt(dplusa*dplusa-4*adminusbc)) / (2*input.accessElement({0,1})*sin2T);
        }
        double oldangle = giveAngle(oldVec,input.r);
        std::cout << "selfmade evd "<< oldangle << " eigen evd " << (retangle<0 ? retangle+3.14159 : retangle ) <<"\n";
#endif
        return retangle<0 ? retangle+3.14159 : retangle; //add pi to get consistent values between 0 and pi, only interested in direction
        
    } else if(ranksum == 4 && false)
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
        Eigen::EigenSolver<Eigen::Matrix3d> solver(mehrabadimatrix,true);
        auto EVvec = solver.eigenvalues().real();
        auto EVec = solver.eigenvectors().real();
        
        //TODO select largest/smallest eigenvalue
        
        //Create Tensor belonging to one of the EVecs, pull one index down which adds factor of sin2T
        Eigen::Matrix2d matrixtoEVec{
            {EVec(0,0),sin2T*EVec(0,1)},
            {EVec(0,1),sin2T*EVec(0,2)}
        };
        
        return 0;
        
        //double ratio = std::abs(realVec.maxCoeff()/realVec.minCoeff());
        
        //return std::isnan(ratio) ? 0 : ratio;
    } else
    {
        std::cerr << "Error: requesting eigenvector direction for ill-defined ranks: rankA rankB = " << input.rankA <<" "<< input.rankB << std::endl;
        throw std::invalid_argument("eigenVecDir not defined for this rank");
    }
}





#endif
