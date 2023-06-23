#ifndef tensorOperations
#define tensorOperations

#include <cmath>
#include <algorithm>
#include <iostream>
#include <type_traits>
#include <vector>
#include <stdexcept>
#include <list>

#include "eigen/Eigen/Eigenvalues"
#include "eigen/Eigen/Dense"
#include "minkTensorIntegrand.hpp"
/** \file tensorOperations.hpp
 * \brief minkTensorStack, trace, eigenvalue quotient: all operations involving tensors
 */

///Save linear combinations of minkTensorIntegrands in one class
struct minkTensorStack
{
    const uint rankA{0}, rankB{0};
    const uint curvIndex{0};
    pointing r{1.5701963268,0};
    uint numnan{0}; ///< tracking how many contributions from masked pixels this tensor contains (in addition to content in nweights)
    uint numnull{0}; ///< tracking how many contributions from empty windows this tensor contains (in addition to content in nweights)
    
    std::list<std::pair<pointing,double>> nweights{}; ///< list normal vectors from which minkTensorIntegrands should be generated, and their respective weights
    
    
    minkTensorStack(minkTensorStack left, minkTensorStack right) : rankA(left.rankA), rankB(left.rankB), curvIndex(left.curvIndex), r(left.r), numnan(left.numnan), numnull(left.numnull), nweights(std::move(left.nweights))
    {
        numnan += right.numnan;
        numnull += right.numnull;
        appendStack_rr(std::move(right));
    }
    
    minkTensorStack(uint rank1, uint rank2, uint curvInd, const pointing& rNew) : rankA(rank1), rankB(rank2), curvIndex(curvInd), r(rNew)
    {
    }
    
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
            nweights.emplace_back(std::make_pair(inp.n,weight));
        }
    }
    
    //Move/copy constructors default
    minkTensorStack(minkTensorStack&& other) = default;
    minkTensorStack(const minkTensorStack& other) = default;
    ~minkTensorStack() = default;
    
    ///Move/copy assignments leave rank untouched (checked with assert)
    minkTensorStack& operator= (const minkTensorStack& other)
    {
        assert(rankA==other.rankA && rankB==other.rankB && curvIndex==other.curvIndex && "Trying to copy assign minkTensorStacks of different rank!");
        nweights = other.nweights;
        r = other.r;
        numnan = other.numnan;
        numnull = other.numnull;
        return *this;
    }
    
    ///Move/copy assignments leave rank untouched (checked with assert)
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
     * Checks whether tensor contains contours/anything non-empty. Empty spots are usually set to NAN in output map
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
        if(isEmpty()) {return 0.;}
        
        double retval = 0.;
        for(const auto& element : nweights)
        {
            minkTensorIntegrand tensorHere(rankA, rankB, curvIndex, r, std::get<0>(element));
            retval+= tensorHere.accessElement(indices)*std::get<1>(element);
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
        for(const auto& element : nweights)
        {
            minkTensorIntegrand tensorHere(rankA, rankB, curvIndex, r, std::get<0>(element));
            retval+= tensorHere.accessElement(indices)*std::get<1>(element);
        }
        return retval*factor;
    }
    
    /** Parallel transport all normal vectors in stack to newR along geodesic
     */
    void moveTo(const pointing& newR)
    {
        for(std::pair<pointing,double>& element : nweights)
        {
            std::get<0>(element) = parallelTransport(r, newR, std::get<0>(element));
        }
        r = newR;
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
            nweights.emplace_back(std::make_pair(newn,weight));
        }
    }
    
    /** Add normal vectors and weights of other stack with same ranks to this stack. Normal vectors are parallel transported to position of this stack if positions differ. Rvalue reference version.
     */
    void appendStack_rr(minkTensorStack&& other)
    {
        assert(rankA==other.rankA && rankB==other.rankB && curvIndex==other.curvIndex && "Trying to append minkTensorStacks of different rank!");
        numnan += other.numnan;
        numnull += other.numnull;
        
        if(arclength(r,other.r)>1e-12)   other.moveTo(r);
        
        auto itEnd = nweights.end();
        nweights.splice(itEnd,std::move(other.nweights));
        
    }
    
    /** Add normal vectors and weights of other stack with same ranks to this stack. Normal vectors are parallel transported to position of this stack if positions differ. Copy version.
     */
    void appendStack(minkTensorStack other)
    {
        appendStack_rr(std::move(other));
        
    }
    
    
    minkTensorStack& operator+= (const minkTensorStack& other)
    {
        appendStack(other);
        return *this;
    }
    
    minkTensorStack& operator+= (minkTensorStack&& other) //rvalue implementation because it can be used by normalHealpixInterface::toHealpix
    {
        appendStack_rr(std::move(other));
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


///minkTensor concept for either minkTensorStack or minkTensorIntegrand
template<typename T>
concept minkTensor = std::is_base_of_v<minkTensorStack,std::remove_cv_t<T>> || std::is_base_of_v<minkTensorIntegrand,std::remove_cv_t<T>>;

template<minkTensor left>
minkTensorStack operator+ (const left& lhs, const minkTensorIntegrand& rhs)
{
    minkTensorStack returnval (lhs);
    returnval.addMinkTensorIntegrand(rhs);
    return returnval;
}

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


///Empty minkTensorStack that represents a spot set to NAN
minkTensorStack nanTensor(uint rank1, uint rank2, uint curvInd, const pointing& rNew)
{
    minkTensorStack thistensor(rank1,rank2,curvInd,rNew);
    thistensor.numnan=1;
    return thistensor;
}

///Empty minkTensorStack that represents a spot without structure
minkTensorStack nullTensor(uint rank1, uint rank2, uint curvInd, const pointing& rNew)
{
    minkTensorStack thistensor(rank1,rank2,curvInd,rNew);
    thistensor.numnull=1;
    return thistensor;
}


/***** Functions that work on tensors ****/

///Calculates trace of a tensor (adding elemets with all indices zero and all one)
template<minkTensor tens>
double trace(const tens& input) //sum of eigenvalues
{
    if(input.isMasked()) {return NAN;}
    if(input.isEmpty()) {return 0.;}
    
    double sinT = sin(input.r.theta);
    if(input.rankA+input.rankB == 0) return input.accessElement_rescaled({});
    
    std::vector<uint_fast8_t> indices(input.rankA+input.rankB,0);
    double summand = input.accessElement_rescaled(indices); //zeroes
    
    for(uint i=0; i<indices.size(); i++) { indices[i] = 1; }
    summand += input.accessElement_rescaled(indices)*pow(sinT, input.rankA+input.rankB); //ones with metric contribution
    
    return summand;
}

///Turn rank 4 tensor into 3x3 matrix analog to description by Mehrabadi et al. (1990), can be used to find eigentensors and eigenvalues. Give sin^2(theta) of input to avoid calculating it twice
template<minkTensor tens>
Eigen::Matrix3d getMehrabadimatrix(const tens& input, double sin2T)
{
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
    return mehrabadimatrix;
}

///Calculates ratio of eigenvalues of tensor (rank 2) / norm of vector of eigenvalues when rank 4 tensor is reduced to 3x3-matrix
template<minkTensor tens>
double eigenValueQuotient(const tens& input)
{
    if(input.isMasked()) {return NAN;}
    if(input.isEmpty()) {return NAN;}
    
    uint ranksum = input.rankA+input.rankB;
    
    switch(ranksum){
        case 0:
            return input.accessElement({});
        case 1:
        {
            double retval = sqrt( pow(input.accessElement({0}),2) + pow(sin(input.r.theta)*input.accessElement({1}),2) ); //Take length of sum of normal vectors
            return retval;
        }
        case 2:
        {
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
            std::cout << "selfmade evq "<< retval << " eigen evq " << ratio <<"\n";
    #endif
            return ratio;
        }
        case 4:
        {
            double sin2T = sin(input.r.theta)*sin(input.r.theta);
            Eigen::Matrix3d mehrabadimatrix = getMehrabadimatrix(input, sin2T);
            Eigen::EigenSolver<Eigen::Matrix3d> solver(mehrabadimatrix,false);
            Eigen::Matrix<double,1,3> EVvec = solver.eigenvalues().real();
            std::sort(EVvec.begin(),EVvec.end()); //smallest to largest
            //Eigen::Vector3d reduced(std::abs(EVvec(0))-0., std::abs(EVvec(1))-0.,std::abs(EVvec(2))-0.);
            //double retval = std::abs(reduced.norm());
            double minmaxratio = std::abs(EVvec(1)); //0 largest 1 mid 2 smallest eigenvalue
            return minmaxratio;
            
            //double ratio = std::abs(realVec.maxCoeff()/realVec.minCoeff());
        }
        default:
            std::cerr << "Error: Eigenvalue quotient not implemented for rank 3 and higher than 4! Trying to calculate rankA rankB = " << input.rankA <<" "<< input.rankB << std::endl;
            throw std::invalid_argument("eigenValueQuotient not implemented for higher ranks");
    } //switch ranksum
}
///Returns normal vector pointing along largest EV of 2x2 Matrix
pointing angleOf2x2Mat(const Eigen::Matrix2d& zahlenfriedhof)
{
    pointing returnVec(0,1);
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
        returnVec.theta = EVec.col(0)[0];
        returnVec.phi = EVec.col(0)[1];
    }
    else
    {
        returnVec.theta = EVec.col(1)[0];
        returnVec.phi = EVec.col(1)[1];
    }
    return returnVec;
}

///Should return direction of eigenvector with highest eigenvalue. Zero means south, pi/2 means east
template<minkTensor tens>
double eigenVecDir(const tens& input)
{
    if(input.isMasked()) {return NAN;}
    if(input.isEmpty()) {return NAN;}
    
    uint ranksum = input.rankA+input.rankB;
    switch(ranksum){
        case 0:
        {
            std::cerr << "Error: Eigenvector direction not defined for rank 0!\n";
            throw std::invalid_argument("eigenVecDir not defined for rank 0");
        }
        case 1:
        {
            return giveAngle(pointing(input.accessElement({0}),input.accessElement({1})), input.r);
        }
        case 2:
        {
            //eigenvalues of matrix (a b, c d)
            double sin2T = sin(input.r.theta)*sin(input.r.theta); //pull down one index = sin^2 (theta) factor wherever left index = 1 (arbitrary choice)
    
            Eigen::Matrix2d zahlenfriedhof{
                {input.accessElement({0,0}),input.accessElement({0,1})*sin2T},
                {input.accessElement({0,1}),input.accessElement({1,1})*sin2T}
            };
            pointing relevantVec = angleOf2x2Mat(zahlenfriedhof);
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
        }
        case 4:
        {
            double sin2T = sin(input.r.theta)*sin(input.r.theta);
            Eigen::Matrix3d mehrabadimatrix = getMehrabadimatrix(input,sin2T);
            Eigen::EigenSolver<Eigen::Matrix3d> solver(mehrabadimatrix,true);
            auto EVvec = solver.eigenvalues().real();
            auto EVec = solver.eigenvectors().real();
            
            //select index of largest eigenvalue
            auto largestIndex = std::max_element(EVvec.begin(),EVvec.end()) - EVvec.begin();
            auto smallestIndex = std::min_element(EVvec.begin(),EVvec.end()) - EVvec.begin();
            //                     is one of the others zero ?          is one of the others one?
            auto midIndex = (largestIndex==0||smallestIndex==0) ? (largestIndex==1||smallestIndex==1 ? 2 : 1) : 0; // dritte zahl von 0,1,2 waehlen
            
            //Create Tensor belonging to EVec, pull one index down which adds factor of sin2T                    ?  --> above calculation gives Minktensor matrix with both indices on top
            Eigen::Matrix2d matrixtoEVec{
                {EVec.col(midIndex)[0],sin2T*EVec.col(midIndex)[1]},
                {EVec.col(midIndex)[1],sin2T*EVec.col(midIndex)[2]}
            };
            
            pointing relevantVec = angleOf2x2Mat(matrixtoEVec);
            double retangle = giveAngle(relevantVec,input.r);
            
            return retangle;
        }
        default:
            std::cerr << "Error: requesting eigenvector direction for ill-defined ranks: rankA rankB = " << input.rankA <<" "<< input.rankB << std::endl;
            throw std::invalid_argument("eigenVecDir not defined for this rank");
    } //switch ranksum
}





#endif
