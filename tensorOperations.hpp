#ifndef tensorOperations
#define tensorOperations

#include "minkTensorStack.hpp"
#include "minkTensorIntegrand.hpp"
#include "irreducibleMinkTens.hpp"

#include "eigen/Eigen/Eigenvalues"
#include "eigen/Eigen/Dense"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <type_traits>
#include <vector>
#include <stdexcept>
#include <list>

/** \file tensorOperations.hpp
 * \brief Functions to compute anisotropy, trace, and preferred direction of a tensor
 */

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

///Calculates the corresponding irreducible tensor, using input.rankB as l
std::complex<double> getPsilm(int m, const minkTensorStack& input)
{
    std::complex<double> psilm = 0.;
    for(const auto& element : input.nweights)
    {
        irreducibleMinkTens tensorHere((int)input.rankB, m, input.r, std::get<0>(element));
        psilm += tensorHere.accessElement()*std::get<1>(element);
    }
    return psilm;
}

///Calls eigenValueQuotient(input) to get a cartesian anisotropy measure
double anisotropy_cart(const minkTensorStack& input)
{
    double aniso = eigenValueQuotient<minkTensorStack>(input);
    return aniso;
}

///Calculates anisotropy by interpreting input as stack of irreducible tensors
double anisotropy_irr(const minkTensorStack& input)
{
    if(input.isMasked()) {return NAN;}
    if(input.isEmpty()) {return NAN;}
    
    std::complex<double> psilm = getPsilm(0, input);
    double retval = std::abs(psilm)*std::abs(psilm);
    
    return retval;
}

///Calls eigenVecDir(input) to get a cartesian direction measure
double direction_cart(const minkTensorStack& input)
{
    double dir = eigenVecDir<minkTensorStack>(input);
    return dir;
}

///Calculates direction by interpreting input as stack of irreducible tensors
double direction_irr(const minkTensorStack& input)
{
    if(input.isMasked()) {return NAN;}
    if(input.isEmpty()) {return NAN;}
    
    std::complex<double> psilm = getPsilm(0, input);
    double retval = std::arg(psilm)*std::arg(psilm);
    
    return retval;
}




#endif
