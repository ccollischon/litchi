#ifndef minkTensorInt
#define minkTensorInt

#include "tensorFamily.hpp"
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <numeric>
#include <cassert>
#include <vector>
#include <cmath>
#include "healpix_cxx/healpix_map.h"

//struct to represent the r^a x n^b part of minkowski tensor calculation, only calculates value when necessary
struct minkTensorIntegrand : tensorFamily {
    pointing n {0,1};
    
    minkTensorIntegrand(uint rank1, uint rank2, uint curvInd = 0) : tensorFamily(rank1, rank2, curvInd), n(0,1)
    { //simple constructor for empty tensor, just give ranks
    }
    
    minkTensorIntegrand(uint rank1, uint rank2, uint curvInd, const pointing& rNew, const pointing& nNew) : tensorFamily(rank1, rank2, curvInd, rNew), n{nNew}
    { //constructor with all necessary data
    }
    
    double accessElement(const std::vector<uint_fast8_t>& inputindices) const override // tensor element with indices i_1 ... i_{a+b} is sum over all permutations of indices of multiplication of elements of r and n at these indices
    {
        std::vector<uint_fast8_t> indices = inputindices; //Copy only in this class, reference in all others
        uint indicesSize = indices.size();
        
        assert((indices.size() == rankA+rankB) &&  ("Error: requesting element with wrong number of indices") );
        
        
        if(indicesSize>0){ //Here Tensor case only
            //Check if all indices are 0 or 1/no index larger than 1
            assert((*(std::max_element(indices.begin(), indices.end()))<2) && ("Error: given tensor index too high, expected 0 or 1") );
        }
        else return 1; //Scalar case, just return 1 and deal with tensors afterwards
        
        
        //// Calculate actual element
        double returnval = 0; 
        
        uint numberOfPermutations = 0;
        
        std::sort (indices.begin(), indices.end());        
        //loop over all permutations (non-degenerate by next_permutation), total number over num of 1
        do //permutations of indices
        {
            double summand = 1;
            uint numberOfROnes = std::accumulate(indices.begin(), indices.begin()+rankA, 0); //Number of ones among indices of r / n
            uint numberOfNOnes = std::accumulate(indices.begin()+rankA, indices.end(), 0);
            
            summand *= pow(r.theta, rankA-numberOfROnes); //different r contributions, depending on index. theta for zeros, phi for ones
            summand *= pow(r.phi, numberOfROnes); // Might need Jacobian, because Polar edit: no, sin is in metric
            
            summand *= pow(n.theta, rankB-numberOfNOnes); //same for n
            summand *= pow(n.phi, numberOfNOnes);
            
            returnval += summand;
            numberOfPermutations++;
        } while ( std::next_permutation(indices.begin(), indices.end()) );
        
        
        //Calculate binomialCoeff, not factorial because next_permutation  does not create degeneracies
        
        
        returnval /= numberOfPermutations; //binomialCoeff(indicesSize, numberOfOnes);
        return returnval;
        
    }
    
    //minkTensorIntegrand& operator= (const minkTensorIntegrand& other)
    //{
        //if( (rankA != other.rankA) || (rankB != other.rankB) || (curvIndex != other.curvIndex))
        //{
            //std::cerr << "Error: trying to use operator = on tensors of different type. Left has (rankA, rankB, curvIndex) =  (" << rankA << "," << rankB << "," << curvIndex << "), Right has ("  << other.rankA << "," << other.rankB << "," << other.curvIndex << "), this is not intended, ranks are const" << std::endl;
            //throw std::invalid_argument( "minkTensorIntegrand: Different ranks in operator =" );
        //}
        //r = other.r;
        //n = other.n;
        //return *this;
    //}
    
};

#endif
