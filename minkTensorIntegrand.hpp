#ifndef minkTensorInt
#define minkTensorInt

//#include "tensorFamily.hpp"
#include "geometryHelpers.hpp"
#include <algorithm>
#include <stdexcept>
//#include <iostream>
#include <numeric>
#include <cassert>
#include <vector>
#include <cmath>
#include "healpix_cxx/healpix_map.h"
/** \file minkTensorIntegrand.hpp
 * \brief minkTensorIntegrand struct contains symmetric tensor product part of Minkowski tensors
 */


///struct to represent the r^a x n^b part of minkowski tensor calculation, only calculates value when necessary
struct minkTensorIntegrand {
    const uint rankA{0}, rankB{0};
    const uint curvIndex{0};
    pointing r{1.5701963268,0};
    pointing n {1,0};
    
    minkTensorIntegrand(uint rank1, uint rank2, uint curvInd = 0) : rankA(rank1), rankB(rank2), curvIndex(curvInd), n(0,1)
    { //simple constructor for empty tensor, just give ranks
    }
    
    minkTensorIntegrand(uint rank1, uint rank2, uint curvInd, const pointing& rNew, const pointing& nNew) : rankA(rank1), rankB(rank2), curvIndex(curvInd), r(rNew), n{nNew}
    { //constructor with all necessary data
    }
    
    /**
     * Access element of this tensor as generated by symmetric product of r^rankA and n^rankB according to ranks
     * \param indices Indices of desired element
     * \return Value of desired element
     */
    double accessElement(const std::vector<uint_fast8_t>& inputindices) const // tensor element with indices i_1 ... i_{a+b} is sum over all permutations of indices of multiplication of elements of r and n at these indices
    {
        std::vector<uint_fast8_t> indices = inputindices; //Copy only in this class, reference in all others
        uint indicesSize = indices.size();
        
        assert((indicesSize == rankA+rankB) &&  ("Error: requesting element with wrong number of indices") );
        
        
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
        
        
        
        returnval /= numberOfPermutations;
        return returnval;
        
    }
    
    /// Parallel transport normal vector to newR along geodesic
    void moveTo(const pointing& newR)
    {
        if(arclength(r,newR)>1e-12)
        {
            n = parallelTransport(r, newR, n);
            r = newR;
        }
    }
    
    
};

#endif
