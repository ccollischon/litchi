#ifndef minkTensorInt
#define minkTensorInt

#include "tensorFamily.hpp"
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <numeric>
#include <vector>
#include <cmath>
#include "healpix_cxx/healpix_map.h"

//struct to represent the r^a x n^b part of minkowski tensor calculation, only calculates value when necessary
struct minkTensorIntegrand : tensorFamily {
    pointing r, n;
    
    minkTensorIntegrand(uint rank1, uint rank2, uint curvInd = 0) : tensorFamily(rank1, rank2, curvInd), r(0,0), n(0,0)
    { //simple constructor for empty tensor, just give ranks
    }
    
    minkTensorIntegrand(uint rank1, uint rank2, uint curvInd, pointing rNew, pointing nNew) : tensorFamily(rank1, rank2, curvInd)
    { //constructor with all necessary data
        r = rNew;
        n = nNew;
    }
    
    double accessElement(std::vector<uint_fast8_t> indices) const // tensor element with indices i_1 ... i_{a+b} is sum over all permutations of indices of multiplication of elements of r and n at these indices
    {
        uint indicesSize = indices.size();
        if(indicesSize != rankA+rankB) //Correct number of indices
        {
            std::cerr << "Error: requesting element with wrong number of indices: " << indicesSize << " expected number: " << rankA+rankB << std::endl;
            throw std::invalid_argument( "minkTensor wrong rank" );
        }
        if(indicesSize>0){ //Here Tensor case only
            //Check if all indices are 0 or 1/no index larger than 1
            std::vector<uint_fast8_t>::iterator maxel = std::max_element(indices.begin(), indices.end());
            if (*maxel>1){
                std::cerr << "Error: given tensor index too high: " << (uint16_t) *maxel << ", expected 0 or 1" << std::endl;
                throw std::invalid_argument( "minkTensor index too high" );
            }
        }
        else return 1; //Scalar case, just return 1 and deal with tensors afterwards
        
        
        //// Calculate actual element
        double returnval = 0; 
        //uint numberOfOnes = std::accumulate(indices.begin(), indices.end(), 0);
        uint numberOfPermutations = 0;
        
        std::sort (indices.begin(), indices.end());        
        //loop over all permutations (non-degenerate by next_permutation), total number over num of 1
        do //permutations of indices
        {
            //for(unsigned char index : indices) std::cout << (int)index;
            //std::cout << std::endl;
            
            double summand = 1;
            uint numberOfROnes = std::accumulate(indices.begin(), indices.begin()+rankA, 0); //Number of ones among indices of r / n
            uint numberOfNOnes = std::accumulate(indices.begin()+rankA, indices.end(), 0);
            
            summand *= pow(r.theta, numberOfROnes); //different r contributions, depending on index
            summand *= pow(r.phi, rankA-numberOfROnes); // Might need Jacobian, because Polar edit: no, sin is in metric
            
            summand *= pow(n.theta, numberOfNOnes); //same for n
            summand *= pow(n.phi, rankB-numberOfNOnes);
            
            returnval += summand;
            numberOfPermutations++;
        } while ( std::next_permutation(indices.begin(), indices.end()) );
        
        
        //Calculate binomialCoeff, not factorial because next_permutation  does not create degeneracies
        
        
        returnval /= numberOfPermutations; //binomialCoeff(indicesSize, numberOfOnes);
        return returnval;
        
    }
    
    minkTensorIntegrand& operator= (const minkTensorIntegrand& other)
    {
        if( (rankA != other.rankA) || (rankB != other.rankB) || (curvIndex != other.curvIndex))
        {
            std::cerr << "Error: trying to use operator = on tensors of different type. Left has (rankA, rankB, curvIndex) =  (" << rankA << "," << rankB << "," << curvIndex << "), Right has ("  << other.rankA << "," << other.rankB << "," << other.curvIndex << "), this is not intended, ranks are const" << std::endl;
            throw std::invalid_argument( "minkTensorIntegrand: Different ranks in operator =" );
        }
        r = other.r;
        n = other.n;
        return *this;
    }
    
};

#endif
