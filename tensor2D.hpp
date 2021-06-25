#ifndef tensor2d
#define tensor2d

#include <stdexcept>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <cassert>
#include "tensorFamily.hpp"
#include "geometryHelpers.hpp"


//struct to represent 2D (Minkowski) tensors of rankA + rankB explicitly
struct tensor2D : tensorFamily
{
    std::vector<double> content {};
    
    tensor2D (uint rank1, uint rank2, uint curvInd = 0) : tensorFamily(rank1, rank2, curvInd) , content(pow(2,rankA+rankB),0)
    { //simple constructor for empty tensor 
    }
    
    explicit tensor2D (const tensorFamily& other) : tensorFamily(other.rankA, other.rankB, other.curvIndex, other.r), content(pow(2,rankA+rankB),0)
    {
        assign(other);
    }
    
    explicit tensor2D (double number,uint curvIndex=0) : tensorFamily{0,0,curvIndex}, content{number}
    {}
    
    tensor2D (const tensor2D& other) = default;
    tensor2D (tensor2D&& other) = default;
    ~tensor2D() = default;
    
    //tensor2D (uint rank1, uint rank2, uint curvInd = 0) : tensorFamily(rank1, rank2, curvInd)
    //{ //simple constructor for empty tensor
        //std::vector<double> newContent(pow(2,rank1+rank2),0);
        //content = newContent;
    //}
    
    //explicit tensor2D (const tensorFamily& other) : tensorFamily(other.rankA, other.rankB, other.curvIndex, other.r)
    //{
        //std::vector<double> newContent(pow(2,rankA+rankB),0);
        //content = newContent;
        //assign(other);
    //}
    
    //explicit tensor2D (double number,uint curvIndex=0) : tensorFamily(0,0,curvIndex)
    //{
        //content.push_back(number);
    //}
    
    
    uint getContentIndex(const std::vector<uint_fast8_t>& indices) const  //Calculate index of content-vector from set of tensor indices
    {
        assert((indices.size() == rankA+rankB) &&  ("Error: requesting element with wrong number of indices") );
        
        uint contentIndex = 0;
        for(auto indexHere : indices)
        {
            assert((indexHere==0 || indexHere==1) && ("Error: given tensor index too high, expected 0 or 1") );
            contentIndex *= 2;
            contentIndex += indexHere;
        }
        return contentIndex;
    }
    
    double accessElement(const std::vector<uint_fast8_t>& indices) const override
    {
        return content.at(getContentIndex(indices));
    }
    
    void writeElement(const std::vector<uint_fast8_t>& indices, double element)
    {
        content.at(getContentIndex(indices)) = element;
    }
    
    void moveTo(const pointing& newR) override
    {
        if(arclength(r,newR)>1e-5)
        {
            r = newR;
        }
        //std::cout << "Trying to move tensor2D, this does nothing to its elements!\n";
    }

/*     /\                   */
    operator double() const 
    {
        if(rankA+rankB == 0)
        {
            return content.at(0);
        }
        else
        {
            std::cerr << "Error: trying to convert non-rank zero tensor to double. Rank is " << rankA+rankB << std::endl;
            throw std::invalid_argument( "tensor2D: Rank too high for double conversion" );
        }
    }
    
    
    tensor2D& operator*= (double multiplier)
    {
        std::transform(content.begin(), content.end(), content.begin(), [&multiplier](double& c){return c*multiplier;});
        return *this;
    }
    
    tensor2D& operator= (const tensorFamily& other) //easily evaluate every element, needs completely symmetric tensorFamily on right
    {
        assign(other);
        return *this;
    }
    
    tensor2D& operator= (const tensor2D& other) = default;
    
    tensor2D& operator= (tensor2D&& other) = default;
    

    void assign(const tensorFamily& other) //easily evaluate every element, needs completely symmetric tensorFamily on right
    {
        if( (rankA != other.rankA) || (rankB != other.rankB) || (curvIndex != other.curvIndex))
        {
            std::cerr << "Error: trying to use operator = on tensors of different type. Left tensor2D has (rankA, rankB, curvIndex) =  (" << rankA << "," << rankB << "," << curvIndex << "), Right tensorFamily has ("  << other.rankA << "," << other.rankB << "," << other.curvIndex << "), this imakes no sense" << std::endl;
            throw std::invalid_argument( "tensor2D: Different ranks in operator =" );
        }
        r = other.r;
        
        std::vector<uint_fast8_t> indices(rankA+rankB,0); //Start with only 0
        double otherElementatZeros = other.accessElement(indices);
        writeElement(indices, otherElementatZeros); //Add element at indices zero 
        
        for(int i=rankA+rankB-1;i>=0;--i) //Add one "1" to vector in every loop from back to front
        {
            std::sort (indices.begin(), indices.end());
            indices.at(i) = 1;
            double elementWithThisManyOnes = other.accessElement(indices); 
            //loop over all permutations (non-degenerate by next_permutation) and write
            do
            {
                writeElement(indices, elementWithThisManyOnes);
            } while ( std::next_permutation(indices.begin(), indices.end()) );
        }
    }
};

#endif

