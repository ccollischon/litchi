// g++ regnitz.cpp -I /userdata/data/collischon/Healpix_3.70/include/ -L /userdata/data/collischon/Healpix_3.70/lib/ -lhealpix_cxx

#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>
#include <cmath>
#include <type_traits>

#include "healpix_cxx/healpix_map.h"
#include "healpix_cxx/healpix_data_io.h"
#include "healpix_cxx/healpix_map_fitsio.h"

#include "tensor2D.hpp"

double pi = 3.14159265358979;

uint binomialCoeff(uint n, uint k) 
{ 
    uint res = 1; 
  
    // Since C(n, k) = C(n, n-k) 
    if (k > n - k) 
        k = n - k; 
  
    // Calculate value of 
    // [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1] 
    for (uint i = 0; i < k; ++i) { 
        res *= (n - i); 
        res /= (i + 1); 
    } 
  
    return res;
}
/*
class tensor {
    private:
    public:
    virtual double accessElement(std::vector<uint_fast8_t> indices) = 0;
    tensor(uint rankA, rankB) : rankA(rankA), rankB(rankB) {}
}*/

// Sum expression template, operator+ of minkTensorIntegrand should return this
template <typename right0, typename left0>
struct minkTensorSum
{
    const right0 rhs;
    const left0 lhs;
    
    template <typename right1, typename left1>
    minkTensorSum(const right1& right,const left1& left)  :
        rhs(right),
        lhs(left)
    {
    }
    
    double accessElement(std::vector<uint_fast8_t> indices) const
    {
        return rhs.accessElement(indices) + lhs.accessElement(indices);
    }
    
};

//like minkTensorSum, but with * instead of + and scalar needs to come after tensor
template <typename tensortype, typename scalar>
struct minkTensorTimes
{
    const tensortype mytensor;
    const scalar myscalar;
    
    template <typename tens, typename scal>
    minkTensorTimes(const tens& te, const scal& sc)  :
        mytensor(te),
        myscalar(sc)
    {
    }
    
    double accessElement(std::vector<uint_fast8_t> indices) const
    {
        return mytensor.accessElement(indices) * myscalar;
    }
    
};

//struct to represent the r^a x n^b part of minkowski tensor calculation, only calculates value when necessary
class minkTensorIntegrand{
    public:
    const uint rankA, rankB;
    const uint curvIndex;
    pointing r, n;
    
    
    minkTensorIntegrand(uint rank1, uint rank2, uint curvInd = 0) : rankA(rank1), rankB(rank2), curvIndex(curvInd) 
    { //simple constructor for empty tensor, just give ranks
    }
    
    minkTensorIntegrand(uint rank1, uint rank2, uint curvInd, pointing rNew, pointing nNew) : rankA(rank1), rankB(rank2), curvIndex(curvInd) 
    { //constructor with all necessary data
        r = rNew;
        n = nNew;
    }
    
    double accessElement(std::vector<uint_fast8_t> indices) const // tensor element with indices i_1 ... i_{a+b} is sum over all permutations of indices of multiplication of 
    {
        uint indicesSize = indices.size();
        if(indicesSize != rankA+rankB) //Correct number of indices
        {
            std::cerr << "Error: requesting element with wrong number of indices: " << indices.size() << " expected number: " << rankA+rankB << std::endl;
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
        uint numberOfOnes = std::accumulate(indices.begin(), indices.end(), 0);
        
        std::sort (indices.begin(), indices.end());        
        //loop over all permutations (non-degenerate by next_permutation), total number over num of 1
        do //permutations of indices
        {
            //for(unsigned char index : indices) std::cout << (int)index;
            //std::cout << std::endl;
            
            double summand = 1;
            uint numberOfROnes = std::accumulate(indices.begin(), indices.begin()+rankA, 0); //Number of ones among indices of r / n
            uint numberOfNOnes = std::accumulate(indices.begin()+rankA, indices.end(), 0);
            //std::cout << numberOfROnes << " " << numberOfNOnes << std::endl;
            summand *= pow(r.theta, numberOfROnes); //different r contributions, depending on index
            summand *= pow(r.phi, rankA-numberOfROnes); // Might need Jacobian, because Polar TODO
            
            summand *= pow(n.theta, numberOfNOnes); //same for n
            summand *= pow(n.phi, rankB-numberOfNOnes);
            
            returnval += summand;
        } while ( std::next_permutation(indices.begin(), indices.end()) );
        
        
        //Calculate binomialCoeff, not factorial because next_permutation  does not create degeneracies
        
        
        returnval /= binomialCoeff(indicesSize, numberOfOnes);
        return returnval;
        
    }
    
};

template<typename left, typename right>
minkTensorSum<left,right> operator + (const left& lhs, const right& rhs)
{
    minkTensorSum<left, right> returnval (lhs, rhs);
    return returnval;
}

//Tensor times scalar, once with scalar on right and once with scalar on left
template<typename left, typename right , typename std::enable_if<std::is_arithmetic<right>::value>::type* = nullptr>
minkTensorTimes<left,right> operator* (const left& lhs, const right& rhs)
{
    return minkTensorTimes<left, right> (lhs, rhs);
}

template<typename left, typename right , typename std::enable_if<std::is_arithmetic<left>::value>::type* = nullptr>
minkTensorTimes<right,left> operator* (const left& lhs, const right& rhs)
{
    return minkTensorTimes<right, left> (rhs, lhs);
}




using namespace std;

int main(int argc,char **argv) {
    //Healpix_Map<double> map ;
    //string inname = "COM_CMB_IQU-smica_2048_R3.00_hm1.fits";
    //Healpix_Map<double> map = read_Healpix_map_from_fits<double>(inname, 1, 2);
    
    pointing r(0.35*pi, 1.5*pi);
    pointing n(0.01*pi, 0.2*pi);
    pointing nr(0.45*pi, 1.7*pi);
    minkTensorIntegrand mytensor(3,4,0,r,n);
    minkTensorIntegrand mytensor2(3,4,0,r,nr);
    cout << mytensor.accessElement({1,1,0,1,1,0,0}) << endl;
    cout << mytensor2.accessElement({1,1,0,1,1,0,0}) << endl;
    //auto timestensor = 3*mytensor2;
    int drei = 3;
    auto newtensor = drei*(3*(mytensor + mytensor));
    //auto newtensor = minkTensorTimes<minkTensorSum<minkTensorIntegrand,minkTensorIntegrand>,int>(minkTensorSum<minkTensorIntegrand,minkTensorIntegrand>(mytensor,mytensor2),3);
    
    cout << newtensor.accessElement({1,1,0,1,1,0,0}) << endl;
    cout << mytensor.accessElement({1,1,0,1,1,0,0}) << endl;
    cout << mytensor2.accessElement({1,1,0,1,1,0,0}) << endl;
    
    /*
    tensor2D testTensor(2,1,0);
    testTensor.writeElement({0,0,1},5);
    tensor2D testTensor2(2,1,0);
    testTensor2.writeElement({0,0,1},5);
    testTensor2.writeElement({1,0,0},3);
    minkTensorSum<tensor2D,tensor2D> addTensor = testTensor + testTensor2;
    //int drei = 3;
    auto addTensor2 = 3*addTensor;
    
    cout << addTensor.accessElement({0,0,1}) << " " << addTensor.accessElement({0,1,0}) << " " << addTensor.accessElement({1,0,0}) << endl;
    cout << addTensor2.accessElement({0,0,1}) << " " << addTensor2.accessElement({0,1,0}) << " " << addTensor2.accessElement({1,0,0}) << endl;
    
*/
    //write_Healpix_map_to_fits("outfile_test.fits", map, PLANCK_FLOAT32);
    
    cout << "done with test program! " << endl;

    return 0;
}
