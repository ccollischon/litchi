// module load gcc/10
// g++ -fsanitize=address -g -Wall  -c regnitz.cpp -I /userdata/data/collischon/Healpix_3.70/include/  -std=c++20 -fopenmp -fconcepts
// g++ -fsanitize=address -g -Wall -o regnitz regnitz.o -L /userdata/data/collischon/Healpix_3.70/lib/ -lhealpix_cxx -std=c++20 -fopenmp -fconcepts
// cmake -D CMAKE_CXX_COMPILER=/software/Ubuntu-20.04/Programming/gcc/10.2.0/bin/g++ ../repo


#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>
#include <cmath>
#include <type_traits>
#include <chrono>
#include <typeinfo>


#include "healpix_cxx/healpix_map.h"
//#include "healpix_cxx/healpix_data_io.h"
#include "healpix_cxx/healpix_map_fitsio.h"

#include "tensor2D.hpp"
#include "tensorOperations.hpp"
#include "tensorFamily.hpp"
#include "minkTensorIntegrand.hpp"
#include "litchi_pulp.hpp"
#include "litchi_peel.hpp"
#include "litchi_eat.hpp"

const double pi = 3.14159265358979;


pointing operator*(const double& left, const pointing& right)
{
    pointing retpoint(left*right.theta, left*right.phi);
    return retpoint;
}
pointing operator*(const pointing& left, const double& right)
{
    pointing retpoint(right*left.theta, right*left.phi);
    return retpoint;
}
/**/





using namespace std;

int main(int argc,char **argv)
{
    string inname = "../litchi/COM_CMB_IQU-smica_2048_R3.00_hm1.fits", outname = "testmap.fits";
    
    struct {
        uint rankA=0, rankB=0, curvIndex=0, numt=1, Nside=0, smooth=0;
        double mint=0, maxt=1;
        bool linThresh=true, forceOutname=false, useTrace=true;
    } params;
    
    
    vector<string> arguments(argv + 1, argv + argc);    
    for(uint i=0; i<arguments.size(); i++)
    {
        string thisArg = arguments.at(i);
        if(thisArg=="--infile" || thisArg=="-i")
        {
            inname = arguments.at(++i);
        }
        else if(thisArg=="--outfile" || thisArg=="-o")
        {
            outname = arguments.at(++i);
        }
        else if(thisArg=="--rankA" || thisArg=="-A")
        {
            int number = stoi(arguments.at(++i));
            if(number>=0) params.rankA = uint(number);
            else
            {
                std::cerr << "Illegal value after rankA: " << number << " , must be positive integer or zero \n";
                return 1;
            }
        }
        else if (thisArg=="--rankB" || thisArg=="-B")
        {
            int number = stoi(arguments.at(++i));
            if(number>=0) params.rankB = uint(number);
            else
            {
                std::cerr << "Illegal value after rankB: " << number << " , must be positive integer or zero \n";
                return 1;
            }
        }
        else if (thisArg=="--curvI" || thisArg=="-C" || thisArg=="-c")
        {
            int number = stoi(arguments.at(++i));
            if(number>=0 && number<=2) params.curvIndex = uint(number);
            else
            {
                std::cerr << "Illegal value after curvI: " << number << " , must be 0, 1, or 2 \n";
                return 1;
            }
        }
        else if (thisArg=="--smooth" || thisArg=="-s")
        {
            int number = stoi(arguments.at(++i));
            if((number&(number-1))==0) params.smooth = uint(number);
            else
            {
                std::cerr << "Illegal value after smooth: " << number << " , must be power of two or zero \n";
                return 1;
            }
        }
        else if (thisArg=="--mint")
        {
            params.mint = stod(arguments.at(++i));
        }
        else if (thisArg=="--maxt")
        {
            params.maxt = stod(arguments.at(++i));
        }
        else if (thisArg=="--numt")
        {
            int number = stoi(arguments.at(++i));
            if(number>=1) params.numt = uint(number);
            else
            {
                std::cerr << "Illegal value after numt: " << number << " , must be positive integer \n";
                return 1;
            }
        }
        else if (thisArg=="--nside" || thisArg=="--Nside")
        {
            int number = stoi(arguments.at(++i));
            if( (number&(number-1))==0 ) params.Nside = uint(number);
            else
            {
                std::cerr << "Illegal value after nside: " << number << " , must be power of two \n";
                return 1;
            }
        }
        else if (thisArg=="--linThresh" || thisArg=="--linthresh")
        {
            params.linThresh = true;
        }
        else if (thisArg=="--logThresh" || thisArg=="--logthresh")
        {
            params.linThresh = false;
        }
        else if (thisArg=="--forceOutname" || thisArg=="--forceoutname")
        {
            params.forceOutname = true;
        }
        else if (thisArg=="--trace")
        {
            params.useTrace = true;
        }
        else if (thisArg=="--EVquotient" || thisArg=="--evquotient")
        {
            params.useTrace = false;
        }
        else
        {
            std::cerr << "Illegal argument: " << thisArg << "\n";
            return 1;
        }
    }
    
    
    auto start = chrono::high_resolution_clock::now();
    
    makeHealpixMinkmap(inname, params, outname); //TODO schauen dass rankB nicht bei curvIndex 0
    
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "Execution took " << duration.count()/1000. << " seconds" << endl;
    
    return 0;
}
