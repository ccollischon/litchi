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

/*
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
*/





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
    
    /*
    Healpix_Map<double> map = read_Healpix_map_from_fits<double>(inname, 1, 2);
    std::vector<double> thresholds = params.linThresh ? makeIntervals_lin(params.mint, params.maxt, params.numt) : makeIntervals_log(params.mint, params.maxt, params.numt);
    std::vector<minkmapSphere> maps;
    for(double thresh : thresholds)
    {
        maps.push_back(minkmapSphere(map, 0, 0, 0, thresh));
    }
    minkmapStack sumOfMaps(maps);
    auto minkmapAverage = sumOfMaps*(1./params.numt);
    auto pix1 = maps.at(0).at(1932);
    auto pix2 = maps.at(1).at(1932);
    auto pix3 = sumOfMaps.at(1932);
    auto pix4 = minkmapAverage.at(1932);
    
    cout << (float)trace(pix1) << " " << trace(pix2) << " " << (float)((1./params.numt)*trace(pix3)) << " " << (float)trace(pix4) << endl;
    
    
    
    
    pointing A(pi/2,0);
    pointing B(pi/4,0);
    pointing C(pi/2,pi/4);
    
    cout << midpoint(C,B)*(1/pi) << endl;
    
    
    flipPointing(A);
    pointing n = getN_cartesian(B,A);
    cout << n << endl;
    minkTensorIntegrand asdf(0,2,0,A,B);
    cout << eigenValueQuotient(asdf);
    
    
    
    cout << sphereArea(A, B, C)/pi << " " << sphereArea(A,C,B)/pi << " " << sphereArea(B,C,A)/pi << " " << sphereArea(B,A,C)/pi << endl;
    
    vector<double> thresholds = makeIntervals_lin(mint,maxt,numt);
    for (double threshs : thresholds) cout << threshs << " ";
    cout << endl;
    
    *//*
    //string inname = "COM_CMB_IQU-smica_2048_R3.00_hm1.fits";
    Healpix_Map<double> map = read_Healpix_map_from_fits<double>(inname, 1, 2);
    map.SetNside(256,RING);
    
    fix_arr<int, 8> allneighbors;
    map.neighbors(144,allneighbors);
    for(uint i=0;i<allneighbors.size();i++)
    {
        cout << allneighbors[i] << endl;
    }
    
    minkmapSphere testmink(map);
    normalHealpixInterface interface(testmink);
    fix_arr<int, 4> neighbors; 
    fix_arr<double, 4> weight; 
    pointing north(0,0);
    map.get_interpol(north,neighbors,weight);
    cout << "north: " << endl;
    for(uint i=0;i<neighbors.size();i++)
    {
        cout << neighbors[i] << " " << interface.ispolar(neighbors[i]) << endl;
    }
    pointing south(pi,0);
    map.get_interpol(south,neighbors,weight);
    cout << "south: " << endl;
    for(uint i=0;i<neighbors.size();i++)
    {
        cout << neighbors[i] << " " << interface.ispolar(neighbors[i]) << endl;
    }
    
    pointing random(4.76,1.345);
    random.normalize();
    map.get_interpol(random,neighbors,weight);
    cout << "random: " << endl;
    for(uint i=0;i<neighbors.size();i++)
    {
        cout << neighbors[i] << " " << interface.ispolar(neighbors[i]) << endl;
    }*/
    /*
    vec3 axis(1,0,0);
    vec3 input(0,1,0);
    vec3 testvec = rotateAroundAxis(axis, input, pi/2);
    cout << testvec << endl;
    
    
     
    minkmapSphere minkmapTest(map,0,0,0, 1e-5);
    minkmapSphere minkmapTest2(map,0,0,0, 2e-5);
    vector<int> pixels{0,1,-1,3};
    auto sum = minkmapTest+minkmapTest2 + minkmapTest;
    tensor2D tensor = sum.at(0);
    
    
    pointing interpolation  = interpPointing(r, 1.5, n, 0.5, 1);
    cout << interpolation << endl;
    interpolation  = interpPointing(n, 0.5, r, 1.5, 1);
    cout << interpolation << endl;
    */
    /*
    pointing r(0.5, 0);
    pointing n(0.5, pi);
    pointing nr(0.45*pi, 0.5*pi);
    minkTensorIntegrand mytensor(0,2,0,r,n);
    minkTensorIntegrand mytensor2(0,2,0,nr,nr);
    tensor2D testtensor = mytensor;
    cout << mytensor.accessElement({1,1}) << endl;
    cout << mytensor2.accessElement({1,1}) << endl;
    cout << testtensor.accessElement({1,1}) << endl;
    
    tensor2D zahlenfriedhof (mytensor*3+mytensor2);
    cout << zahlenfriedhof.accessElement({1,1}) << endl;
    zahlenfriedhof = mytensor + mytensor2;
    cout << zahlenfriedhof.accessElement({1,1}) << endl;
    
    auto timestensor = 3*mytensor2;
    cout << "\n" << typeid(timestensor).name() << "\n";
    
    int drei = 3;
    auto newtensor = drei*(3*(mytensor + mytensor2));
    //auto newtensor = minkTensorTimes<minkTensorSum<minkTensorIntegrand,minkTensorIntegrand>,int>(minkTensorSum<minkTensorIntegrand,minkTensorIntegrand>(mytensor,mytensor2),3);
    
    cout << newtensor.accessElement({1,1,0,1,1,0,0}) << endl;
    cout << mytensor.accessElement({1,1,0,1,1,0,0}) << endl;
    cout << mytensor2.accessElement({1,1,0,1,1,0,0}) << endl;
    
    
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
    

    return 0;
}
