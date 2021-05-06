#ifndef litchi_peel
#define litchi_peel

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>

#include "healpix_cxx/healpix_map.h"
#include "litchi_pulp.hpp"

std::vector<double> makeIntervals_lin(double mint, double maxt, uint numt)
{
    std::vector<double> thresholds;
    double interval = (maxt-mint)/(numt-1); //numt-1 because we want mint and maxt to be in thresholds. If numt=1 this is inf, but for loop below is ignored
    
    thresholds.push_back(mint);
    for(uint i=1; i<numt; ++i)
    {
        thresholds.push_back(mint+i*interval);
    }
    return thresholds;
}

std::vector<double> makeIntervals_log(double mint, double maxt, uint numt)
{
    std::vector<double> thresholds;
    double interval = (log(maxt)-log(mint))/(numt-1); //numt-1 because we want mint and maxt to be in thresholds. If numt=1 this is inf, but for loop below is ignored
    
    thresholds.push_back(mint);
    for(uint i=1; i<numt; ++i)
    {
        thresholds.push_back(exp(log(mint)+i*interval));
    }
    return thresholds;
}

template <typename maptype, typename std::enable_if_t<std::is_base_of<minkmapFamily,maptype>::value>* = nullptr >
struct normalHealpixInterface //class only for giving minkmaps normal pixel numbering
{
    maptype& baseminkmap;
    //normalHealpixInterface(Healpix_Map<double>& map, uint rank1 = 0, uint rank2 = 0, uint curvind = 0) : minkmapFamily(map, rank1, rank2, curvind) {}
    
    explicit normalHealpixInterface(maptype& othermap) : baseminkmap(othermap) {}
    
    tensor2D at(int pixnum) const
    {
        fix_arr<int, 8> neighbors; //neighbors of this pixel
        baseminkmap.originalMap.neighbors(pixnum,neighbors);
        std::vector<int> westernNeighborship{pixnum, neighbors[0],neighbors[1],neighbors[2]}; // non-polar: {E, SW, W, NW} corners, N-polar: {E, S, notacorner, W} corners, S-polar: {E, W, notacorner, N} in minkmap, replace notacorner
        
        
        uint pole = ispolar(pixnum);
        switch(pole) //for polar pixels replace useless neighbors[1] with pole pixnum
        {
            case 1: 
                westernNeighborship.at(2) = -11; //11ORTH POLE
                break;
            case 2: 
                westernNeighborship.at(2) = -5; //5OUTH POLE
        }
        
        if(neighbors[1]==-1) //west does not exist for some pixels, need north instead
        {
            westernNeighborship.at(2) = neighbors[3];
        }
                
        tensor2D output(baseminkmap.rankA,baseminkmap.rankB,baseminkmap.curvIndex);
        for(int minkpix : westernNeighborship)
        {
            if(minkpix != -1)
            {
                output = output + baseminkmap.at(minkpix); //TODO parallel transport, not just add. baseminkmap-pixels are already weighted with 1/nr of times they appear here
            }
        }
        return output;
    }
    
    uint ispolar(int pixnum) const //return 0 if not polar, 1 if north, 2 if south
    {
        int nside = baseminkmap.originalMap.Nside();
        int npix = baseminkmap.originalMap.Npix();
        
        if(baseminkmap.originalMap.Scheme()==RING)
        {
            if(pixnum<=3 && pixnum>=0) return 1;
            else if(pixnum<npix && pixnum>=npix-4) return 2;
        }
        else
        {
            int nsidesquared = nside*nside;
            if(pixnum==(nsidesquared-1) || pixnum==(nsidesquared*2-1) || pixnum==(nsidesquared*3-1) || pixnum==(nsidesquared*4-1)) return 1;
            else if(pixnum==(nsidesquared*8) || pixnum==(nsidesquared*9) || pixnum==(nsidesquared*10) || pixnum==(nsidesquared*11)) return 2;
        }
        return 0;
    }
};

template <typename tensortype>
Healpix_Map<double> HealpixFromMinkmap(const normalHealpixInterface<auto>& input, double func(tensortype), uint smooth) // generates scalar Healpix-type map from Minkmap(sum) via specified function
{
    int outputNside = input.baseminkmap.originalMap.Nside();
    //if smoothing: reduce resolution of outputmap
    if(smooth>1)
    {
        if((int)smooth>outputNside)
        {
            std::cerr<< "Error: smooth > Nside of input, this is not possible! smooth=" << smooth << ", Nside="<< outputNside << std::endl;
            throw std::invalid_argument( "HealpixFromMinkmap: Invalid smooth" );
        }
        outputNside /= smooth;
    }
    Healpix_Map<double> map(outputNside, input.baseminkmap.originalMap.Scheme(), SET_NSIDE);
    
    
    #pragma omp parallel for
    for(int pixel=0; pixel<map.Npix(); ++pixel)
    {
        if(!(pixel%10000))
        {
            std::cout << "Converting pixel " << pixel << "..." << std::endl;
        }
        //smooth input
        //pointing thiscenter = map.pix2ang(pixel);
        //double radius = 1.; //TODO properly based on smooth, smoothing radius evtl larger than outputpixel
        if(smooth>1)
        {
            std::vector<vec3> corners;
            map.boundaries(pixel, 1, corners); //find corners of pixel in outputmap (larger pixels)
            std::vector<pointing> cornersPoint;
            for(auto corner : corners) {cornersPoint.push_back(pointing(corner));}
        
            rangeset<int> pixelsNearbyRange = input.baseminkmap.originalMap.query_polygon(cornersPoint); //pixels in original map contained in pixel of outputmap
            std::vector<int> pixelsNearby = pixelsNearbyRange.toVector();
        
            tensor2D tensorHere(input.baseminkmap.rankA, input.baseminkmap.rankB, input.baseminkmap.curvIndex);
            for(auto pixelToAdd : pixelsNearby)
            {
                tensorHere = tensorHere+input.at(pixelToAdd); //TODO parallel transport, not just add
            }
            double norm = double(smooth*smooth)/(pixelsNearby.size());//normalize such that sum over all pixels remains same
            tensorHere *= 1./norm;
            map[pixel] = func( tensorHere );
        }
        else
        {
            tensor2D tensorHere = input.at(pixel);
            map[pixel] = func( tensorHere );
        }
    }
    return map;
}

#endif
