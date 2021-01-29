#ifndef litchi_peel
#define litchi_peel

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <filesystem>

#include "healpix_cxx/healpix_map.h"
#include "healpix_cxx/healpix_map_fitsio.h"
#include "litchi_pulp.hpp"

std::vector<double> makeIntervals_lin(double mint, double maxt, uint numt)
{
    std::vector<double> thresholds;
    double interval = (maxt-mint)/(numt-1); //numt-1 because we want mint and maxt to be in thresholds. If numt=1 this is inf, but for loop below is ignored
    
    thresholds.push_back(mint);
    for(uint i=1; i<numt; i++)
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
    for(uint i=1; i<numt; i++)
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
    
    normalHealpixInterface(maptype& othermap) : baseminkmap(othermap) {}
    
     //TODO: adapt to 3-pixel regions, poles. idea: give poles own (negative) pixnum, check for poles and replace pix in westernNeighborship with polepixel
    tensor2D at(int pixnum) const
    {
        fix_arr<int, 8> neighbors; //neighbors of this pixel
        baseminkmap.originalMap.neighbors(pixnum,neighbors); // pixels with these numbers in minkmap are positioned around original pixel with this number
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
            std::cout << "here" << std::endl;
        }
                
        tensor2D output(baseminkmap.rankA,baseminkmap.rankB,baseminkmap.curvIndex);
        for(int minkpix : westernNeighborship)
        {
            if(minkpix != -1)
            {
                output.assign(output + baseminkmap.at(minkpix)); //TODO parallel transport, not just add. baseminkmap-pixels are already weighted with 1/nr of times they appear here
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
            if(pixnum<=4 && pixnum>=0) return 1;
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
Healpix_Map<double> HealpixFromMinkmap(const normalHealpixInterface<auto>& input, double func(tensortype)) // generates scalar Healpix-type map from Minkmap(sum) via specified function
{
    Healpix_Map<double> map(input.baseminkmap.originalMap.Nside(), input.baseminkmap.originalMap.Scheme(), SET_NSIDE);
    
    //if smoothing: reduce resolution of outputmap
    
    //smooth input
    
    //TODO parallelize here
    for(int pixel=0; pixel<map.Npix(); pixel++)
    {
        if(!(pixel%10000)) 
        {
            std::cout << "Converting pixel " << pixel << "..." << std::endl;
        }
        tensor2D tensorHere = input.at(pixel);
        map[pixel] = func( tensorHere ); //actually need to average over .at(_) of western neighborhood, because .at() refers to vertex east of pixel, so map should be interfaces
    }
    return map;
}

//C++20 solution for checking if params-struct has required members
/*
bool has_x(const auto &obj) {
    if constexpr (requires {obj.x;}) {
      return true;
    } else
      return false;
}
*/

template <typename tensortype, typename paramtype>
void makeHealpixMinkmap(Healpix_Map<double>& map, paramtype params, double func(tensortype), std::string outname)
{
    std::cout << "TODO: check params for members with C++20" << std::endl;
    if(params.Nside)
    {
        Healpix_Map<double> degradedMap(params.Nside, map.Scheme(), SET_NSIDE);
        degradedMap.Import_degrade(map);
        map = degradedMap;
    }
    std::vector<double> thresholds = params.linThresh ? makeIntervals_lin(params.mint, params.maxt, params.numt) : makeIntervals_log(params.mint, params.maxt, params.numt);
    std::vector<minkmapSphere> maps;
    for(double thresh : thresholds)
    {
        maps.push_back(minkmapSphere(map, params.rankA, params.rankB, params.curvIndex, thresh));
    }
    minkmapStack sumOfMaps(maps);
    
    //Probably should smooth before adding, so can parallel transport n only
    
    auto minkmapAverage = sumOfMaps*(1./params.numt);
    normalHealpixInterface interface(minkmapAverage);
    Healpix_Map<double> outputmap = HealpixFromMinkmap(interface,func);
    
    
    std::filesystem::path f{outname};
    if (std::filesystem::exists(f))
    {
        std::cout << "File already exists, deleting old file ..." << std::endl;
        fitshandle handle;
        handle.delete_file(outname);
    }
    
    write_Healpix_map_to_fits(outname, outputmap, PLANCK_FLOAT32);
}

template <typename tensortype, typename paramtype>
//void makeHealpixMinkmap(std::string inname, uint rankA, uint rankB, uint curvIndex, uint numt, double mint, double maxt, bool linthresh, double func(tensortype), std::string outname)
void makeHealpixMinkmap(std::string inname, paramtype params, double func(tensortype), std::string outname)
{
    Healpix_Map<double> map = read_Healpix_map_from_fits<double>(inname, 1, 2);
    makeHealpixMinkmap(map, params, func, outname);
}

#endif
