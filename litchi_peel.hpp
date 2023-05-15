#ifndef litchi_peel
#define litchi_peel

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>

#include "healpix_cxx/healpix_map.h"
#include "litchi_pulp.hpp"

/** \file litchi_peel.hpp
 * \brief Everything between minkmap and Healpix map, as well as helper functions for creating vectors with numbers at constant intervals and masking
 */

/// Create a vector of numt equally, linearly spaced doubles between mint and maxt, including mint and maxt
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

/// Create a vector of numt equally, logarithmically spaced doubles between mint and maxt, including mint and maxt
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

/** Apply mask to input image using threshold. All masked pixels are set to NAN
 * \param map Input image
 * \param mask Mask image, should be of same size and numbering scheme as map. May contain any double values; threshold will be applied. Convention: 1 = pixel not masked, 0 = pixel masked
 * \param thresh Threshold to be applied to mask. Any pixel that is below this value in the mask will be set to NAN in the input image
 * 
 */
void maskMap(Healpix_Map<double>& map, const Healpix_Map<double>& mask, double thresh = 0.9)
{
    auto nside = map.Nside();
    if(map.Nside() != mask.Nside())
    {
        std::cerr<< "Error: Nside of mask not equal to Nside of input!" << std::endl;
        throw std::invalid_argument( "maskMap: Incompatible Nside!" );
    }
    if(map.Scheme()!=mask.Scheme())
    {
        std::cerr<< "Error: Scheme of mask not equal to Scheme of input!" << std::endl;
        throw std::invalid_argument( "maskMap: Incompatible scheme!" );
    }
    
    for(int i=0; i<12*nside*nside; ++i)
    {
        if(mask[i]<thresh) map[i] = NAN;
    }
    
}

///Struct for giving minkmaps normal pixel numbering
template <minkmapFamilyType maptype>
struct normalHealpixInterface 
{
    maptype& baseminkmap;
    
    explicit normalHealpixInterface(maptype& othermap) : baseminkmap(othermap) {}
    
    /**
     * Pixel value as a linear combination of tensors at given Healpix pixel, interpolated from surrounding minkmap pixels
     * \brief Tensor value at given Healpix pixel
     * \param pixnum Pixel number
     * \return minkTensorStack with linear combination of minkmap pixels
     */
    minkTensorStack at(int pixnum) const;
    
    
    /**
     * Actual conversion of whole minkmap into Healpix map
     */
    template <minkTensor tensortype>
    Healpix_Map<double> toHealpix(double (&func)(const tensortype&), double smoothRad, int outputNside) const;
    
    ///Return 0 if pixnum not polar, 1 if north, 2 if south
    uint ispolar(int pixnum) const 
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

template <minkmapFamilyType maptype>
minkTensorStack normalHealpixInterface<maptype>::at(int pixnum) const
{
    #ifdef THISISPYTHON
        if (PyErr_CheckSignals() != 0)
            throw pybind11::error_already_set();
    #endif
    
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
            
    minkTensorStack output(baseminkmap.rankA,baseminkmap.rankB,baseminkmap.curvIndex, baseminkmap.originalMap.pix2ang(pixnum));
    for(int minkpix : westernNeighborship)
    {
        if(minkpix != -1)
        {
            output += baseminkmap.at(minkpix); //parallel transport, not just add. baseminkmap-pixels are already weighted with 1/nr of times they appear here, DONE in minkTensorStack +=
        }
    }
    return output;
}

/** Generates scalar Healpix-type map from baseminkmap via specified function
 * \param input normalHealpixInterface containing desired minkmap
 * \param func Function accepting tensor and returning scalar, e.g. trace, eigenValueQuotient, or eigenVecDir
 * \param smoothRad Window radius [rad] for input pixels included in each output pixel
 * \param outputNside Nside of output map
 * \return Healpix_Map of desired minkmap ready for saving to file
 */
template <minkmapFamilyType maptype>
template <minkTensor tensortype>
Healpix_Map<double> normalHealpixInterface<maptype>::toHealpix(double (&func)(const tensortype&), double smoothRad, int outputNside) const
{
    Healpix_Map<double> map(outputNside, baseminkmap.originalMap.Scheme(), SET_NSIDE);
    
    double smoothFactor = (baseminkmap.originalMap.Nside()/outputNside)*(baseminkmap.originalMap.Nside()/outputNside); //Amount of input pixels covered by one output pixel
    
    auto npix = map.Npix();
    int step = (outputNside <= 16) ? npix/16 : npix/64;
    
    #pragma omp parallel for
    for(int pixel=0; pixel<npix; ++pixel)
    {
        
        if(!(pixel%step))
        {
            std::cout << "Converting pixel " << pixel << "/" << npix << "...\n";
        }
        
        if(smoothRad>0)
        {
            auto pixelsNearbyRange = baseminkmap.originalMap.query_disc(map.pix2ang(pixel), smoothRad);
            std::vector<int> pixelsNearby = pixelsNearbyRange.toVector();
        
            minkTensorStack tensorHere(baseminkmap.rankA, baseminkmap.rankB, baseminkmap.curvIndex, map.pix2ang(pixel));
            for(auto pixelToAdd : pixelsNearby)
            {
                tensorHere += at(pixelToAdd); //parallel transport, not just add, DONE in minkTensorStack +=
            }
            double norm = smoothFactor/(pixelsNearby.size());//normalize such that sum over all pixels remains same
            if(func == trace<minkTensorStack>) tensorHere *= norm; //TODO check if this makes sense
            map[pixel] = func( tensorHere );
        }
        else
        {
            minkTensorStack tensorHere = at(pixel);
            map[pixel] = func( tensorHere );
        }
    }
    return map;
}

#endif
