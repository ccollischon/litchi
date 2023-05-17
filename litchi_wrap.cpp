#include "litchi_eat.hpp"
#include "paramStruct.hpp"

#include <vector>


/*!
 * Wrapper function that calls creation of actual minkmap from given input filename, params, applying mask
 * \param inname Filename of input Healpix map
 * \param params Struct containing Minkowski map generation parameters
 * \param outname Path to and file prefix of outputfile
 */
void makeSingleMinkmap(std::string inname, paramStruct params, std::string outname)
{
    params.sequence = false; //should already be false when calling function, can mess up outname otherwise
    Healpix_Map<double> map = read_Healpix_map_from_fits<double>(inname, 1, 2);
    
    prepareMap(map,params.Nside,params.maskname,params.maskThresh);
    
    makeHealpixMinkmap(map, params, outname);
}

/*!
 * Wrapper function that calls creation of sequence of minkmaps from given input filename, params, applying mask; using numt as the number of minkmaps at thresholds between mint and maxt instead of averaging one map over several thresholds
 * \param inname Filename of input Healpix map
 * \param params Struct containing Minkowski map generation parameters
 * \param outname Path to and file prefix of outputfile
 */
void makeSequence(std::string inname, paramStruct params, std::string outname)
{
    Healpix_Map<double> map = read_Healpix_map_from_fits<double>(inname, 1, 2);
    params.sequence = true; //should already be true when calling this function, but just to be sure. Can mess up outname otherwise
    
    prepareMap(map,params.Nside,params.maskname,params.maskThresh);
    
    const std::vector<double> thresholds = params.linThresh ? makeIntervals_lin(params.mint, params.maxt, params.numt) : makeIntervals_log(params.mint, params.maxt, params.numt);
    for(uint i=0; i<params.numt; ++i)
    {
        paramStruct paramsHere(params);
        paramsHere.numt = 1;
        paramsHere.mint = thresholds.at(i);
        makeHealpixMinkmap(map, paramsHere, outname, i);
    }
}

/*!
 * Wrapper function that calls makeSequence or makeSingleMinkmap for given input filename, params
 * \param inname Filename of input Healpix map
 * \param params Struct containing Minkowski map generation parameters
 * \param outname Path to and file prefix of outputfile
 */
void makeMinkmap(std::string inname, paramStruct params, std::string outname)
{
    if(params.sequence) makeSequence(inname, params, outname);
    else makeSingleMinkmap(inname, params, outname);
}


