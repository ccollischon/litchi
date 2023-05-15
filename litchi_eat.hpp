#ifndef litchi_eat
#define litchi_eat


#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <cstdlib> // for sprintf

#include "healpix_cxx/healpix_map_fitsio.h"
#include "healpix_cxx/healpix_map.h"

#include "litchi_peel.hpp"

/** \file litchi_eat.hpp
 * \brief Higher-level functions to create and save Healpix-minkmaps generated from input map
 */

//C++20 solution for checking if params-struct has required members
/*
bool has_x(const auto &obj) {
    if constexpr (requires {obj.x;}) {
      return true;
    } else
      return false;
}
*/

///Helper struct containing all necessary minkmap generation parameters
struct paramStruct{
        uint rankA{0}, rankB{0}, curvIndex{0}, numt{1}, Nside{0}, smooth{0}, NsideOut{0};
        double mint{0}, maxt{1}, maskThresh{0.9}, smoothRad{0};
        bool linThresh{true}, forceOutname{false}, sequence{false};
        std::string function{"trace"}, maskname{""};
};

///Sanity check for all parameters. Throws std::invalid_argument if something goes wrong
void checkParams(const auto &obj) {
    
    if(obj.curvIndex==0 && obj.rankB)
    {
        std::cerr << "Error: rankB > 0 not possible when curvIndex = 0. rankB = " << obj.rankB << std::endl;
        throw std::invalid_argument("rankB && curvIndex!=0");
    }
    if(obj.curvIndex==0 && obj.rankA)
    {
        std::cerr << "Error: rankA > 0 not implemented when curvIndex = 0. rankA = " << obj.rankA << std::endl;
        throw std::invalid_argument("rankA && curvIndex!=0");
    }
    if(obj.rankA)
    {
        std::cerr << "Warning: rankA > 0 not implemented properly" << std::endl;
    }
    if(obj.curvIndex > 2)
    {
        std::cerr << "Illegal value for curvIndex: " << obj.curvIndex << " , must be 0, 1, or 2 \n";
        throw std::invalid_argument("curvIndex > 2");
    }
    if((obj.smooth&(obj.smooth-1)) != 0)
    {
        std::cerr << "Illegal value for smooth: " << obj.smooth << " , must be power of two or zero \n";
        throw std::invalid_argument("smooth not power of two or zero");
    }
    if((obj.Nside&(obj.Nside-1)) != 0)
    {
        std::cerr << "Illegal value for Nside: " << obj.Nside << " , must be power of two or zero \n";
        throw std::invalid_argument("Nside not power of two or zero");
    }
    if((obj.NsideOut&(obj.NsideOut-1)) != 0)
    {
        std::cerr << "Illegal value for NsideOut: " << obj.NsideOut << " , must be power of two or zero \n";
        throw std::invalid_argument("NsideOut not power of two or zero");
    }
    if(obj.NsideOut > obj.Nside)
    {
        std::cerr << "NsideOut larger than Nside; NsideOut = " << obj.NsideOut << ", Nside = " << obj.Nside << " , this makes no sense! \n";
        throw std::invalid_argument("NsideOut larger than Nside");
    }
    if(obj.smoothRad < 0)
    {
        std::cerr << "Illegal value for smoothRad: " << obj.smoothRad << " , must be non-negative \n";
        throw std::invalid_argument("smoothRad negative");
    }
    if(obj.smoothRad == 0. && (obj.NsideOut<obj.Nside) )
    {
        std::cerr << "If smoothRad==0 (no smoothing), but NsideOut < Nside, you're gonna have a bad time. smoothRad = " << obj.smoothRad << " , NsideOut = " << obj.NsideOut << ", Nside = " << obj.Nside << "\n";
        throw std::invalid_argument("no smoothing but smoothing");
    }
    if(obj.mint==0. && !obj.linThresh)
    {
        std::cerr << "Minimal threshold zero not possible with logThresh!\n";
        throw std::invalid_argument("mint==0 && logThresh");
    }
    if(obj.function!="trace" && obj.function!="EVQuo" && obj.function!="EVDir")
    {
        std::cerr << "Invalid tensor-to-scalar function given, only permits trace, EVQuo, EVDir!, but have "+obj.function+"\n";
        throw std::invalid_argument("function invalid");
    }
}

///If params.smooth is set, overwrite params.NsideOut and params.smoothRad to corresponding values. Through this, user can use either smooth or smoothRad
void backwardsCompatibilitySmooth(paramStruct& params)
{
    if(params.smooth > 0) //User has used deprecated feature
    {
        std::cerr << "Warning: using params.smooth instead of setting smoothRad and NsideOut. If you have also set NsideOut and smoothRad, they are overwritten. \n";
        
        //Set NsideOut
        if(params.smooth>params.Nside)
        {
            std::cerr<< "Error: smooth > Nside of input, this is not possible! smooth=" << params.smooth << ", Nside="<< params.Nside << std::endl;
            throw std::invalid_argument( "smooth larger than Nside" );
        }

        params.NsideOut = params.Nside/params.smooth;
        
        //Set smoothRad
        Healpix_Map<int> tempmap(params.NsideOut, RING, SET_NSIDE);
        params.smoothRad = 1.5*tempmap.max_pixrad(); //Take distance pixel center-corners in new map. In smoothed map, consider all input map pixels up to 1.5* that distance
    }
}


/**
 * Formats given outname to contain all relevant parameters if params.forceOutname is false; adds a counter if params.forceOutname is true and params.sequence is true
 * \param outname Path and desired file prefix with or without .fits ending. Parameters are added accordingly
 * \param params Struct containing Minkowski map generation parameters
 * \param counter Number of image in sequence if sequence of images is generated
 */
void formatOutname(std::string& outname, const paramStruct& params, const int counter=0)
{
    if (!params.forceOutname) //generate filename with all parameters
    {
        std::size_t fitspos = outname.find(".fit");
        if(fitspos!=std::string::npos)
        {
            outname = outname.substr(0,fitspos); //remove given .fit(s) ending, will be added later again
        }
        char mintmaxtnumt[59];
        if(params.numt==1) //if just one threshold write that, else write mint_maxt_numt. Always add smoothRad
        {
            sprintf(mintmaxtnumt,"%.3e_rad=%.3e", params.mint,params.smoothRad); //printf %g for nicer formatting
        }
        else
        {
            sprintf(mintmaxtnumt,"%.3e_%.3e_%d_%s_rad=%.3e", params.mint,params.maxt,params.numt, params.linThresh ? "lin" : "log", params.smoothRad); //printf %g for nicer formatting
        }
        std::string funString = (params.function=="trace") ? "_tr" : (params.function=="EVQuo") ? "_evq" : (params.function=="EVDir") ? "_evd" : "_error";
        char maskString[20];
        (params.maskname!="") ? sprintf(maskString,"_mask_%4.2f", params.maskThresh) : sprintf(maskString,"_nomask");
        outname = outname +"_"+ std::to_string(params.rankA) +"-"+ std::to_string(params.rankB) +"-"+ std::to_string(params.curvIndex) + funString + "_Nside="+std::to_string(params.Nside) + "_NsideOut="+std::to_string(params.NsideOut) + "_thresh="+mintmaxtnumt + maskString + ".fits";
    }
    else if(params.sequence) //for sequence cannot leave the outname unchanged, or else will overwrite one file over and over
    {
        std::cout << "Warning: creating sequence with --forceOutname active, adding counter not to overwrite everything\n";
        char buf[9];
        sprintf(buf,"%03d", counter);
        
        std::size_t fitspos = outname.find(".fit");
        if(fitspos!=std::string::npos)
        {
            outname = outname.substr(0,fitspos); //remove given .fit(s) ending, will be added later again
            outname = outname + buf + ".fits";
        }
        else
        {
            outname = outname + buf;
        }
    }
}

/**
 * Write given map to file specified by outname containing params in header.\n
 * Checks if given path for outputfile exists and creates it if not. Files with same name are overwritten. Uses PLANCK_FLOAT32 as output data type
 * \param outputmap Healpix map to be saved to file
 * \param params Struct containing Minkowski map generation parameters to write into header
 * \param outname Name of output file
 */
void writeToFile(const Healpix_Map<double>& outputmap, const paramStruct& params, std::string outname)
{
    fitshandle handle;
    std::filesystem::path f{outname};

    
    if(!std::filesystem::exists(f.parent_path()) && f.parent_path()!="")
    {
        std::cout << "Path does not exist: " << f.parent_path() << ", creating..." << std::endl;
        std::filesystem::create_directories(f.parent_path());
    }
    if (std::filesystem::exists(f))
    {
        std::cout << "File already exists, deleting old file ..." << std::endl;
        handle.delete_file(outname);
    }
    
    handle.create(outname);
    handle.add_comment("Minkowski map created using litchi with the following parameters: ");
    handle.set_key("rankA", (int)params.rankA, "First rank of Minkowski tensor (r)");
    handle.set_key("rankB", (int)params.rankB, "Second rank of Minkowski tensor (n)");
    handle.set_key("curvIndex ", (int)params.curvIndex, "Curvature (bottom) index of Minkowski tensor");
    handle.set_key("Nside ", (int)params.Nside, "Nside of input map");
    handle.set_key("NsideOut ", (int)params.NsideOut, "Nside of output map");
    handle.set_key("smooth ", (int)params.smooth, "Factor by which output was smoothed (optional)");
    handle.set_key("smoothRad ", params.smoothRad, "Smoothing window radius in rad");
    handle.set_key("mint ", params.mint, "Min threshold");
    handle.set_key("maxt ", params.maxt, "Max threshold");
    handle.set_key("numt ", (int)params.numt, "Number of thresholds (mint used if 1)");
    if(params.linThresh)
    {
        handle.set_key("lin/logThresh ", std::string("lin"), "linear or logarithmic spacing");
    }
    else
    {
        handle.set_key("lin/logThresh ", std::string("log"), "linear or logarithmic spacing");
    }
    
    if(params.maskname!="")
    {
        handle.set_key("mask", params.maskname, "mask used before generating map");
        handle.set_key("maskThresh", params.maskThresh, "mask pixels above this treated as not masked");
    }
    
    std::vector<std::string> funStrings; //Description of function used to generate scalar from tensor
    if(params.function=="trace"){
        funStrings = {"trace", "Trace of Tensor"};
    }
    else if (params.function=="EVQuo"){
        funStrings = {"EV quotient", "Eigenvalue quotient of tensor"};
    }
    else if (params.function=="EVDir"){
        funStrings = {"EV direction", "Direction of eigenvec with largest eigenval"};
    }
    else funStrings = {params.function,"unknown function"};
    
    handle.set_key("Function ", funStrings.at(0), funStrings.at(1));
    
    std::cout << "Writing file " << outname << std::endl;
    write_Healpix_map_to_fits(handle, outputmap, PLANCK_FLOAT32);
}

/*!
 * Function that creates actual minkmap from given inputmap, params
 * \param map Input Healpix map
 * \param params Struct containing Minkowski map generation parameters
 * \param outname Path to and file prefix of outputfile
 * \param counter Number of image in sequence if sequence of images is generated. Used for output name generation
 */
void makeHealpixMinkmap(Healpix_Map<double>& map, paramStruct params, std::string outname, const int counter=0)
{
    if(!(int)params.Nside)       params.Nside = (uint)map.Nside();
    if(!(int)params.NsideOut) params.NsideOut = params.Nside;
    backwardsCompatibilitySmooth(params);
    checkParams(params);
    
    if((int)params.Nside != map.Nside())
    {
        Healpix_Map<double> degradedMap(params.Nside, map.Scheme(), SET_NSIDE);
        degradedMap.Import_degrade(map);
        map = degradedMap;
    }
    
    
    const std::vector<double> thresholds = params.linThresh ? makeIntervals_lin(params.mint, params.maxt, params.numt) : makeIntervals_log(params.mint, params.maxt, params.numt);
    std::vector<minkmapSphere> maps;
    for(double thresh : thresholds)
    {
        maps.push_back(minkmapSphere(map, params.rankA, params.rankB, params.curvIndex, thresh));
    }
    const minkmapStack sumOfMaps(maps);
    
    
    const auto minkmapAverage = sumOfMaps*(1./params.numt);
    const normalHealpixInterface interface(minkmapAverage);
    Healpix_Map<double> outputmap;
    if(params.function=="trace")
    {
        outputmap = interface.toHealpix(trace<minkTensorStack>,params.smoothRad, (int)params.NsideOut);
    }
    else if(params.function=="EVQuo")
    {
        outputmap = interface.toHealpix(eigenValueQuotient<minkTensorStack>,params.smoothRad, (int)params.NsideOut);
    }
    else if(params.function=="EVDir")
    {
        outputmap = interface.toHealpix(eigenVecDir<minkTensorStack>,params.smoothRad, (int)params.NsideOut);
    }
    else
    {
        std::cerr << "Invalid tensor-to-scalar function given, only permits trace, EVQuo, EVDir!, but have "+params.function+"\n This should not happen as it is already checked in checkParams\n";
        throw std::invalid_argument("function invalid");
    }
    
    /*  Map is generated, now create outname  */
    formatOutname(outname, params, counter);
    
    /* create Fitsfile with params in Header */
    writeToFile(outputmap,params,outname);
    
}


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
    
    if(params.maskname!="")
    {
        Healpix_Map<double> mask = read_Healpix_map_from_fits<double>(params.maskname, 1, 2);
        maskMap(map, mask, params.maskThresh);
    }
    
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
    
    if(params.maskname!="")
    {
        Healpix_Map<double> mask = read_Healpix_map_from_fits<double>(params.maskname, 1, 2);
        maskMap(map, mask, params.maskThresh);
    }
    
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


#endif
