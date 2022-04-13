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
 * \brief Higher-level functions to create and save Healpix-Minkmaps generated from input map
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

///Helper struct containing all necessary Minkmap generation parameters
struct paramStruct{
        uint rankA{0}, rankB{0}, curvIndex{0}, numt{1}, Nside{0}, smooth{0};
        double mint{0}, maxt{1}, maskThresh{0.9};
        bool linThresh{true}, forceOutname{false}, sequence{false};
        std::string function{"trace"}, maskname{""};
};

void checkParams(const auto &obj) {
    
    if(obj.curvIndex==0 && obj.rankB)
    {
        std::cerr << "Error: rankB > 0 not possible when curvIndex=0. rankB = " << obj.rankB << std::endl;
        throw std::invalid_argument("rankB && curvIndex!=0");
    }
    if(obj.curvIndex==0 && obj.rankA)
    {
        std::cerr << "Error: rankA > 0 not implemented when curvIndex=0. rankA = " << obj.rankA << std::endl;
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
/**
 * Formats given outname to contain all relevant parameters if params.forceOutname is false
 * \param outname Path and desired file prefix with or without .fits ending. Parameters are added accordingly
 * \param params Struct containing Minkowski map generation parameters to write into header: Nside, rankA, rankB, curvIndex, mint, maxt, numt, smooth, linThresh, function, forceOutname, sequence
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
        char mintmaxtnumt[44];
        if(params.numt==1) //if just one threshold write that, else write mint_maxt_numt
        {
            sprintf(mintmaxtnumt,"%.3e", params.mint); //printf %g for nicer formatting
        }
        else
        {
            sprintf(mintmaxtnumt,"%.3e_%.3e_%d_%s", params.mint,params.maxt,params.numt, params.linThresh ? "lin" : "log"); //printf %g for nicer formatting
        }
        std::string funString = (params.function=="trace") ? "_tr" : (params.function=="EVQuo") ? "_evq" : (params.function=="EVDir") ? "_evd" : "_error";
        char maskString[20];
        if(params.maskname!="") sprintf(maskString,"_mask_%4.2f", params.maskThresh);
        outname = outname +"_"+ std::to_string(params.rankA) +"-"+ std::to_string(params.rankB) +"-"+ std::to_string(params.curvIndex) + funString + "_Nside="+std::to_string(params.Nside) + "_smooth="+std::to_string(params.smooth) + "_thresh="+mintmaxtnumt + maskString + ".fits";
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
 * Write given map to file specified by outname caontaining params in header
 * Function checks if given path for outputfile exists and creates it if not. Files with same name are overwritten.
 * \param outputmap Healpix map to be saved to file
 * \param params Struct containing Minkowski map generation parameters to write into header: Nside, rankA, rankB, curvIndex, mint, maxt, numt, smooth, function, useTrace, forceOutname, sequence
 * \param outname Name of output file
 */
void writeToFile(const Healpix_Map<double>& outputmap, const paramStruct& params, std::string outname)
{
    fitshandle handle;
    std::filesystem::path f{outname};
    
    if(!std::filesystem::exists(f.parent_path()))
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
    handle.add_comment("Minkowski map created using Litchi with the following parameters: ");
    handle.set_key("rankA", (int)params.rankA, "First rank of Minkowski tensor (r)");
    handle.set_key("rankB", (int)params.rankB, "Second rank of Minkowski tensor (n)");
    handle.set_key("curvIndex ", (int)params.curvIndex, "Curvature (bottom) index of Minkowski tensor");
    handle.set_key("Nside ", (int)params.Nside, "Nside of input map");
    handle.set_key("smooth ", (int)params.smooth, "Factor by which output was smoothed");
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
        handle.set_key("maskThresh", params.maskThresh, "All mask pixels above this threshold were treated as not masked");
    }
    
    std::vector<std::string> funStrings; //Description of function used to generate scalar from tensor
    if(params.function=="trace"){
        funStrings = {"trace", "Trace of Tensor"};
    }
    else if (params.function=="EVQuo"){
        funStrings = {"EV quotient", "Eigenvalue quotient of tensor"};
    }
    else if (params.function=="EVDir"){
        funStrings = {"EV direction", "Direction of eigenvector with largest eigenvalue"};
    }
    else funStrings = {params.function,"unknown function"};
    
    handle.set_key("Function ", funStrings.at(0), funStrings.at(1));
    
    std::cout << "Writing file " << outname << std::endl;
    write_Healpix_map_to_fits(handle, outputmap, PLANCK_FLOAT32);
}

/*!
 * "Wrapper"  function that creates actual minkmap from given inputmap, params. Warning: degrades input map if Nside parameter set
 * \param map Input Healpix map
 * \param params Struct containing Minkowski map generation parameters: Nside, rankA, rankB, curvIndex, mint, maxt, numt, smooth, linThresh, function, forceOutname, sequence
 * \param outname Path to and file prefix of outputfile
 */
void makeHealpixMinkmap(Healpix_Map<double>& map, paramStruct params, std::string outname, const int counter=0)
{
    checkParams(params);
    if(!(int)params.Nside) params.Nside = (uint)map.Nside();
    
    if((int)params.Nside != map.Nside())
    {
        Healpix_Map<double> degradedMap(params.Nside, map.Scheme(), SET_NSIDE);
        degradedMap.Import_degrade(map);
        map = degradedMap;
    }
    
    if(params.maskname!="")
    {
        Healpix_Map<double> mask = read_Healpix_map_from_fits<double>(params.maskname, 1, 2);
        maskMap(map, mask, params.maskThresh);
    }
    
    const std::vector<double> thresholds = params.linThresh ? makeIntervals_lin(params.mint, params.maxt, params.numt) : makeIntervals_log(params.mint, params.maxt, params.numt);
    std::vector<minkmapSphere> maps;
    for(double thresh : thresholds)
    {
        maps.push_back(minkmapSphere(map, params.rankA, params.rankB, params.curvIndex, thresh));
    }
    const minkmapStack sumOfMaps(maps);
    
    //Probably should smooth before adding, so can parallel transport n only
    
    const auto minkmapAverage = sumOfMaps*(1./params.numt);
    const normalHealpixInterface interface(minkmapAverage);
    Healpix_Map<double> outputmap;
    if(params.function=="trace")
    {
        outputmap = interface.toHealpix(trace<minkTensorStack>,params.smooth);
    }
    else if(params.function=="EVQuo")
    {
        outputmap = interface.toHealpix(eigenValueQuotient<minkTensorStack>,params.smooth);
    }
    else if(params.function=="EVDir")
    {
        outputmap = interface.toHealpix(eigenVecDir<minkTensorStack>,params.smooth);
    }
    
    /*  Map is generated, now create outname  */
    formatOutname(outname, params, counter);
    
    /* create Fitsfile with params in Header */
    writeToFile(outputmap,params,outname);
    
}


/*!
 * "Wrapper"  function that creates actual minkmap from given input filename, params
 * \param inname Filename of input Healpix map
 * \param params Struct containing Minkowski map generation parameters: Nside, rankA, rankB, curvIndex, mint, maxt, numt, smooth, linThresh, useTrace, forceOutname, sequence
 * \param outname Path to and file prefix of outputfile
 */
void makeSingleMinkmap(std::string inname, paramStruct params, std::string outname)
{
    params.sequence = false; //should already be false when calling function, can mess up outname otherwise
    Healpix_Map<double> map = read_Healpix_map_from_fits<double>(inname, 1, 2);
    makeHealpixMinkmap(map, params, outname);
}

/*!
 * "Wrapper"  function that creates sequence of minkmaps from given input filename, params using numt as the number of minkmaps at thresholds between mint and maxt instead of averaging one map over several thresholds
 * \param inname Filename of input Healpix map
 * \param params Struct containing Minkowski map generation parameters: Nside, rankA, rankB, curvIndex, mint, maxt, numt, smooth, linThresh, useTrace, forceOutname, sequence
 * \param outname Path to and file prefix of outputfile
 */
void makeSequence(std::string inname, paramStruct params, std::string outname)
{
    Healpix_Map<double> map = read_Healpix_map_from_fits<double>(inname, 1, 2);
    params.sequence = true; //should already be true when calling this function, but just to be sure. Can mess up outname otherwise
    
    const std::vector<double> thresholds = params.linThresh ? makeIntervals_lin(params.mint, params.maxt, params.numt) : makeIntervals_log(params.mint, params.maxt, params.numt);
    for(uint i=0; i<params.numt; ++i)
    {
        paramStruct paramsHere(params);
        paramsHere.numt = 1;
        paramsHere.mint = thresholds.at(i);
        makeHealpixMinkmap(map, paramsHere, outname, i);
    }
}


#endif
