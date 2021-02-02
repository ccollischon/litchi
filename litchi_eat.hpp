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
    
    
    /**  Map is generated, now write outfile (and parameterfile)  **/
    if (!params.forceOutname) //generate filename with all parameters
    {
        std::size_t fitspos = outname.find(".fits");
        if(fitspos!=std::string::npos)
        {
            outname = outname.substr(0,fitspos); //remove given .fits ending, will be added later again
        }
        char mintmaxtnumt[44];
        if(params.numt==1) //if just one threshold write that, else write mint_maxt_numt
        {
            sprintf(mintmaxtnumt,"%g", params.mint); //printf %g for nicer formatting
        }
        else
        {
            sprintf(mintmaxtnumt,"%g_%g_%d_%s", params.mint,params.maxt,params.numt, params.linThresh ? "lin" : "log"); //printf %g for nicer formatting
        }
        //TODO add smooth option
        outname = outname +"_"+ std::to_string(params.rankA) +"-"+ std::to_string(params.rankB) +"-"+ std::to_string(params.curvIndex)+ "_Nside="+std::to_string(params.Nside) + "_smooth="+"0_thresh=" + mintmaxtnumt + ".fits";
    }
    else //generate textfile with params and use given filename
    {
        std::size_t fitspos = outname.find(".fits");
        std::string paramfilename;
        if(fitspos!=std::string::npos)
        {
            paramfilename = outname.substr(0,fitspos)+"_params.txt";
        } else{
            paramfilename = outname+"_params.txt";
            outname += ".fits";
        }
        
        std::ofstream file(paramfilename);
        file << "rankA " << params.rankA << "\n";
        file << "rankB " << params.rankB << "\n";
        file << "curvIndex " << params.curvIndex << "\n";
        file << "Nside " << params.Nside << "\n";
        file << "smooth " << "TODO" << "\n"; //TODO
        file << "mint " << params.mint << "\n";
        file << "maxt " << params.maxt << "\n";
        file << "numt " << params.numt << "\n";
        file << "lin/logThresh " << (params.linThresh ? "lin" : "log") << "\n";
        file.close();
    }
    
    std::filesystem::path f{outname};
    
    if(!std::filesystem::exists(f.parent_path()))
    {
        std::cout << "Path does not exist: " << f.parent_path() << ", creating..." << std::endl;
        std::filesystem::create_directories(f.parent_path());
    }
    if (std::filesystem::exists(f))
    {
        std::cout << "File already exists, deleting old file ..." << std::endl;
        fitshandle handle;
        handle.delete_file(outname);
    }
    
    std::cout << "Writing file " << outname << std::endl;
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
