#ifndef litchi_wrap
#define litchi_wrap

#include "paramStruct.hpp"

#include <string>

/** \file litchi_wrap.hpp
 * \brief Wrapper functions to create and save Healpix-minkmaps for given name and parameters
 */

/*!
 * Wrapper function that calls creation of actual minkmap from given input filename, params, applying mask
 * \param inname Filename of input Healpix map
 * \param params Struct containing Minkowski map generation parameters
 * \param outname Path to and file prefix of outputfile
 */
void makeSingleMinkmap(std::string inname, paramStruct params, std::string outname);

/*!
 * Wrapper function that calls creation of sequence of minkmaps from given input filename, params, applying mask; using numt as the number of minkmaps at thresholds between mint and maxt instead of averaging one map over several thresholds
 * \param inname Filename of input Healpix map
 * \param params Struct containing Minkowski map generation parameters
 * \param outname Path to and file prefix of outputfile
 */
void makeSequence(std::string inname, paramStruct params, std::string outname);

/*!
 * Wrapper function that calls makeSequence or makeSingleMinkmap for given input filename, params
 * \param inname Filename of input Healpix map
 * \param params Struct containing Minkowski map generation parameters
 * \param outname Path to and file prefix of outputfile
 */
void makeMinkmap(std::string inname, paramStruct params, std::string outname);


#endif
