#ifndef litchi_pulp
#define litchi_pulp

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>

#include "healpix_cxx/healpix_map.h"
#include "litchi_kernel.hpp"

#include "tensorOperations.hpp"
#include "minkTensorIntegrand.hpp"
#include "geometryHelpers.hpp"

/** \file litchi_pulp.hpp
 * \brief minkmapSphere class, which does the heavy lifting when creating a minkmap
 */

extern const double pi;


///Contains "raw" (marching-square-level) minkowski tensors, returns on demand and doesn't save them because of space. 
struct minkmapSphere :  minkmapFamily{ 
    const double thresh{};
    
    
    explicit minkmapSphere(Healpix_Map<double>& map) : minkmapFamily(map) {}
    
    minkmapSphere(Healpix_Map<double>& map, uint rank1, uint rank2, uint curvind, double threshold) : minkmapFamily(map, rank1, rank2, curvind), thresh(threshold) {
        if (curvind > 2) {
            std::cerr << "Error: invalid curvIndex: " << curvind << " , this makes no sense. In 2D only up to two!" << std::endl;
            throw std::invalid_argument( "minkmapSphere: Weird curvIndex" );
        }
    }
    
    pointing pix2ang(int pixnum) const
    {
        if(pixnum>=originalMap.Npix())
        {
            std::cerr << "Error: requesting pixnum too large for Healpix_Map, pixnum is " << pixnum << ", Npix is " << originalMap.Npix() << std::endl;
            throw std::invalid_argument("minkmapSphere::pix2ang: pixnum too high");
        }
        else if(pixnum<0 && pixnum!=-5 && pixnum!= -11)
        {
            std::cerr << "Error: requesting negative undefined pixnum, pixnum is " << pixnum << ", Npix is " << originalMap.Npix() << std::endl;
            throw std::invalid_argument("minkmapSphere::pix2ang: pixnum invalid");
        }
        
        if(pixnum==-5) return pointing(0,pi);
        if(pixnum==-11) return pointing(0,0);
        
        std::vector<vec3> corners;
        originalMap.boundaries(pixnum, 1, corners); //find corners of original pixel in original map
        return pointing(corners.at(3)); //position east of pixel is center of vertex
    }
    
    ///returns Healpix-pixel numbers of pixels surrounding given minkmap pixel, -11 = north pole, -5 = south pole
    std::vector<int> minkmapPixelNeighbors(int pixnum) const
    {
        if(pixnum>=originalMap.Npix())
        {
            std::cerr << "Error: requesting pixnum too large for Healpix_Map, pixnum is " << pixnum << ", Npix is " << originalMap.Npix() << std::endl;
            throw std::invalid_argument("minkmapSphere::minkmapPixelNeighbors: pixnum too high");
        }
        else if(pixnum<0 && pixnum!=-5 && pixnum!= -11)
        {
            std::cerr << "Error: requesting negative undefined pixnum, pixnum is " << pixnum << ", Npix is " << originalMap.Npix() << std::endl;
            throw std::invalid_argument("minkmapSphere::minkmapPixelNeighbors: pixnum invalid");
        }
        
        if(pixnum==-5) //5OUTH pole
        {
            std::vector<int> southPolarCap;
            if(originalMap.Scheme()==RING)
            {
                int npix = originalMap.Npix();
                southPolarCap = {npix-4,npix-3,npix-2,npix-1};
            }
            else
            {
                int nsidesquared = originalMap.Nside()*originalMap.Nside();
                southPolarCap = {nsidesquared*8, nsidesquared*9, nsidesquared*10, nsidesquared*11};
            }
            return southPolarCap;
        }
        else if(pixnum==-11) //11ORTH pole
        {
            std::vector<int> northPolarCap;
            if(originalMap.Scheme()==RING)
            {
                northPolarCap = {0,1,2,3};
            }
            else
            {
                int nsidesquared = originalMap.Nside()*originalMap.Nside();
                northPolarCap = {nsidesquared-1, nsidesquared*2-1, nsidesquared*3-1, nsidesquared*4-1};
            }
            return northPolarCap;
        }
        
        fix_arr<int, 8> neighbors; //neighbors of this pixel
        originalMap.neighbors(pixnum,neighbors);
        std::vector<int> easternNeighborship{pixnum, neighbors[4],neighbors[5],neighbors[6]}; //neighbors east of this pixel and this pixel
        return easternNeighborship;
    }
    
    ///Returns linear combination Minkowski tensors at given Minkmap-pixel (vertex east of Healpix pixel with same number) 
    minkTensorStack at(int pixnum) const override
    {
        auto neighbors = minkmapPixelNeighbors(pixnum);
        //calculate one marching square/triangle
        return integrateMinktensor(neighbors);
    }
    
    
    ///Returns MT at marching square defined by surrounding Healpix pixels given by neighborship
    minkTensorStack integrateMinktensor(std::vector<int>& neighborship) const
    {
        //marching square (which above, below thresh)
        std::vector<double> values;
        double area=0;              //Use these for adding everything up, need to think about curvature TODO
        double length=0;
        double curvature=0;
        double weight=0.25; //most squares appear in 4 actual pixels, need to divide by this
        
        for (int pixnum : neighborship)//values-creation-loop
        {
            if(pixnum!=-1) values.push_back(originalMap[pixnum]);
        }
        const uint valuesSize = values.size();
        
        uint caseindex = 0; //caseindex: every corner gets a position in a 4 bit number, bit set to 1 if corner>thresh 0 else. Gives number of each case. If diagonal above/below, check overall average to see whether connected
        //If one of the pixels is masked, i.e. equal to nan, treat the whole area as masked
        for(uint i=0;i<valuesSize;i++)
        {
            if(std::isnan(values[i])) {
                return nanTensor(rankA,rankB,curvIndex,pointing(1.5701963268,0));
            }
            if(values[i]>=thresh) caseindex += pow(2,i);
            
        }
        
        minkTensorStack integralNumbers(rankA,rankB,curvIndex,pointing(1.5701963268,0)); 
        
        if(valuesSize==3)
        {
            //do triangle things
            weight = 1.d/3.; //triangles appear in 3 final map pixels
            integralNumbers = std::move(threeCornerCases(neighborship, values, caseindex, area, length, curvature)*weight); 
        } else if (valuesSize==4)
        {
            //do 4 corner things
            integralNumbers = std::move(fourCornerCases(neighborship, values, caseindex, area, length, curvature)*weight);
        }
        else{
            std::cerr << "Error: neighborhood has neither 3 nor 4 corners! Number of corners: " << valuesSize << " , this makes no sense. Eastern end of neighborhood: px number " << neighborship.at(0) << std::endl;
            throw std::invalid_argument( "minkmapSphere::integrateMinktensor: Weird number of corners" );
        }
        
        return integralNumbers;
    }
    
    //corners: two interpolated corners such that corners[0] cross corners[1] points away from body, corners[2] is within triangle, corners[3] and corners[4] needed for giveCurv and should be given in appropriate order (outBodyA, outBodyB)
    minkTensorStack oneCornerOver(double& newlength, double& area, double& newcurv, const std::vector<pointing>& corners, bool ranksum) const
    {
        assert(corners.size()==5 && "minkmapSphere::oneCornerOver: corners-vector has wrong size!");
        
        
        newlength = arclength(corners[0],corners[1]);
        double newarea = sphereArea(corners[2],corners[0],corners[1]);
        area += newarea;
        double& factor{ curvIndex==0 ? newarea : curvIndex==1 ? newlength : newcurv };
        if(curvIndex>1)
        { //Find exterior angle between normal vector and hypothetical curve perpendicular to edge of cell, should be positive if convex
            vec3 dirN = crossprod(corners[0].to_vec3(), corners[1].to_vec3());
            dirN.Normalize();
            newcurv = giveCurv(dirN, corners[2], corners[3], corners[2], corners[4]);
        }
        if(ranksum)
        {
            pointing n = getN_rotation(corners[0],corners[1]);
            return (minkTensorIntegrand(rankA, rankB, curvIndex, corners[0], n)* factor);
        }
        return minkTensorIntegrand(rankA, rankB, curvIndex, corners[0], pointing(1,0))*factor; //return one for scalar case
    }
    
    //corners: two interpolated corners such that corners[0] cross corners[1] points away from body, corners[2] and [4] are within polytope, corners[3] and corners[5] needed for giveCurv and should be given in appropriate order (outBodyA, outBodyB)
    minkTensorStack twoCornersOver(double& newlength, double& area, double& newcurv, const std::vector<pointing>& corners, bool ranksum) const
    {
        assert(corners.size()==6 && "minkmapSphere::twoCornersOver: corners-vector has wrong size!");
        
        double& factor{ curvIndex==0 ? area : curvIndex==1 ? newlength : newcurv };
        
        newlength = arclength(corners[0],corners[1]);
        area += sphereArea(corners[4],corners[2],corners[0]);
        area += sphereArea(corners[2],corners[1],corners[0]);
        if(curvIndex>1)
        { //Find exterior angle between normal vector and hypothetical curve perpendicular to edge of cell, should be positive if convex
            vec3 dirN = crossprod(corners[0].to_vec3(), corners[1].to_vec3());
            dirN.Normalize();
            newcurv = giveCurv(dirN, corners[2], corners[3], corners[4], corners[5]);
        }
        if(ranksum)
        {
            pointing n = getN_rotation(corners[0],corners[1]);
            return (minkTensorIntegrand(rankA, rankB, curvIndex, corners[0], n)* factor);
        }
        return minkTensorIntegrand(rankA, rankB, curvIndex, corners[0], pointing(1,0))*factor; //return one for scalar case
    }
    
    minkTensorStack fourCornerCases(const std::vector<int>& neighborship,const std::vector<double>& values,const uint& caseindex, double& area, double& length, double& curvature) const
    {
        std::vector<pointing> positions;
        for(auto pixnum : neighborship) {positions.push_back(originalMap.pix2ang(pixnum));}
        
        pointing oneCorner;
        pointing otherCorner;
        pointing n(1,0);
        minkTensorStack returntens = nullTensor(rankA, rankB, curvIndex, positions[0]);
        double newlength{0}, newcurv{0}, mean; //newlength: one-corner cases can't just use *length because of cases 5 and 10 (would multiply return tensor with other segment also)
        
        /* Corner numeration: 
         * 
         *    1
         * 0     2
         *    3
         * 
         * Pay attention to direction of n! In getN(first, second) (first x second) should point away from body
         * caseindex: every corner gets a position in a 4 bit number, bit set to 1 if corner>thresh, 0 else
         */
        double& factor{ curvIndex==0 ? area : curvIndex==1 ? newlength : newcurv };
        
        const uint ranksum = rankA+rankB;
        switch (caseindex)
        {
            case 0: //nix
                break;
            case 1: //one corner (pixel 0)
                oneCorner = interpPointing(positions.at(0),values.at(0),positions.at(1),values.at(1),thresh);
                otherCorner = interpPointing(positions.at(0),values.at(0),positions.at(3),values.at(3),thresh);
                
                returntens = oneCornerOver(newlength, area, newcurv, {oneCorner, otherCorner, positions[0], positions[3], positions[1]}, (bool)ranksum);
                length += newlength;
                curvature += newcurv;
                
                break;
            case 2: //other corner (pixel 1)
                oneCorner = interpPointing(positions.at(1),values.at(1),positions.at(2),values.at(2),thresh);
                otherCorner = interpPointing(positions.at(1),values.at(1),positions.at(0),values.at(0),thresh);
                
                returntens = oneCornerOver(newlength, area, newcurv, {oneCorner, otherCorner, positions[1], positions[0], positions[2]}, (bool)ranksum);
                length += newlength;
                curvature += newcurv;
                
                break;
            case 3: //0 and 1 over
                oneCorner = interpPointing(positions.at(0),values.at(0),positions.at(3),values.at(3),thresh);
                otherCorner = interpPointing(positions.at(1),values.at(1),positions.at(2),values.at(2),thresh);
                
                returntens = twoCornersOver(newlength, area, newcurv, {otherCorner,oneCorner, positions.at(0), positions.at(3), positions.at(1), positions.at(2)}, (bool)ranksum);
                length += newlength;
                curvature += newcurv;
                
                break;
            case 4: //other corner (pixel 2)
                oneCorner = interpPointing(positions.at(2),values.at(2),positions.at(3),values.at(3),thresh);
                otherCorner = interpPointing(positions.at(1),values.at(1),positions.at(2),values.at(2),thresh);
                
                returntens = oneCornerOver(newlength, area, newcurv, {oneCorner, otherCorner, positions[2], positions[1], positions[3]}, (bool)ranksum);
                length += newlength;
                curvature += newcurv;
                
                break;
            case 5: //2 and 0 over, check if average above or below thresh and view as connected or not accordingly
                mean = std::accumulate(values.begin(), values.end(), 0.)/values.size();
                if (mean>thresh) //larger mean = connected hexagon
                {
                    oneCorner = interpPointing(positions.at(0),values.at(0),positions.at(1),values.at(1),thresh);
                    otherCorner = interpPointing(positions.at(0),values.at(0),positions.at(3),values.at(3),thresh);
                    pointing thirdCorner = interpPointing(positions.at(2),values.at(2),positions.at(1),values.at(1),thresh);
                    pointing fourthCorner = interpPointing(positions.at(2),values.at(2),positions.at(3),values.at(3),thresh);
                    
                    
                    area += sphereArea(positions.at(0),otherCorner,oneCorner);
                    area += sphereArea(positions.at(2),thirdCorner,fourthCorner);
                    area += sphereArea(thirdCorner,oneCorner,fourthCorner);
                    area += sphereArea(otherCorner,fourthCorner,oneCorner);
                    double length1 = arclength(oneCorner,thirdCorner);
                    double length2 = arclength(fourthCorner,otherCorner); //getN, integrieren
                    length += length1+length2;
                    newlength = length;
                    
                    pointing n1 = getN_rotation(oneCorner, thirdCorner);
                    pointing n2 = getN_rotation(fourthCorner,otherCorner);
                    n2 = parallelTransport(fourthCorner, oneCorner, n2); 
                    
                    double curv1{0}, curv2{0};
                    if(curvIndex>1)
                    { //Find exterior angle between normal vector and hypothetical curve perpendicular to edge of cell, should be positive if convex
                        vec3 dirN = crossprod(oneCorner.to_vec3(), thirdCorner.to_vec3());
                        dirN.Normalize();
                        curv1 = giveCurv(dirN, positions.at(2), positions.at(1), positions.at(0), positions.at(1));
                        curvature += curv1;
                        
                        dirN = crossprod(fourthCorner.to_vec3(), otherCorner.to_vec3());
                        dirN.Normalize();
                        curv2 = giveCurv(dirN, positions.at(0), positions.at(3), positions.at(2), positions.at(3));
                        curvature += curv2;
                        newcurv = curvature;
                        return (minkTensorStack(minkTensorIntegrand(rankA, rankB, curvIndex, oneCorner, n1),curv1) + minkTensorStack(minkTensorIntegrand(rankA, rankB, curvIndex, oneCorner, n2),curv2)); //each section separately, here with curvature
                    }
                    else if(curvIndex==1)
                    {
                        return (minkTensorStack(minkTensorIntegrand(rankA, rankB, curvIndex, oneCorner, n1),length1) + minkTensorStack(minkTensorIntegrand(rankA, rankB, curvIndex, oneCorner, n2),length2)); //each section separately, here with length only
                    }
                    else
                    {
                        return minkTensorStack(minkTensorIntegrand(rankA, rankB, curvIndex, oneCorner, n1),area); //only one section with area
                    }
                    
                }
                else //smaller mean = disconnected triangles
                {
                    return( fourCornerCases(neighborship,values,1,area,length,curvature) + fourCornerCases(neighborship,values,4,area,length,curvature) );
                }
                break;
            case 6: //2 and 1 over
                oneCorner = interpPointing(positions.at(2),values.at(2),positions.at(3),values.at(3),thresh);
                otherCorner = interpPointing(positions.at(1),values.at(1),positions.at(0),values.at(0),thresh);
                
                returntens = twoCornersOver(newlength, area, newcurv, {oneCorner,otherCorner, positions.at(1), positions.at(0), positions.at(2), positions.at(3)}, (bool)ranksum);
                length += newlength;
                curvature += newcurv;
                
                break; 
            case 7: //all except pixel 3
                oneCorner = interpPointing(positions.at(2),values.at(2),positions.at(3),values.at(3),thresh);
                otherCorner = interpPointing(positions.at(3),values.at(3),positions.at(0),values.at(0),thresh);
                newlength = arclength(oneCorner,otherCorner);
                length += newlength;
                area += sphereArea(positions.at(0),otherCorner,positions.at(1));
                area += sphereArea(positions.at(1),otherCorner,oneCorner);
                area += sphereArea(positions.at(2),positions.at(1),oneCorner);
                if(curvIndex>1)
                { //Find exterior angle between normal vector and hypothetical curve perpendicular to edge of cell, should be positive if convex
                    vec3 dirN = crossprod(oneCorner.to_vec3(), otherCorner.to_vec3());
                    dirN.Normalize();
                    newcurv = giveCurv(dirN, positions.at(0), positions.at(3), positions.at(2), positions.at(3));
                    curvature += newcurv;
                }
                if(ranksum)
                {
                    n = getN_rotation(oneCorner,otherCorner);
                }
                return (minkTensorIntegrand(rankA, rankB, curvIndex, oneCorner, n)* factor);
                break;
            case 8: //other corner (pixel 3)
                oneCorner = interpPointing(positions.at(2),values.at(2),positions.at(3),values.at(3),thresh);
                otherCorner = interpPointing(positions.at(3),values.at(3),positions.at(0),values.at(0),thresh);
                
                returntens = oneCornerOver(newlength, area, newcurv, {otherCorner, oneCorner, positions[3], positions[2], positions[0]}, (bool)ranksum);
                length += newlength;
                curvature += newcurv;
                
                break;
            case 9: //3 and 0 over
                oneCorner = interpPointing(positions.at(2),values.at(2),positions.at(3),values.at(3),thresh);
                otherCorner = interpPointing(positions.at(1),values.at(1),positions.at(0),values.at(0),thresh);
                
                returntens = twoCornersOver(newlength, area, newcurv, {otherCorner,oneCorner, positions.at(3), positions.at(2), positions.at(0), positions.at(1)}, (bool)ranksum);
                length += newlength;
                curvature += newcurv;
                
                break;
            case 10: //3 and 1 over, check if average above or below thresh and view as connected or not accordingly
                mean = std::accumulate(values.begin(), values.end(), 0.)/values.size();
                if(mean>thresh) //larger mean = connected hexagon
                {
                    oneCorner = interpPointing(positions.at(0),values.at(0),positions.at(1),values.at(1),thresh);
                    otherCorner = interpPointing(positions.at(0),values.at(0),positions.at(3),values.at(3),thresh);
                    pointing thirdCorner = interpPointing(positions.at(2),values.at(2),positions.at(1),values.at(1),thresh);
                    pointing fourthCorner = interpPointing(positions.at(2),values.at(2),positions.at(3),values.at(3),thresh);
                    area += sphereArea(positions.at(1),oneCorner,thirdCorner);
                    area += sphereArea(positions.at(3),fourthCorner,otherCorner);
                    area += sphereArea(thirdCorner,oneCorner,fourthCorner);
                    area += sphereArea(otherCorner,fourthCorner,oneCorner);
                    double length1 = arclength(otherCorner,oneCorner);
                    double length2 = arclength(thirdCorner,fourthCorner); //getN, integrieren
                    length += length1+length2;
                    
                    pointing n1 = getN_rotation(otherCorner,oneCorner);
                    pointing n2 = getN_rotation(thirdCorner,fourthCorner);
                    n2 = parallelTransport(thirdCorner, otherCorner, n2);
                    
                    double curv1{0}, curv2{0};
                    if(curvIndex>1)
                    { //Find exterior angle between normal vector and hypothetical curve perpendicular to edge of cell, should be positive if convex
                        vec3 dirN = crossprod(otherCorner.to_vec3(), oneCorner.to_vec3());
                        dirN.Normalize();
                        curv1 = giveCurv(dirN, positions.at(1), positions.at(0), positions.at(3), positions.at(0));
                        curvature += curv1;
                        
                        dirN = crossprod(thirdCorner.to_vec3(), fourthCorner.to_vec3());
                        dirN.Normalize();
                        curv2 = giveCurv(dirN, positions.at(3), positions.at(2), positions.at(1), positions.at(2));
                        curvature += curv2;
                        return (minkTensorStack(minkTensorIntegrand(rankA, rankB, curvIndex, oneCorner, n1),curv1) + minkTensorStack(minkTensorIntegrand(rankA, rankB, curvIndex, oneCorner, n2),curv2)); //each section separately, here with curvature
                    }
                    else if(curvIndex==1)
                    {
                        return(minkTensorIntegrand(rankA, rankB, curvIndex, otherCorner, n1)*length1 + minkTensorIntegrand(rankA, rankB, curvIndex, otherCorner, n2)*length2 ); //each section separately, here with length
                    }
                    else
                    {
                        return minkTensorIntegrand(rankA, rankB, curvIndex, otherCorner, n1)*area; //only one section with area
                    }
                }
                else //smaller mean = disconnected triangles
                {
                    return( fourCornerCases(neighborship,values,8,area,length,curvature) + fourCornerCases(neighborship,values,2,area,length,curvature) );
                }
                break;
            case 11: //all except pixel 2
                oneCorner = interpPointing(positions.at(2),values.at(2),positions.at(3),values.at(3),thresh);
                otherCorner = interpPointing(positions.at(1),values.at(1),positions.at(2),values.at(2),thresh);
                newlength = arclength(oneCorner,otherCorner);
                length += newlength;
                area += sphereArea(positions.at(3),oneCorner,positions.at(0));
                area += sphereArea(positions.at(0),oneCorner,otherCorner);
                area += sphereArea(positions.at(1),positions.at(0),otherCorner);
                if(curvIndex>1)
                { //Find exterior angle between normal vector and hypothetical curve perpendicular to edge of cell, should be positive if convex
                    vec3 dirN = crossprod(otherCorner.to_vec3(), oneCorner.to_vec3());
                    dirN.Normalize();
                    newcurv = giveCurv(dirN, positions.at(3), positions.at(2), positions.at(1), positions.at(2));
                    curvature += newcurv;
                }
                if(ranksum)
                {
                    n = getN_rotation(otherCorner,oneCorner);
                }
                return (minkTensorIntegrand(rankA, rankB, curvIndex, otherCorner, n)* factor);
                break;
            case 12: //2 and 3 over
                oneCorner = interpPointing(positions.at(0),values.at(0),positions.at(3),values.at(3),thresh);
                otherCorner = interpPointing(positions.at(1),values.at(1),positions.at(2),values.at(2),thresh);
                
                returntens = twoCornersOver(newlength, area, newcurv, {oneCorner,otherCorner, positions.at(2), positions.at(1), positions.at(3), positions.at(0)}, (bool)ranksum);
                length += newlength;
                curvature += newcurv;
                
                break;
            case 13: //all except pixel 1
                oneCorner = interpPointing(positions.at(1),values.at(1),positions.at(2),values.at(2),thresh);
                otherCorner = interpPointing(positions.at(1),values.at(1),positions.at(0),values.at(0),thresh);
                newlength = arclength(oneCorner,otherCorner);
                length += newlength;
                area += sphereArea(positions.at(3),positions.at(2),oneCorner);
                area += sphereArea(positions.at(3),oneCorner,otherCorner);
                area += sphereArea(positions.at(3),otherCorner,positions.at(0));
                if(curvIndex>1)
                { //Find exterior angle between normal vector and hypothetical curve perpendicular to edge of cell, should be positive if convex
                    vec3 dirN = crossprod(otherCorner.to_vec3(), oneCorner.to_vec3());
                    dirN.Normalize();
                    newcurv = giveCurv(dirN, positions.at(2), positions.at(1), positions.at(0), positions.at(1));
                    curvature += newcurv;
                }
                if(ranksum)
                {
                    n = getN_rotation(otherCorner,oneCorner);
                }
                return (minkTensorIntegrand(rankA, rankB, curvIndex, otherCorner, n)* factor);
                break;
            case 14: //all except pixel 0
                oneCorner = interpPointing(positions.at(0),values.at(0),positions.at(1),values.at(1),thresh);
                otherCorner = interpPointing(positions.at(0),values.at(0),positions.at(3),values.at(3),thresh);
                newlength = arclength(oneCorner,otherCorner);
                length += newlength;
                area += sphereArea(positions.at(2),positions.at(1),oneCorner);
                area += sphereArea(positions.at(2),oneCorner,otherCorner);
                area += sphereArea(positions.at(2),otherCorner,positions.at(3));
                if(curvIndex>1)
                { //Find exterior angle between normal vector and hypothetical curve perpendicular to edge of cell, should be positive if convex
                    vec3 dirN = crossprod(otherCorner.to_vec3(), oneCorner.to_vec3());
                    dirN.Normalize();
                    newcurv = giveCurv(dirN, positions.at(1), positions.at(0), positions.at(3), positions.at(0));
                    curvature += newcurv;
                }
                if(ranksum)
                {
                    n = getN_rotation(otherCorner,oneCorner);
                }
                return (minkTensorIntegrand(rankA, rankB, curvIndex, otherCorner, n)* factor);
                break;
            case 15: //alles
                area += sphereArea(positions.at(0),positions.at(2),positions.at(1));
                area += sphereArea(positions.at(0),positions.at(3),positions.at(2));
                if(curvIndex==0)
                {
                    return minkTensorIntegrand(rankA, rankB, curvIndex, positions.at(0),n)*area;
                }
                break;
            default:
                std::cerr << "Error: invalid case number: " << caseindex << " , this makes no sense. Eastern end of neighborhood: px number " << neighborship.at(0) << std::endl;
                throw std::invalid_argument( "minkmapSphere: Weird caseindex" );
        }

        return returntens;
    }
    
    minkTensorStack threeCornerCases(std::vector<int>& neighborship, std::vector<double>& values, uint caseindex, double& area, double& length, double& curvature) const
    {
        std::vector<pointing> positions;
        
        for(auto pixnum : neighborship) 
        {
            if(pixnum!=-1)  positions.push_back(originalMap.pix2ang(pixnum));
        }
        pointing oneCorner;
        pointing otherCorner;
        //pointing n;
        minkTensorStack returntens = nullTensor(rankA, rankB, curvIndex, positions[0]);
        double newlength{0}, newcurv{0};
        
        //double& factor{ curvIndex==0 ? area : curvIndex==1 ? newlength : newcurv };
        
        /* Corner numeration:
         * 
         *     1
         * 0
         *      2
         * 
         */
        //tensor2D theTensor(rankA, rankB, curvIndex); //TODO hier integrieren
        
        const uint ranksum = rankA+rankB;
        switch (caseindex)
        {
            case 0: //nix
                break;
            case 1: //one corner (pixel 0)
                oneCorner = interpPointing(positions.at(0),values.at(0),positions.at(1),values.at(1),thresh);
                otherCorner = interpPointing(positions.at(0),values.at(0),positions.at(2),values.at(2),thresh);
                
                returntens = oneCornerOver(newlength, area, newcurv, {oneCorner, otherCorner, positions[0], positions[2], positions[1]}, (bool)ranksum);
                length += newlength;
                curvature += newcurv;
                
                break;
            case 2: //other corner (pixel 1)
                oneCorner = interpPointing(positions.at(0),values.at(0),positions.at(1),values.at(1),thresh);
                otherCorner = interpPointing(positions.at(1),values.at(1),positions.at(2),values.at(2),thresh);
                
                returntens = oneCornerOver(newlength, area, newcurv, {otherCorner, oneCorner, positions[1], positions[0], positions[2]}, (bool)ranksum);
                length += newlength;
                curvature += newcurv;
                
                break;
            case 3: //all except pixel 2
                oneCorner = interpPointing(positions.at(0),values.at(0),positions.at(2),values.at(2),thresh);
                otherCorner = interpPointing(positions.at(1),values.at(1),positions.at(2),values.at(2),thresh);
                returntens = twoCornersOver(newlength, area, newcurv, {otherCorner,oneCorner, positions.at(0), positions.at(2), positions.at(1), positions.at(2)}, (bool)ranksum);
                length += newlength;
                curvature += newcurv;
                
                break;
            case 4: //other corner (pixel 2)
                oneCorner = interpPointing(positions.at(0),values.at(0),positions.at(2),values.at(2),thresh);
                otherCorner = interpPointing(positions.at(1),values.at(1),positions.at(2),values.at(2),thresh);
                                
                returntens = oneCornerOver(newlength, area, newcurv, {oneCorner, otherCorner, positions[2], positions[1], positions[0]}, (bool)ranksum);
                length += newlength;
                curvature += newcurv;
                
                break;
            case 5: //all except pixel 1
                oneCorner = interpPointing(positions.at(0),values.at(0),positions.at(1),values.at(1),thresh);
                otherCorner = interpPointing(positions.at(1),values.at(1),positions.at(2),values.at(2),thresh);
                returntens = twoCornersOver(newlength, area, newcurv, {oneCorner,otherCorner, positions.at(2), positions.at(1), positions.at(0), positions.at(1)}, (bool)ranksum);
                length += newlength;
                curvature += newcurv;
                
                break;
            case 6: //all except pixel 0
                oneCorner = interpPointing(positions.at(0),values.at(0),positions.at(1),values.at(1),thresh);
                otherCorner = interpPointing(positions.at(0),values.at(0),positions.at(2),values.at(2),thresh);
                returntens = twoCornersOver(newlength, area, newcurv, {otherCorner,oneCorner, positions.at(1), positions.at(0), positions.at(2), positions.at(0)}, (bool)ranksum);
                length += newlength;
                curvature += newcurv;
                
                break;
            case 7: //alles
                area += sphereArea(positions.at(0),positions.at(2),positions.at(1));
                if(curvIndex==0)
                {
                    return minkTensorIntegrand(rankA, rankB, curvIndex, positions.at(0),pointing(1,0))*area;
                }
                break;
            default:
                std::cerr << "Error: invalid case number: " << caseindex << " , this makes no sense. Eastern end of neighborhood: px number " << neighborship.at(0) << std::endl;
                throw std::invalid_argument( "minkmapSphere: Weird caseindex" );
        }
        return returntens;
    }
    
    /*
     * Gives total exterior angle for one cell. Angle is calculated to hypothetical contour continuing perpendicular to cell wall
     * pointings should be given such that inBodyA (cross) outBodyA and outBodyB (cross) inBodyB points into cell ("right side first")
     */
    double giveCurv(vec3 dirAwayFromBody, const pointing& inBodyA, const pointing& outBodyA, const pointing& inBodyB, const pointing& outBodyB) const
    {
        double retCurv = 0;
        vec3 dir03 = crossprod(inBodyA.to_vec3(), outBodyA.to_vec3());
        dir03.Normalize();
        double ang1 = asin(dotprod(dirAwayFromBody,dir03)); //arcsin because dir03 is perpendicular to direction whose angle we want
        retCurv += ang1;
        
        vec3 dir10 = crossprod(outBodyB.to_vec3(), inBodyB.to_vec3()); //other way round to take care of sign, should point towards center of cell
        dir10.Normalize();
        double ang2 = asin(dotprod(dirAwayFromBody,dir10)); //arcsin because dir10 is perpendicular to direction whose angle we want
        retCurv += ang2;
        return retCurv;
    }
    
};




#endif
