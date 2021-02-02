#ifndef litchi_pulp
#define litchi_pulp

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>

#include "healpix_cxx/healpix_map.h"
#include "litchi_kernel.hpp"

#include "tensorFamily.hpp"
#include "tensor2D.hpp"
#include "tensorOperations.hpp"
#include "minkTensorIntegrand.hpp"
#include "geometryHelpers.hpp"


extern const double pi;



struct minkmapSphere :  minkmapFamily{ //should contain "raw" (marching-square-level) minkowski tensor functions. Don't save all because of space. 
    const double thresh;
    
    minkmapSphere(Healpix_Map<double>& map, uint rank1 = 0, uint rank2 = 0, uint curvind = 0, double threshold = 0) : minkmapFamily(map, rank1, rank2, curvind), thresh(threshold) {
        if (curvind > 2) {
            std::cerr << "Error: invalid curvIndex: " << curvind << " , this makes no sense. In 2D only up to two!" << std::endl;
            throw std::invalid_argument( "minkmapSphere: Weird curvIndex" );
        }
    }
    
    pointing pix2ang(int pixnum) const //TODO poles
    {
        std::vector<vec3> corners;
        originalMap.boundaries(pixnum, 1, corners); //find corners of original pixel in original map
        return pointing(corners.at(3)); //position east of pixel is center of vertex
    }
    
    //template <typename tensortype>
    tensor2D at(int pixnum) const
    {
        if(pixnum>=originalMap.Npix())
        {
            std::cerr << "Error: requesting pixnum too large for Healpix_Map, pixnum is " << pixnum << ", Npix is " << originalMap.Npix() << std::endl;
            throw std::invalid_argument("minkmapSphere::at: pixnum too high");
        }
        else if(pixnum<0 && pixnum!=-5 && pixnum!= -11)
        {
            std::cerr << "Error: requesting negative undefined pixnum, pixnum is " << pixnum << ", Npix is " << originalMap.Npix() << std::endl;
            throw std::invalid_argument("minkmapSphere::at: pixnum invalid");
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
            return integrateMinktensor(southPolarCap);
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
            return integrateMinktensor(northPolarCap);
        }
        
        fix_arr<int, 8> neighbors; //neighbors of this pixel
        originalMap.neighbors(pixnum,neighbors);
        std::vector<int> easternNeighborship{pixnum, neighbors[4],neighbors[5],neighbors[6]}; //neighbors east of this pixel and this pixel
        //calculate one marching square/triangle
        tensor2D tensorHere = integrateMinktensor(easternNeighborship);
        return tensorHere;
    }
    
    
    //template <typename tensortype>
    tensor2D integrateMinktensor(std::vector<int>& neighborship) const
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
        uint valuesSize = values.size();
        
        uint caseindex = 0; //Number of case (pattern above/below thresh). If diagonal above/below, check overall average to see whether connected
        for(uint i=0;i<valuesSize;i++)
        {
            if(values.at(i)>=thresh) caseindex += pow(2,i);
        }
        
        //Check if 3 or 4 corners and check every case
        tensor2D integralNumbers(rankA, rankB, curvIndex); //TODO evtl dont use tensor2d
        if(valuesSize==3)
        {
            //do triangle things
            weight = 1.d/3.; //triangles appear in 3 final map pixels
            integralNumbers = threeCornerCases(neighborship, values, caseindex, area, length, curvature)*weight; 
        } else if (valuesSize==4)
        {
            //do 4 corner things
            integralNumbers = fourCornerCases(neighborship, values, caseindex, area, length, curvature)*weight;
        }
        else{
            std::cerr << "Error: neighborhood has neither 3 nor 4 corners! Number of corners: " << valuesSize << " , this makes no sense. Eastern end of neighborhood: px number " << neighborship.at(0) << std::endl;
            throw std::invalid_argument( "minkmapSphere::integrateMinktensor: Weird number of corners" );
        }
        
        
        //define minkTensorIntegrand (find normal vector), curvature for good point
        std::vector<vec3> corners;
        originalMap.boundaries(neighborship.at(0), 1, corners); //find corners of original pixel in original map
        pointing n(pi/2,0);
        pointing r(corners.at(3)); //position east of pixel is center of vertex
        
        //do the integration (parallel transport, multiply with length)
        if(rankA+rankB == 0) //scalar functional
        {
            switch (curvIndex)
            {
                case 0:
                    return tensor2D(area*weight,curvIndex);
                case 1:
                    return tensor2D(length*weight,curvIndex);
                case 2:
                    std::cout << "Warning! Curvature functional not implemented in integrateMinktensor!" << std::endl;
                    return tensor2D(curvature*weight,curvIndex);
                default:
                    std::cerr << "Error: invalid curvIndex: " << curvIndex << " , this makes no sense. In 2D only up to two!" << std::endl;
                    throw std::invalid_argument( "minkmapSphere::integrateMinktensor: Weird curvIndex" );
            }
                    
        }
        return integralNumbers;
    }
    
    //template<typename tensortype>
    minkTensorIntegrand fourCornerCases(const std::vector<int>& neighborship,const std::vector<double>& values,const uint& caseindex, double& area, double& length, double& curvature) const
    {
        std::vector<pointing> positions;
        for(auto pixnum : neighborship) {positions.push_back(originalMap.pix2ang(pixnum));}
        pointing oneCorner;
        pointing otherCorner;
        pointing n;
        double mean;
        minkTensorIntegrand theTensor(rankA, rankB, curvIndex); //TODO hier integrieren
        
        /* Corner numeration: 
         * 
         *    1
         * 0     2
         *    3
         * 
         * Pay attention to direction of n! In getN(first, second) (first x second) should point away from body
         */
        
        uint ranksum = rankA+rankB;
        
        switch (caseindex)
        {
            case 0: //nix
                break;
            case 1: //one corner (pixel 0)
                oneCorner = interpPointing(positions.at(0),values.at(0),positions.at(1),values.at(1),thresh);
                otherCorner = interpPointing(positions.at(0),values.at(0),positions.at(3),values.at(3),thresh);
                if(ranksum)
                {
                    n = getN_cartesian(oneCorner,otherCorner);
                }
                area += sphereArea(positions.at(0),oneCorner,otherCorner);
                length += arclength(oneCorner,otherCorner);
                break;
            case 2: //other corner (pixel 1)
                oneCorner = interpPointing(positions.at(1),values.at(1),positions.at(2),values.at(2),thresh);
                otherCorner = interpPointing(positions.at(1),values.at(1),positions.at(0),values.at(0),thresh);
                if(ranksum)
                {
                    n = getN_cartesian(oneCorner,otherCorner);
                }
                area += sphereArea(positions.at(1),oneCorner,otherCorner);
                length += arclength(oneCorner,otherCorner);
                break;
            case 3: //0 and 1 over
                oneCorner = interpPointing(positions.at(0),values.at(0),positions.at(3),values.at(3),thresh);
                otherCorner = interpPointing(positions.at(1),values.at(1),positions.at(2),values.at(2),thresh);
                if(ranksum)
                {
                    n = getN_cartesian(otherCorner,oneCorner);
                }
                area += sphereArea(positions.at(0),oneCorner,positions.at(1));
                area += sphereArea(positions.at(1),oneCorner,otherCorner);
                length += arclength(oneCorner,otherCorner);
                break;
            case 4: //other corner (pixel 2)
                oneCorner = interpPointing(positions.at(2),values.at(2),positions.at(3),values.at(3),thresh);
                otherCorner = interpPointing(positions.at(1),values.at(1),positions.at(2),values.at(2),thresh);
                if(ranksum)
                {
                    n = getN_cartesian(oneCorner,otherCorner);
                }
                area += sphereArea(positions.at(2),oneCorner,otherCorner);
                length += arclength(oneCorner,otherCorner);
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
                    length += arclength(oneCorner,thirdCorner);
                    length += arclength(otherCorner,fourthCorner); //TODO getN, integrieren
                }
                else //smaller mean = disconnected triangles
                {
                    auto theSum = fourCornerCases(neighborship,values,1,area,length,curvature) + fourCornerCases(neighborship,values,4,area,length,curvature);
                }
                break; 
            case 6: //2 and 1 over
                oneCorner = interpPointing(positions.at(2),values.at(2),positions.at(3),values.at(3),thresh);
                otherCorner = interpPointing(positions.at(1),values.at(1),positions.at(0),values.at(0),thresh);
                if(ranksum)
                {
                    n = getN_cartesian(oneCorner,otherCorner);
                }
                area += sphereArea(positions.at(2),positions.at(1),oneCorner);
                area += sphereArea(positions.at(1),otherCorner,oneCorner);
                length += arclength(oneCorner,otherCorner);
                break; 
            case 7: //all except pixel 3
                oneCorner = interpPointing(positions.at(2),values.at(2),positions.at(3),values.at(3),thresh);
                otherCorner = interpPointing(positions.at(3),values.at(3),positions.at(0),values.at(0),thresh);
                if(ranksum)
                {
                    n = getN_cartesian(oneCorner,otherCorner);
                }
                area += sphereArea(positions.at(0),otherCorner,positions.at(1));
                area += sphereArea(positions.at(1),otherCorner,oneCorner);
                area += sphereArea(positions.at(2),positions.at(1),oneCorner);
                length += arclength(oneCorner,otherCorner);
                break;
            case 8: //other corner (pixel 3)
                oneCorner = interpPointing(positions.at(2),values.at(2),positions.at(3),values.at(3),thresh);
                otherCorner = interpPointing(positions.at(3),values.at(3),positions.at(0),values.at(0),thresh);
                if(ranksum)
                {
                    n = getN_cartesian(otherCorner,oneCorner);
                }
                area += sphereArea(positions.at(3),oneCorner,otherCorner);
                length += arclength(oneCorner,otherCorner);
                break;
            case 9: //3 and 0 over
                oneCorner = interpPointing(positions.at(2),values.at(2),positions.at(3),values.at(3),thresh);
                otherCorner = interpPointing(positions.at(1),values.at(1),positions.at(0),values.at(0),thresh);
                if(ranksum)
                {
                    n = getN_cartesian(otherCorner,oneCorner);
                }
                area += sphereArea(positions.at(3),oneCorner,otherCorner);
                area += sphereArea(positions.at(0),positions.at(3),otherCorner);
                length += arclength(oneCorner,otherCorner);
                break;
            case 10: //3 and 1 over, check if average above or below thresh and view as connected or not accordingly
                mean = std::accumulate(values.begin(), values.end(), 0.)/values.size();
                if(mean>thresh) //larger mean = connected hectagon
                {
                    oneCorner = interpPointing(positions.at(0),values.at(0),positions.at(1),values.at(1),thresh);
                    otherCorner = interpPointing(positions.at(0),values.at(0),positions.at(3),values.at(3),thresh);
                    pointing thirdCorner = interpPointing(positions.at(2),values.at(2),positions.at(1),values.at(1),thresh);
                    pointing fourthCorner = interpPointing(positions.at(2),values.at(2),positions.at(3),values.at(3),thresh);
                    area += sphereArea(positions.at(1),oneCorner,thirdCorner);
                    area += sphereArea(positions.at(3),fourthCorner,otherCorner);
                    area += sphereArea(thirdCorner,oneCorner,fourthCorner);
                    area += sphereArea(otherCorner,fourthCorner,oneCorner);
                    length += arclength(oneCorner,thirdCorner);
                    length += arclength(otherCorner,fourthCorner); //TODO getN, integrieren
                }
                else //smaller mean = disconnected triangles
                {
                    auto theSum = fourCornerCases(neighborship,values,8,area,length,curvature) + fourCornerCases(neighborship,values,2,area,length,curvature);
                }
                break;
            case 11: //all except pixel 2
                oneCorner = interpPointing(positions.at(2),values.at(2),positions.at(3),values.at(3),thresh);
                otherCorner = interpPointing(positions.at(1),values.at(1),positions.at(2),values.at(2),thresh);
                if(ranksum)
                {
                    n = getN_cartesian(otherCorner,oneCorner);
                }
                area += sphereArea(positions.at(3),oneCorner,positions.at(0));
                area += sphereArea(positions.at(0),oneCorner,otherCorner);
                area += sphereArea(positions.at(1),positions.at(0),otherCorner);
                length += arclength(oneCorner,otherCorner);
                break;
            case 12: //2 and 3 over
                oneCorner = interpPointing(positions.at(0),values.at(0),positions.at(3),values.at(3),thresh);
                otherCorner = interpPointing(positions.at(1),values.at(1),positions.at(2),values.at(2),thresh);
                if(ranksum)
                {
                    n = getN_cartesian(oneCorner,otherCorner);
                }
                area += sphereArea(positions.at(3),otherCorner,oneCorner);
                area += sphereArea(positions.at(3),positions.at(2),otherCorner);
                length += arclength(oneCorner,otherCorner);
                break;
            case 13: //all except pixel 1
                oneCorner = interpPointing(positions.at(1),values.at(1),positions.at(2),values.at(2),thresh);
                otherCorner = interpPointing(positions.at(1),values.at(1),positions.at(0),values.at(0),thresh);
                if(ranksum)
                {
                    n = getN_cartesian(otherCorner,oneCorner);
                }
                area += sphereArea(positions.at(3),positions.at(2),oneCorner);
                area += sphereArea(positions.at(3),oneCorner,otherCorner);
                area += sphereArea(positions.at(3),otherCorner,positions.at(0));
                length += arclength(oneCorner,otherCorner);
                break;
            case 14: //all except pixel 0
                oneCorner = interpPointing(positions.at(0),values.at(0),positions.at(1),values.at(1),thresh);
                otherCorner = interpPointing(positions.at(0),values.at(0),positions.at(3),values.at(3),thresh);
                if(ranksum)
                {
                    n = getN_cartesian(otherCorner,oneCorner);
                }
                area += sphereArea(positions.at(2),positions.at(1),oneCorner);
                area += sphereArea(positions.at(2),oneCorner,otherCorner);
                area += sphereArea(positions.at(2),otherCorner,positions.at(3));
                length += arclength(oneCorner,otherCorner);
                break;
            case 15: //alles
                area += sphereArea(positions.at(0),positions.at(2),positions.at(1));
                area += sphereArea(positions.at(0),positions.at(3),positions.at(2));
                break;
            default:
                std::cerr << "Error: invalid case number: " << caseindex << " , this makes no sense. Eastern end of neighborhood: px number " << neighborship.at(0) << std::endl;
                throw std::invalid_argument( "minkmapSphere: Weird caseindex" );
        }
        
        return theTensor;
    }
    
    minkTensorIntegrand threeCornerCases(std::vector<int>& neighborship, std::vector<double>& values, uint caseindex, double& area, double& length, double& curvature) const
    {
        std::vector<pointing> positions;
        
        for(auto pixnum : neighborship) 
        {
            if(pixnum!=-1)  positions.push_back(originalMap.pix2ang(pixnum));
        }
        pointing oneCorner;
        pointing otherCorner;
        pointing n;
        
        /* Corner numeration: 
         * 
         *     1
         * 0
         *      2
         * 
         */
        
        uint ranksum = rankA+rankB;
        switch (caseindex)
        {
            case 0: //nix
                break;
            case 1: //one corner (pixel 0)
                oneCorner = interpPointing(positions.at(0),values.at(0),positions.at(1),values.at(1),thresh);
                otherCorner = interpPointing(positions.at(0),values.at(0),positions.at(2),values.at(2),thresh);
                if(ranksum)
                {
                    n = getN_cartesian(oneCorner,otherCorner);
                }
                area += sphereArea(positions.at(0),otherCorner,oneCorner);
                length += arclength(oneCorner,otherCorner);
                break;
            case 2: //other corner (pixel 1)
                oneCorner = interpPointing(positions.at(0),values.at(0),positions.at(1),values.at(1),thresh);
                otherCorner = interpPointing(positions.at(1),values.at(1),positions.at(2),values.at(2),thresh);
                if(ranksum)
                {
                    n = getN_cartesian(otherCorner,oneCorner);
                }
                area += sphereArea(positions.at(1),oneCorner,otherCorner);
                length += arclength(oneCorner,otherCorner);
                break;
            case 3: //all except pixel 2
                oneCorner = interpPointing(positions.at(0),values.at(0),positions.at(2),values.at(2),thresh);
                otherCorner = interpPointing(positions.at(1),values.at(1),positions.at(2),values.at(2),thresh);
                if(ranksum)
                {
                    n = getN_cartesian(otherCorner,oneCorner);
                }
                area += sphereArea(positions.at(1),oneCorner,otherCorner);
                area += sphereArea(positions.at(1),positions.at(0),oneCorner);
                length += arclength(oneCorner,otherCorner);
                break;
            case 4: //other corner (pixel 2)
                oneCorner = interpPointing(positions.at(0),values.at(0),positions.at(2),values.at(2),thresh);
                otherCorner = interpPointing(positions.at(1),values.at(1),positions.at(2),values.at(2),thresh);
                if(ranksum)
                {
                    n = getN_cartesian(oneCorner,otherCorner);
                }
                area += sphereArea(positions.at(2),otherCorner,oneCorner);
                length += arclength(oneCorner,otherCorner);
                break;
            case 5: //all except pixel 1
                oneCorner = interpPointing(positions.at(0),values.at(0),positions.at(1),values.at(1),thresh);
                otherCorner = interpPointing(positions.at(1),values.at(1),positions.at(2),values.at(2),thresh);
                if(ranksum)
                {
                    n = getN_cartesian(oneCorner,otherCorner);
                }
                area += sphereArea(positions.at(2),otherCorner,oneCorner);
                area += sphereArea(positions.at(2),oneCorner,positions.at(0));
                length += arclength(oneCorner,otherCorner);
                break;
            case 6: //all except pixel 0
                oneCorner = interpPointing(positions.at(0),values.at(0),positions.at(1),values.at(1),thresh);
                otherCorner = interpPointing(positions.at(0),values.at(0),positions.at(2),values.at(2),thresh);
                if(ranksum)
                {
                    n = getN_cartesian(otherCorner,oneCorner);
                }
                area += sphereArea(positions.at(1),oneCorner,otherCorner);
                area += sphereArea(positions.at(1),otherCorner,positions.at(2));
                length += arclength(oneCorner,otherCorner);
                break;
            case 7: //alles
                area += sphereArea(positions.at(0),positions.at(2),positions.at(1));
                break;
            default:
                std::cerr << "Error: invalid case number: " << caseindex << " , this makes no sense. Eastern end of neighborhood: px number " << neighborship.at(0) << std::endl;
                throw std::invalid_argument( "minkmapSphere: Weird caseindex" );
        }
        minkTensorIntegrand theTensor(rankA, rankB, curvIndex); //TODO hier integrieren
        return theTensor;
    }
    
    
};




#endif
