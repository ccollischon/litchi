#ifndef minkmapSphere
#define minkmapSphere

struct minkmapSphere { //should contain "raw" (marching-square-level) minkowski tensor functions. Don't save all because of space. 
    Healpix_Map<double>& originalMap;
    uint rankA, rankB, curvIndex;
    double thresh;
    
    minkmapSphere(Healpix_Map<double>& map) : originalMap(map) {}
    
    minkmapSphere(Healpix_Map<double>& map, uint rank1, uint rank2, uint curvind, double threshold) : originalMap(map), rankA(rank1), rankB(rank2), curvIndex(curvind), thresh(threshold) {}
    
    //template <typename tensortype>
    tensor2D integrateMinktensor(std::vector<int>& neighborship) const
    {
        //marching square (which above, below thresh)
        std::vector<double> values; 
        std::for_each(neighborship.begin(), neighborship.end(), [&](const auto& pixnum) {  //add value of neighborhood pixel to values-vector (-1 in neighbors() means no neighbor in this direction)
                                                                    if(pixnum!=-1) values.push_back(originalMap[pixnum]);
                                                                    } );
        uint valuesSize = values.size();
        
        uint caseindex = 0; //Number of case (pattern above/below thresh). If diagonal above/below, check overall average to see whether connected
        for(uint i=0;i<valuesSize;i++)
        {
            if(values.at(i)>=thresh) caseindex += pow(2,i);
        }
        
        //Check if 3 or 4 corners and check every case        
        std::cout << "Warning! Area not properly implemented!" << std::endl;
        minkTensorIntegrand integrand(rankA, rankB, curvIndex);
        if(valuesSize==3)
        {
            //do triangle things
            integrand = threeCornerCases(neighborship, values, caseindex); //TODO needs probably edge lengths or something that come from case analysis
        } else if (valuesSize==4)
        {
            //do 4 corner things
            integrand = fourCornerCases(neighborship, values, caseindex);
        }
        else{
            std::cerr << "Error: neighborhood has neither 3 nor 4 corners! Number of corners: " << valuesSize << " , this makes no sense. Eastern end of neighborhood: px number " << neighborship.at(0) << std::endl;
            throw std::invalid_argument( "minkmapSphere: Weird number of corners" );
        }
        
        
        //define minkTensorIntegrand (find normal vector), curvature for good point
        std::vector<vec3> corners;
        originalMap.boundaries(neighborship.at(0), 1, corners); //find corners of original pixel in original map
        pointing n(pi/2,0);
        pointing r(corners.at(3)); //position east of pixel is center of vertex
        
        //do the integration (parallel transport, multiply with length)
        return tensor2D(0,0,0);
    }
    
    minkTensorIntegrand fourCornerCases(std::vector<int>& neighborship, std::vector<double>& values, uint caseindex) const
    {
        switch (caseindex)
        {
            case 0: //nix
                break;
            case 15: //alles
                break; //TODO area
            default:
                std::cerr << "Error: invalid case number: " << caseindex << " , this makes no sense. Eastern end of neighborhood: px number " << neighborship.at(0) << std::endl;
                throw std::invalid_argument( "minkmapSphere: Weird caseindex" );
        }
        minkTensorIntegrand theTensor(rankA, rankB, curvIndex);
        return theTensor;
    }
    
    minkTensorIntegrand threeCornerCases(std::vector<int>& neighborship, std::vector<double>& values, uint caseindex) const
    {
        switch (caseindex)
        {
            case 0: //nix
                break;
            case 7: //alles
                break; //TODO area
            default:
                std::cerr << "Error: invalid case number: " << caseindex << " , this makes no sense. Eastern end of neighborhood: px number " << neighborship.at(0) << std::endl;
                throw std::invalid_argument( "minkmapSphere: Weird caseindex" );
        }
        minkTensorIntegrand theTensor(rankA, rankB, curvIndex);
        return theTensor;
    }
    
    
    pointing interpPointing(pointing A, double valA, pointing B, double valB) // Do the Mantz et al 2008 interpolation between 2 lattice points
    {
        // TODO proper interpolation

        return pointing (pi/2,0);
    }
    
    //template <typename tensortype>
    tensor2D at(int pixnum) const //TODO: poles
    {
        fix_arr<int, 8> neighbors; //neighbors of this pixel
        originalMap.neighbors(pixnum,neighbors);
        std::vector<int> easternNeighborship{pixnum, neighbors[4],neighbors[5],neighbors[6]}; //neighbors east of this pixel and this pixel
        //calculate one marching square/triangle
        tensor2D tensorHere = integrateMinktensor(easternNeighborship);
        return tensorHere;
    }
    
};


template <typename ltype, typename rtype> //For sums of minkmaps
struct minkmapSum
{
    Healpix_Map<double>& originalMap;
    uint rankA, rankB, curvIndex;
    
    const rtype* rhs;
    const ltype* lhs;
    
    minkmapSum(const ltype* left,const rtype* right) : rhs(right), lhs(left), originalMap(left->originalMap), rankA(left->rankA), rankB(left->rankB), curvIndex(left->curvIndex)
    {
        if( (right->rankA != left->rankA) || (right->rankB != left->rankB) || (right->curvIndex != left->curvIndex))
        {
            std::cerr << "Error: trying to add maps of tensors of different type. Left has (rankA, rankB, curvIndex) =  (" << left->rankA << "," << left->rankB << "," << left->curvIndex << "), Right has ("  << right->rankA << "," << right->rankB << "," << right->curvIndex << "), this makes no sense" << std::endl;
            throw std::invalid_argument( "minkmapSum: Different parameters in addition" );
        }
    }
    
    template <typename tensortype>
    tensortype at() const
    {
        return (rhs->at()) + (lhs->at());
    }
    
};


#endif
