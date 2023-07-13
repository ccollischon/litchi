#ifndef minkTensorStack_hpp
#define minkTensorStack_hpp


/** \file minkTensorStack.hpp
 * \brief minkTensorStack and its operators
 */


#include "minkTensorIntegrand.hpp"

#include "healpix_cxx/pointing.h"

#include <cmath>
#include <algorithm>
#include <type_traits>
#include <vector>
#include <stdexcept>
#include <list>
#include <memory_resource>



///Save stacks of normal vectors and weights in one class
struct minkTensorStack
{
    const uint rankA{0}, rankB{0};
    const uint curvIndex{0};
    pointing r{1.5701963268,0};
    uint numnan{0}; ///< tracking how many contributions from masked pixels this tensor contains (in addition to content in nweights)
    uint numnull{0}; ///< tracking how many contributions from empty windows this tensor contains (in addition to content in nweights)
    
    using weightedN = std::pair<pointing,double>;
    
    std::pmr::monotonic_buffer_resource buffer{};
    std::pmr::polymorphic_allocator<weightedN> pa{&buffer};
    std::pmr::list<weightedN> nweights{pa}; ///< list normal vectors from which minkTensorIntegrands should be generated, and their respective weights
    
    
    minkTensorStack(minkTensorStack left, minkTensorStack right) : rankA(left.rankA), rankB(left.rankB), curvIndex(left.curvIndex), r(left.r), numnan(left.numnan), numnull(left.numnull), nweights(std::move(left.nweights),left.nweights.get_allocator())
    {
        numnan += right.numnan;
        numnull += right.numnull;
        appendStack_rr(std::move(right));
    }
    
    minkTensorStack(uint rank1, uint rank2, uint curvInd, const pointing& rNew, uint capacity=4) : rankA(rank1), rankB(rank2), curvIndex(curvInd), r(rNew),
                                                                                                   buffer{std::pmr::monotonic_buffer_resource(8*capacity*sizeof(weightedN))}, 
                                                                                                   pa{&buffer}, nweights{ pa }
    {
    }
    
    explicit minkTensorStack(const minkTensorIntegrand& inp, double weight=1) : rankA(inp.rankA), rankB(inp.rankB), curvIndex(inp.curvIndex), r(inp.r)
    {
        if(std::isnan(weight))
        {
            numnan = 1;
        }
        else if (std::abs(weight)<1e-15)
        {
            numnull = 1;
        }
        else
        {
            nweights.emplace_back(std::make_pair(inp.n,weight));
        }
    }
    
    //Move/copy constructors default
    minkTensorStack(const minkTensorStack& other) : rankA(other.rankA), rankB(other.rankB), curvIndex(other.curvIndex), r(other.r), numnan(other.numnan), numnull(other.numnull), buffer{}, pa{&buffer}, nweights(other.nweights,pa)
    {
    }
    
    
    minkTensorStack(minkTensorStack&& other) = default;
    ~minkTensorStack() = default;
    
    ///Move/copy assignments leave rank untouched (checked with assert)
    minkTensorStack& operator= (const minkTensorStack& other)
    {
        assert(rankA==other.rankA && rankB==other.rankB && curvIndex==other.curvIndex && "Trying to copy assign minkTensorStacks of different rank!");
        nweights = other.nweights;
        r = other.r;
        numnan = other.numnan;
        numnull = other.numnull;
        return *this;
    }
    
    ///Move/copy assignments leave rank untouched (checked with assert)
    minkTensorStack& operator= (minkTensorStack&& other)
    {
        assert(rankA==other.rankA && rankB==other.rankB && curvIndex==other.curvIndex && "Trying to move assign minkTensorStacks of different rank!");
        nweights = std::move(other.nweights);
        r = std::move(other.r);
        numnan = other.numnan;
        numnull = other.numnull;
        return *this;
    }
    
    /**
     * Checks whether the amont of masked/NAN contributions to this stack is too high.
     * Allows for a relative amount of masked pixels up to 1/16
     */
    bool isMasked() const
    {
        return 15*numnan > nweights.size()+numnull;
    }
    
    /**
     * Checks whether tensor contains contours/anything non-empty. Empty spots are usually set to NAN in output map
     */
    bool isEmpty() const
    {
        return nweights.size() == 0;
    }
    
    /**
     * Access element of linear combination of tensors determined by weights
     * \param indices Indices of desired element
     * \return Value of desired element
     */
    double accessElement(const std::vector<uint_fast8_t>& indices) const
    {
        if(isMasked()) {return NAN;}
        if(isEmpty()) {return 0.;}
        
        double retval = 0.;
        for(const auto& element : nweights)
        {
            minkTensorIntegrand tensorHere(rankA, rankB, curvIndex, r, std::get<0>(element));
            retval+= tensorHere.accessElement(indices)*std::get<1>(element);
        }
        return retval;
    }
    
    /**
     * Access element of linear combination of tensors and rescales with total/masked # pixels
     * \param indices Indices of desired element
     * \return Value of desired element
     */
    double accessElement_rescaled(const std::vector<uint_fast8_t>& indices) const
    {
        double factor = 1.*(nweights.size()+numnull+numnan)/(1.*(nweights.size()+numnull));
        
        double retval = 0.;
        for(const auto& element : nweights)
        {
            minkTensorIntegrand tensorHere(rankA, rankB, curvIndex, r, std::get<0>(element));
            retval+= tensorHere.accessElement(indices)*std::get<1>(element);
        }
        return retval*factor;
    }
    
    
    /** Parallel transport all normal vectors in stack to newR along geodesic
     */
    void moveTo(const pointing& newR)
    {
        for(std::pair<pointing,double>& element : nweights)
        {
            std::get<0>(element) = parallelTransport(r, newR, std::get<0>(element));
        }
        r = newR;
    }
    
    
    /** Add normal vector of minkTensorIntegrand with same ranks to stack. Normal vector is parallel transported to position of stack if positions differ
     */
    void addMinkTensorIntegrand(const minkTensorIntegrand& tens, double weight=1)
    {
        assert(rankA==tens.rankA && rankB==tens.rankB && curvIndex==tens.curvIndex && "Trying to addMinkTensorIntegrand to minkTensorStack of different rank!");
        if(std::isnan(weight))
        {
            numnan += 1;
        }
        else if (std::abs(weight)<1e-15) 
        {
            numnull += 1;
        }
        else
        {
            pointing newn = tens.n;
            if(arclength(r,tens.r)>1e-12)   newn = parallelTransport(tens.r,r,newn);
            nweights.emplace_back(std::make_pair(newn,weight));
        }
    }
    
    /** Add normal vectors and weights of other stack with same ranks to this stack. Normal vectors are parallel transported to position of this stack if positions differ. Rvalue reference version.
     */
    void appendStack_rr(minkTensorStack&& other)
    {
        assert(rankA==other.rankA && rankB==other.rankB && curvIndex==other.curvIndex && "Trying to append minkTensorStacks of different rank!");
        numnan += other.numnan;
        numnull += other.numnull;
        
        if(arclength(r,other.r)>1e-12)   other.moveTo(r);
        
        auto itEnd = nweights.end();
        nweights.insert(itEnd, other.nweights.begin(), other.nweights.end());
        
    }
    
    /** Add normal vectors and weights of other stack with same ranks to this stack. Normal vectors are parallel transported to position of this stack if positions differ. Copy version.
     */
    void appendStack(minkTensorStack other)
    {
        appendStack_rr(std::move(other));
        
    }
    
    
    minkTensorStack& operator+= (const minkTensorStack& other)
    {
        appendStack(other);
        return *this;
    }
    
    minkTensorStack& operator+= (minkTensorStack&& other) //rvalue implementation because it can be used by normalHealpixInterface::toHealpix
    {
        appendStack_rr(std::move(other));
        return *this;
    }
    
    minkTensorStack& operator+= (const minkTensorIntegrand& other)
    {
        addMinkTensorIntegrand(other);
        return *this;
    }
    
    minkTensorStack& operator*= (double other)
    {
        if(std::isnan(other))
        {
            numnan += nweights.size();
            nweights.clear();
        }
        else if(std::abs(other)<1e-15)
        {
            numnull += nweights.size();
            nweights.clear();
        }
        else
        {
            std::for_each(nweights.begin(), nweights.end(), [&other](auto& inp){std::get<1>(inp)*=other;});
        }
        return *this;
    }
};

minkTensorStack operator+ (minkTensorStack lhs, const minkTensorStack& rhs)
{
    minkTensorStack returnval (std::move(lhs));
    returnval += rhs;
    return returnval;
}

minkTensorStack operator+ (const minkTensorIntegrand& lhs, const minkTensorStack& rhs)
{
    minkTensorStack returnval (rhs);
    returnval.addMinkTensorIntegrand(lhs);
    return returnval;
}


///minkTensor concept for either minkTensorStack or minkTensorIntegrand
template<typename T>
concept minkTensor = std::is_base_of_v<minkTensorStack,std::remove_cv_t<T>> || std::is_base_of_v<minkTensorIntegrand,std::remove_cv_t<T>>;

template<minkTensor left>
minkTensorStack operator+ (const left& lhs, const minkTensorIntegrand& rhs)
{
    minkTensorStack returnval (lhs);
    returnval.addMinkTensorIntegrand(rhs);
    return returnval;
}

template<minkTensor left>
minkTensorStack operator* (const left& lhs, double rhs)
{
    minkTensorStack retStack(lhs);
    retStack*=rhs;
    return retStack;
}

template<minkTensor right>
minkTensorStack operator* (double lhs, const right& rhs)
{
    minkTensorStack retStack(rhs);
    retStack*=lhs;
    return retStack;
}


///Empty minkTensorStack that represents a spot set to NAN
minkTensorStack nanTensor(uint rank1, uint rank2, uint curvInd, const pointing& rNew)
{
    minkTensorStack thistensor(rank1,rank2,curvInd,rNew);
    thistensor.numnan=1;
    return thistensor;
}

///Empty minkTensorStack that represents a spot without structure
minkTensorStack nullTensor(uint rank1, uint rank2, uint curvInd, const pointing& rNew)
{
    minkTensorStack thistensor(rank1,rank2,curvInd,rNew);
    thistensor.numnull=1;
    return thistensor;
}




#endif
