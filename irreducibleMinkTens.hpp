#ifndef irreducible
#define irreducible

#include "healpix_cxx/pointing.h"

#include <complex>
#include <cmath>

struct irreducibleMinkTens {
    const int l_, m_;
    pointing r_{0,0};
    pointing n_{0,0};
    
    irreducibleMinkTens(int l, int m, pointing r, pointing n) : l_{l}, m_{m}, r_{std::move(r)}, n_{std::move(n)}
    {
    }
    
    //~ std::complex<double> accessElement() const // Y*lm(n)
    //~ {
        //~ double thetapart = std::sph_legendre(l_,m_,n_.theta);
        //~ return thetapart*std::exp(std::complex<double>(0.0, -1.0) * (double)m_ * n_.phi);
    //~ }
    
    std::complex<double> accessElement() const // exp(i l phi(n))
    {
        double phi = giveAngle(n_,r_);
        return std::exp(std::complex<double>(0.0, +1.0) * (double)l_ * phi);
    }
    
    
};



#endif
