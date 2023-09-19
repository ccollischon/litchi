/*
 * This file is part of litchi, a lightweight C++ library
 * for Minkowski analysis
 * 
 * Copyright (C) 2021-2023 Caroline Collischon <caroline.collischon@fau.de>
 * 
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */


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
